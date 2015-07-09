function log_like = calc_like(obs,par,opt)

   method = opt.method;
   pr = calc_prior(par);
   log_like = log(pr);
%   fprintf('Log(L) prior: %f\t',log_like);
   if(pr~=0)
      switch lower(method)
         %case 'fast'
         %   log_like = log_like+calc_fast_like(obs,par,opt);
         %case 'best'
         %   log_like = log_like+calc_best_like(10000,obs,par);
         case 'full'
            log_like = log_like+calc_full_like(obs,par,opt);
         case 'prior'
      log_like = log_like;
         case 'full_noprior'
            log_like = calc_full_like(obs,par,opt);
         otherwise
            error('calc_like.m: method not recognized.');
      end
   end
%   fprintf('Log(L): %f\n',log_like);
return

function prior = calc_prior(par)
   %implement various priors.

   %prior #1
   %gaussian priors around a central value

   %convert normalizations to prior pivot
   yBeta = par.yBeta+par.yAlpha*log(par.prior_pivot(1)/par.mass_pivot);
   nBeta = par.nBeta+par.nAlpha*log(par.prior_pivot(2)/par.mass_pivot);
   val = [par.yAlpha, yBeta, par.ySigma, ...
          par.nAlpha, nBeta, par.nSigma, par.rhoYN];

   %val = [par.yAlpha, par.yBeta, par.ySigma, ...
   %      par.nAlpha, par.nBeta, par.nSigma, par.rhoYN];

   %if prior value is nan there is no prior
   ii = ~isnan(par.prior_cent_val);
   val = val(ii);
   prior_cent_val = par.prior_cent_val(ii);
   prior = normal(val-prior_cent_val,par.prior_cov_mat);

   %prior 2
   %require that sigmas be >0 and rho must be between -1 and 1
   if( (par.ySigma < 3e-2) | (par.nSigma < 3e-2) )
      prior = 0;
   elseif( (par.rhoYN > 1) | (par.rhoYN < -1) )
      prior = 0;
   end

   %prior 3
   %check and see if there is a flat prior on the variance of the scatter
   %if so, multiply prior by the value. gaussian prior on scatter overrides
   %this
   if( (par.flat_var_priors(1) ~= 0 ) & (isnan(par.prior_cent_val(3))) )
      prior = prior * par.ySigma;
   end

   if( (par.flat_var_priors(2) ~= 0) & (isnan(par.prior_cent_val(6))) )
      prior = prior * par.nSigma;
   end

return

function log_like = calc_full_like(obs,par,opt)
  %returns log(L) for a given set of observations and parameters
  %opt is a structure with various options

  %integration limits
  mmin = log(opt.MMIN);
  mmax = log(opt.MMAX);

  lam_min = min(obs.logn);
  lam_max = log(opt.LMAX/60);

  ymin = log(opt.YMIN);
  ymax = log(opt.YMAX);

  %some cosmological values used by all integrals
  aobs = median(1 ./ (obs.zobs+1));
  vol = (obs.area/4/pi).*co_vol([obs.zmin,obs.zmax],par.cosm,opt,'Mpc');

  %set up weights and nodes for numerical integration
  n = opt.NNODES;
  h = (lam_max - lam_min)/(n-1);
  loglam = (lam_min:h:lam_max).';
  wl = ones(1,n); wl(2:2:n-1) = 4; wl(3:2:n-2) = 2; wl=wl*h/3; %weights
  h = (mmax - mmin)/(n-1);
  logm = (mmin:h:mmax);
  wm = ones(1,n); wm(2:2:n-1) = 4; wm(3:2:n-2) = 2; wm=wm*h/3; %weights
  h = (ymax-ymin)/(n-1);
  logy = (ymin:h:ymax).';
  wy = ones(1,n); wy(2:2:n-1) = 4; wy(3:2:n-2) = 2; wy=wy*h/3;

  [logm2d,logy2d] = meshgrid(logm,logy); logm2d=logm2d(:); logy2d=logy2d(:);
  wmy2d = wy.'*wm; wmy2d=wmy2d(:).';

  [logm2d,loglam2d] = meshgrid(logm,loglam);
  logm2d=logm2d(:); loglam2d=loglam2d(:);
  wmlam2d = wl.'*wm; wmlam2d=wmlam2d(:).';

  %construct or load the mass function. need to change this to always
  %construct if cosmology changes
  %mass function is in units of 1e14 Msun
  m = exp(logm); mf = 0.*m;
  m2d = exp(logm2d); mf2d = 0.*logm2d;
  y2d = exp(logy2d);
  if(exist(opt.mf_fn,'file'))
    mf = load(opt.mf_fn,'-ascii','mf');
  else
    fprintf('Mass function file doesnt exist. Creating it...\n');
    for j=1:length(m)
      mf(j) = mass_function(m(j),aobs,par.cosm,opt);
    end
    mf = mf(:); mf2d=mf2d(:);
    save(opt.mf_fn,'-ascii','mf');
  end
  for j=1:length(m)
    mf2d(m2d==m(j))=mf(j);
  end

  %calculate the number of expected clusters at a range of loglam
  ncl_expect = vol * wmlam2d * poisson_int2(logm2d,loglam2d,par,mf2d);

  %fprintf('ncl_expect: %f\n',ncl_expect);
  %expect 0 for the likelihood
  poi_log_like = -ncl_expect;

  %do marginalization over the mass function evaluated at observed Y and lambda
  cl_integral = zeros(obs.ncl,1);
  for i=1:obs.ncl
    if ( any(obs.missing_y == i) )
       %do dn\dlam integral, not dn\dlam dY
       cl_integral(i) = vol*wm*cl_int_noY(logm,m,obs.logn(i),par,mf);
    else
       cl_integral(i) = vol*wmy2d*cl_int2(logm2d,logy2d,m2d,y2d,...
          obs.logn(i),obs.logy(i),obs.logy_err(i),par,mf2d);
    end
  end

  %working in log space, so take the log and sum them up
  cl_log_like = sum(log(cl_integral));

  %return log(L)
  log_like = poi_log_like+cl_log_like;
return

function pois_int = poisson_int2(logm,loglam,par,mf)
  %integrate over the mass function weighting by P(ln lam| ln M)
  %to get dn/dlnlam. Then integrate over ln lam from the minimum
  %lambda to ~infty to get expected number of cls

  %calculate P(ln lam|ln M) -- a gaussia in ln lam
  logn_m = par.nAlpha*(logm-log(par.mass_pivot))+par.nBeta;
  sig2 = par.nSigma.^2;
  logx = loglam-logn_m;
  p_lam = exp(-0.5*logx.*logx/sig2)/sqrt(2*pi*sig2);

  %the extra exp(logm) is because mf is dn/dM, not dn/dlnM
  pois_int = p_lam.*mf.*exp(logm);
return

function clus_int = cl_int_noY(logm,m,logn_i,par,mf)
   n_i = exp(logn_i);
   logn_m = par.nAlpha.*(logm-log(par.mass_pivot))+par.nBeta;
   q_ny = -(0.5/par.nSigma/par.nSigma)*(logn_i - logn_m).^2;


   %integrand is dn/dlogm*P(log lam| log M)
   clus_int = exp(q_ny) .* m .* mf' ./ sqrt(2*pi) ./ n_i ./ par.nSigma;
   clus_int = clus_int';
return

function clus_int = cl_int2(logm,logy,m,y,logn_i,logy_i,logy_i_err,par,mf)
  %integrate the mass function over ln M weighting by the lognormal
  %P(ln lam, ln Y|ln M) to find dN/dlamdY. Then marginalize over the
  %gaussian distribution P(Yobs|Ytrue).
  %
  %pass exp(logy) and y to save the computation time

  y_i = exp(logy_i);
  n_i = exp(logn_i);
  y_i_err2 = (y_i*logy_i_err).^2;

  %here are the observables given a mass
  logn_m = par.nAlpha.*(logm-log(par.mass_pivot))+par.nBeta;
  logy_m = par.yAlpha.*(logm-log(par.mass_pivot))+par.yBeta;

  logx = [ (logn_i - logn_m), (logy-logy_m) ];

%  C = [par.nSigma.^2, par.rhoYN*par.nSigma*par.ySigma; ...
%       par.rhoYN*par.nSigma*par.ySigma, par.ySigma.^2];
%  invC = inv(C);
%  detC = det(C);

  %have unrolled the gaussian computations for speed
  t = 1-(par.rhoYN^2);

  %the exponent in the P(Yobs|Ytrue)
  q_y = -0.5*((y_i-y).^2)/y_i_err2;
  denom_y = sqrt(2*pi*y_i_err2);

  %the exponent in the P(ln Y, ln lam|ln M)
  q_ny = -(0.5/t)*( (logx(:,1)/par.nSigma ).^2 - ...
    2*par.rhoYN*logx(:,1).*logx(:,2)/par.nSigma/par.ySigma + ...
    (logx(:,2)/par.ySigma).^2 );

  %integrand is dn/dlogYdloglam*P(Yobs|Ytrue)*P(ln Y, ln lam|ln M)
  %because P(Yobs|Ytrue) is a Gaussian and we're integrating over
  %log space, we get an extra factor of Y in the numerator. This
  %cancels with one substituting the 1/Y from converting from
  %dn/dlamdY to dn/dloglamdlogY The factor of 1/n_i comes from
  %the same conversion
  clus_int = exp(q_y+q_ny) .* mf .* m ./ ...
  (2*pi*n_i*par.ySigma*par.nSigma*sqrt(t)*denom_y);
return


%%%%%deprecated functions%%%%%%%%%

function pos_int = poisson_integrand(loglam,aobs,par,opt,mf)
   %poisson integrand is mf.*vol.*P(lam|mass);

   lam = exp(loglam);
   dndlam = calc_dndobs(loglam,aobs,par.nAlpha,par.nBeta,par.nSigma,...
        par.cosm,opt,mf);
   pos_int = lam(:).*dndlam(:);
return


function clus_int = cl_integrand(logm,logn_i,logy_i,logy_i_err,aobs,par,opt,mf)
   %for each cluster integrand, it is mf .* vol .* P(lam,Y|mass);
   %also need to integrate over error in y

   mass = exp(logm);
   y_i = exp(logy_i);
   y_i_err2 = (y_i*logy_i_err).^2;
   p_lam_ysz = 0.*logm;

   %calculate P(ln lam,ln Yobs|ln mass)
   %involves an integral over Y
   ymin = log(opt.YMIN); ymax = log(opt.YMAX);
   n = opt.NNODES;
   h = (ymax-ymin)/(n-1);
   logy = (ymin:h:ymax).';
   ysz = exp(logy);
   wy = ones(1,n); wy(2:2:n-1) = 4; wy(3:2:n-2) = 2; wy=wy*h/3;

   %P(Yobs|Ytrue)
   %p_yobs_ytrue = ysz.*exp(-0.5*((y_i-ysz).^2)/y_i_err2)/sqrt(2*pi*y_i_err2);
   p_yobs_ytrue = exp(-0.5*((y_i-ysz).^2)/y_i_err2)/sqrt(2*pi*y_i_err2);

   logn_m = par.nAlpha.*logm'+par.nBeta;
   logy_m = par.yAlpha.*logm'+par.yBeta;
   C = [par.nSigma.^2, par.rhoYN*par.nSigma*par.ySigma; ...
    par.rhoYN*par.nSigma*par.ySigma, par.ySigma.^2];
   invC = inv(C);

   p_lam_ysz_true = 0.*logy;
   logn = logn_i*ones(size(logy));

   %generate a set of nx2xn arrays where the (:,:,i) entries are
   % (logn-logn_m(i) logy-logy_m(i))
   logx = repmat([logn logy],[1 1 length(logm)]) - ...
  permute(repmat([logn_m logy_m],[1,1,length(logm)]),[3,2,1]);

   %foreach nx2 array in the 3rd dim of logx, calculate p_lam_ysz
   %    p_lam_ysz is calculated by taking each 1x2 vector x in the 1st
   %  dim and multiplying it by x*invC*x'.
   p_lam_ysz_true = (logx(:,1,:).^2.*invC(1,1))+...
      (2*invC(1,2)*logx(:,1,:).*logx(:,2,:))+...
      (logx(:,2,:).^2.*invC(2,2));
   p_lam_ysz_true = exp(-0.5*p_lam_ysz_true)./(2*pi*sqrt(det(C)))./exp(logn_i);
   p_lam_ysz_int = repmat(p_yobs_ytrue,[1 1 length(logm)]).*p_lam_ysz_true;
   wy2 = repmat(wy',[1,length(logm)]);
   p_lam_ysz = sum(wy2.*squeeze(p_lam_ysz_int));
   p_lam_ysz = p_lam_ysz';

   clus_int = mass' .* mf .* p_lam_ysz;
return
