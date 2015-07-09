function obs = gen_fake_obs(obs,par,opt)
%outputs ncl_expect
%zobs
%aobs
%vol
%logn
%logy
%logy_err
%logm

obs.zobs = mean([obs.zmin,obs.zmax]).*ones(obs.ncl,1);
obs.aobs = 1./(1+obs.zobs);

%find mass bins
n = opt.NUMPTS;
M = logspace(log10(opt.MMIN),log10(opt.MMAX),n-1);
dM = diff(M);
M_cent = M(1:end-1)+dM/2;
%M_scal = 0.*M_cent;

dn_dM = 0.*M_cent;
for i=1:length(M_cent)
  dn_dM(i) = mass_function(M_cent(i),median(obs.aobs),obs.simcosm,opt);
  %M_crit_i = conv_Mdelta(M_cent(i),obs.zobs(i),500,500,'mean','crit',cosm);
  %M_scal(i) = M_cent(i)/M_crit_i;
end

%convert from 500m to 500c definition
%dn_dM./M_scal;
%M = M./M_scal;
%dM = diff(M);

%calculate expected value of clusters per bin
dN_dM = obs.vol.*dn_dM;
N = dM .* dN_dM;

%poisson deviate around mean expect N
N_obs = random('Poisson',N,1,length(N));

lcl = [];

%get all the cluster masses needed
for i=1:length(N_obs)
  lcl = [lcl; ones(N_obs(i),1)*M_cent(i)];
end

%convert to mass pivot point
lcl = lcl/obs.simpivot;

%work in logm
lcl = log(lcl);

C = [par.nSigma^2 par.rhoYN*par.ySigma*par.nSigma; ...
   par.rhoYN*par.nSigma*par.ySigma par.ySigma^2];

%get mean Y (1e-5 Mpc^2)
logy_true = par.yAlpha * lcl + par.yBeta;

%get mean N (N/60)
logn = par.nAlpha * lcl + par.nBeta;

%get scatter using covariance matrix
scatters = corr_rand(length(lcl),C);
logy_true = logy_true + scatters(:,2);
logn = logn + scatters(:,1);

logy_err = obs.sim_yerr*ones(size(logy_true));
%logy_err = obs.sim_yerr*randn(size(logy_true)); %some random measurement error on Y

%logy = logy_true+logy_err;
logy = logy_true;

%pick the richest clusters
[~, ii] = maxk(logn,obs.ncl,1);

obs.logn = logn(ii); %logn is log(N/60) around simpivot
obs.logy = logy(ii); %logy is log(Y/1e-5 Mpc^2) around simpivot
obs.logy_err = abs(logy_err(ii));
obs.logm = lcl(ii); %logm is log(simpivot x 1e14 Msun)`

if(obs.missing_y)
	obs.logy(obs.missing_y) = nan;
  obs.logy_err(obs.missing_y) = nan;
%	obs.logy_err(obs.missing_y) = 100*obs.logy_err(obs.missing_y);
end


ind = strfind(opt.fn,'.');
ind = max(ind);
fn_rt=opt.fn(1:ind-1);
sim_fn= strcat(fn_rt,'.sim');
fid=fopen(sim_fn,'w');
fprintf(fid,'#lambda  ysz  ysz_err/ysz  mass\n');
for i=1:obs.ncl
  fprintf(fid,'%.2f %.6f %.2f %.3f \n',exp(obs.logn(i))*60,...
		exp(obs.logy(i))*1e-5,obs.logy_err(i),exp(obs.logm(i))*par.mass_pivot);
end
fclose(fid);

return
	
