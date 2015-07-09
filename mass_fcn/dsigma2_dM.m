function dsdm = dsigma2_dM(m,a,cosm,opt)
%
%calculate the derivative of the var of the
%matter density field

%rc = rho_crit(a,cosm,'tinker','expan');
rc = rho_crit(1.0,cosm,'tinker','expan');

%temp=(3.0*m)/(4.0*pi*cosm.Omega_m*rc);
%temp=(temp^(1.0/3.0));

r = ( (3.0.*m)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);


%switch opt.INTEGRATOR
%   case 'quad'
%      dsdm = quad(@(x)dsigma_dR_integrand(x,a,r,cosm,opt),-0.95,0.95,1e-6);
%   case 'composite'
      n = opt.NNODES; %number of nodes for integration
      h = (0.95+0.95)/(n-1); %spacing for nodes
      x = (-0.95:h:0.95).'; %node locations
      wt = ones(1,n); wt(2:2:n-1) = 4; wt(3:2:n-2)= 2; wt=wt*h/3; %weights
      dsdmi = wt * dsigma_dR_integrand(x,a,r,cosm,opt); %integral
%   otherwise
%     error('mass_function.m: function dsdm -- integration method not specified');
%end

%mlo=0.99*m;
%mhi=1.01*m;

%rlo=( (3.0.*mlo)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);
%rhi=( (3.0.*mhi)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);

%slo=sigma_R(rlo,a,cosm,opt);
%shi=sigma_R(rhi,a,cosm,opt);

%dsdm = (shi-slo)./(mhi-mlo);

%sm = sigma_R(r,a,cosm,opt);
dsdm = (1.0/3.0).*r.*dsdmi./m;
%dsdm = dsdm .* r.^(11/3)./m./sm/6;
%dsdm =dsdm.*(1.0/6.0).*(1./m./sm).*r^(11/3);

return


function dsdRi = dsigma_dR_integrand(x,a,R,cosm,opt)

%k=exp( x./(1.0-x.*x) );
%dm=Delta_m(k,a,cosm,opt);

%kR=k.*R;

%w = tophat_window(kR);

%dwdr = 3.*k.*((kR.^2-3).*sin(kR)+3.*kR.*cos(kR))./(kR.^4);
%dw2dr = 2.*w.*dwdr;

%dsdRi = dm .* dw2dr .* (1.0+x.*x)./(1.0-x.*x)./(1.0-x.*x);



k=exp( x./(1.0-x.*x) );
val=Delta_m(k,a,cosm,opt);

kR=k.*R;
wf=tophat_window(kR);
dwf=(3.0./kR).*( sin(kR)./kR - wf );% /* [window function] */
val=val.*2.0.*wf.*dwf.*k;

val=val.*(1.0+x.*x)./(1.0-x.*x)./(1.0-x.*x); %/* Measure for integration variable. */
dsdRi = val;

return;
