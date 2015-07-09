%cov with gaussian priors
function plot_maxbcg_mcmc(in_fn,out_fn,mass_pivot)

[data,labels]=read_markov_output('greer/mf/like_test7.out');
[data] = unwrap_markov_data(data);
[data,labels]=plot_mcmc_cov(data,labels,50,0);
[obs,par,opt]=read_mcmc_conf('greer/mf/test1.in');
n=101;
for i=find(~isnan(par.prior_cent_val))
  x0 = par.prior_cent_val(i);
  sig2 = (x0*par.prior_widths(i))^2;
  subplot(7,7,1+(i-1)*8);
  hold on;
  xl = xlim; 
  h = (xl(2)-xl(1))/(n-1);
  x = xl(1):h:xl(2);
  gauss = exp( -0.5 * ((x-x0).^2)./sig2 )./sqrt(2*pi*sig2);
  hl=line([x0 x0], [0 1]);
  set(hl,'Color','r');
  plot(x,gauss/max(gauss),'r');
  hold off;
end

for i=1:7;
 x0 = par.init(i);
 subplot(7,7,1+(i-1)*8);
 hold on;
 hl=line([x0 x0], [0 1]);
 set(hl,'Color','g');
 hold off;
end

%%%%%%%%%%%%%%%%%%%


%find covariance matrix for parameters

[data,labels]= read_markov_out('greer/mf/like_test3.out',50,1);    
data_uw=unwrap_markov_data(data);
par_cov = cov(data_uw);

[V,D]=eig(par_cov);
V


%%%%%%%%%%%%%%%


%play with lognormal dists

x=0.01:0.01:10;
mu=5;
sigma=1;
y = exp(-0.5 * (x-mu).^2/2/sigma/sigma)/sqrt(2*pi)/sigma;

plot(log(x),y)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc num expected above certain lambda:

clear;
[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');
obs = gen_fake_obs(obs,par,opt);

opt.MMIN=0.100;
opt.LMIN=16.2628;
opt.LMAX=1500;
n = 901;
zobs=0.25; aobs = 1./(1+zobs);

vol = (obs.area/4/pi).*co_vol([obs.zmin,obs.zmax],par.cosm,opt,'Mpc');

alpha = par.init(4); beta=par.init(5); sigma=par.init(6);
cosm = par.cosm; 

lognmin = log(opt.LMIN/60);
lognmax = log(opt.LMAX/60);
h = (lognmax-lognmin)/(n-1);
logn = (lognmin:h:lognmax);
wx = ones(1,n); wx(2:2:n-1) = 4; wx(3:2:n-2) = 2; wx=wx*h/3; %weights

dndlam = calc_dndobs(logn,aobs,alpha,beta,sigma,par.cosm,opt);
ncl_i = exp(logn).*dndlam;
ncl_expect = vol.* wx * ncl_i'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot P(lam/mass)

[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');

opt.MMIN=0.1;
opt.LMIN=5;
opt.LMAX=1500;
nstep = 201;
zobs=0.25;
alpha = par.init(4); beta=par.init(5); sigma=par.init(6);
cosm = get_default_cosm('lcdm',opt);

mmin=log(opt.MMIN);
mmax=log(opt.MMAX);
nstep=201;
h = (mmax-mmin)/(nstep-1);
logm = mmin:h:mmax;

opt.LMIN=5;
opt.LMAX=1500;

lognmin = log(opt.LMIN/60);
lognmax = log(opt.LMAX/60);
h = (lognmax-lognmin)/(nstep-1);
logn = (lognmin:h:lognmax);

[logm2d,logn2d] = meshgrid(logm,logn); %logm2d=logm2d(:); logn2d=logn2d(:);

logmobs = (logn2d-beta)/alpha;

x=logmobs-logm2d;
p_lam = exp(-x.*x/2/sigma/sigma)/sqrt(2*pi*sigma*sigma);


p_lam_1=p_lam(100,:);
plot(logn,p_lam_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate dn/dlambda
%3/20/12

[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');

opt.MMIN=0.01;
opt.MMAX=500;
opt.NNODES=299;
opt.LMIN=10;
opt.LMAX=1000;
nstep = 401;
zobs=0.25;
alpha = par.init(4); beta=par.init(5); sigma=par.init(6);
cosm = get_default_cosm('lcdm',opt);

%work in units of lambda/60
lognmin = log(opt.LMIN/60);
lognmax = log(opt.LMAX/60);
h = (lognmax-lognmin)/(nstep-1);
logn = (lognmin:h:lognmax).';

[x1,x2] = calc_dndobs(logn,zobs,alpha,beta,sigma,cosm,opt);
figure(4)
loglog(exp(logn)*60,x1,exp(logn)*60,x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m mf] = compute_halo_mf(0.0,[]);


par_init = [1.78, log(2.925e-5*0.6145) , 0.2, 1/1.06, -0.95/1.06, 0.45, 0.191];
par_var  = [0.07,  0.07, 0.04,   0.08,   0.08, 0.04, 0.04];


par_y.alpha = par_init(1);
par_y.beta  = par_init(2);
par_y.sigma = par_init(3);

par_n.alpha = par_init(4);
par_n.beta  = par_init(5);
par_n.sigma = par_init(6);
r = par_init(7);

%generate test data
C = [par_y.sigma^2 r*par_y.sigma*par_n.sigma; ...
        r*par_y.sigma*par_n.sigma par_n.sigma^2];
[V,D] = eigs(C);
U = (V*sqrt(D))';
[y N] = gen_obs(25, par_y, par_n, r, m, mf,U);
[par_obs.c par_obs.d par_obs.Sy] = calc_yN_fit(y,N);

%----------------------

[m mf] = compute_halo_mf(0.0,[]);
cmf=cumsum(mf)/sum(mf);

m_lim=1;
[m1 cmf1] = cum_mass_fcn(m_lim,0.0,m,mf);

loglog(m,cmf,m1,cmf1);
axis([1 10 0.01 1.1])


%-----------------------

[m mf] = compute_halo_mf(0.0,[]);

frac_count_mf = cumsum(mf)/sum(mf);

lcl = get_cl_samp3(10000,1,m,mf);
cl  = 1e2*exp(lcl);
binned_cl = hist(cl,m);
frac_count_cl = cumsum(binned_cl)/sum(binned_cl);

loglog(m,frac_count_mf,'b',m,frac_count_cl,'g.');
axis([0.9 1e2 0.01 1.1]);
xlabel('Mass (10^{12} M_{sun})');
ylabel('Fractional Count');

%-------------------
[m mf] = compute_halo_mf(0.0,[]);

numpts = 2000;
minm = log(1);
maxm = log(1e4);
delta=(maxm-minm)/(1.0*(numpts-1.0));
%m = minm:delta:maxm;
dm = diff( exp(minm:delta:(maxm+delta)) );

cdf = cumsum(mf)/sum(mf);
%pdf = diff(cdf)./diff(m);
pdf = mf./dm./sum(mf);

m_pdf = m(2:length(m));
pdf1=diff(cdf)./diff(m);

loglog(m,cdf,m,pdf,m_pdf,pdf1);


%----------------------------

[m mf] = compute_halo_mf(0.0,[]);

m_lim = 1;
ncl = 10000;

ii = find(m >= m_lim);
m_in = m(ii);
mf_in = mf(ii);
m_in = log(m_in/1e2);

%m_out = randsample(m_in,ncl,true,mf_in);
p = mf_in(:)' / sum(mf_in);
edges = min([0 cumsum(p)],1);
edges(end) = 1;
[~, y] = histc(rand(ncl,1),edges);
y = m_in(y);


%------test mf sample draws like in mcmc_best-----
[m mf] = compute_halo_mf(0.0,[]); %m is in units of 1e12 Msun.
m_lim = 5e1; %only consider masses above 5e13 Msun.
ii = find(m >= m_lim); %don't carry around stuff we don't need
m = m(ii);
mf = mf(ii);
m = log(m/1e2); %lets work in log 1e14 Msun units cause everyone else does
frac_count_mf = cumsum(mf)/sum(mf);


%set up probability for drawing from mf
%hacked from randsample
p = mf(:)'/sum(mf);
mf = min([0 cumsum(p)],1);
mf(end) = 1;
rands = rand(5000,1);
[~, y] = histc(rands,mf);
lcl = m(y);

binned_cl = hist(lcl,m);
frac_count_cl = cumsum(binned_cl)/sum(binned_cl);

loglog(exp(m),frac_count_mf,'b',exp(m),frac_count_cl,'g.');
axis([0.5 10 0.1 1.1]);
xlabel('Mass (10^{14} M_{sun})');
ylabel('Fractional Count');


