
%test the power spectrum against camb

[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');
cosm=get_default_cosm;

[kc,psc]=textread('/home/greer//src/camb/test_matterpower.dat','%f %f');

%output from camb is in k/h and h/Mpc units
%remove the h dependence
kc = kc*cosm.h;
psc = psc/(cosm.h^3);

%get delta m = k^3*P(k)/2/pi/pi
dm = Delta_m(kc,1.0,cosm,opt);

%convert to power spectrum 
psg = 2*pi*pi*dm./(kc.^3);

semilogx(kc,(psc-psg)./psc);
%loglog(kc,psc,kc,psg);
%legend('CAMB','Greer');
xlabel('k [Mpc^{-1}]');
ylabel('Relative Difference');

%%%%%%%%%%%%%%%%

%test mass function code against tinker and eduardo

[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');
cosm=get_default_cosm;
[mg,mfg]=compute_halo_mf(0.25,cosm);
[mr,mfr]=textread('greer/mf/rozo/mf_z=0.0.dat','%f %f');
[mt,mft,rt,smt,fst,dsdmt]=textread('greer/mf/tinker/test.dndM','%f %f %f %f %f %f');

%my numbers are in 10^14 Msun:
mg = mg*1e14;
mfg = mfg/1e14;

%tinker is in h^-1 MPc etc.
mt=mt/cosm.h;
mft = mft * (cosm.h)^4;

%spline eduardo and tinker onto my mass values
mfrs=spline(mr,mfr,mg);
mfts=spline(mt,mft,mg);

%plot

loglog(mg,mfg,mg,mfrs,mg,mfts);
title('raw dn/dM for various methods');
legend('Greer z=0.25','Rozo z=0','Tinker z=0.25');
xlabel('Mass (Msol)')
ylabel('dn/dM')


%semilogx(mg,(mfrs-mfts)./mfrs);
semilogx(mg,(mfts-mfg)*100./mfts);
ylabel('Percent deviation.');
xlabel('Mass.');
title('Compare Tinker - Greer at z=0.25');


%%%%%%%%%%%%%%%%%%%%%%%

%compute sigma(M)

clear
cosm=get_default_cosm;
a=0.8; nstep=2001; mmin=0.05/cosm.h; mmax=100/cosm.h;
dlogm = (log(mmax) - log(mmin))/(nstep-1);
[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');

mg = zeros(nstep,1);
smg=0*mg; fsg=0*mg; dsdmg=0*mg; mfg=0*mg;

for i=1:nstep
  mg(i) = exp((i-1)*dlogm)*mmin;
  smg(i) = sigma_M(mg(i),a,cosm,opt);
  fsg(i) = fsigma(mg(i),a,cosm,opt);
  dsdmg(i) = dsigma2_dM(mg(i),a,cosm,opt);
  mfg(i) = mass_function(mg(i),a,cosm,opt);
end

[mt,mft,rt,smt,fst,dsdmt] = textread('greer/mf/tinker/test.dndM',...
	'%f %f %f %f %f %f');

dsdmg = dsdmg./smg/2;

%mass in 10^14
mg=mg*1e14; dsdmg = dsdmg/1e14; mfg=mfg/1e14;

%tinkers in h^-1 Msun, etc.
mt = mt/cosm.h; rt = rt/cosm.h; mft = mft * cosm.h^4; dsdmt = dsdmt*cosm.h;

subplot(2,2,1)
semilogx(mg,(smg-smt)./smg); axis tight;
xlabel('M'); ylabel('Relative Difference'); title('\sigma(M)');
subplot(2,2,2)
semilogx(mg,(fsg-fst)./fsg); 
xlabel('M'); ylabel('Relative Difference'); title('f(\sigma)');
subplot(2,2,3)
semilogx(mg,(dsdmg-dsdmt)./dsdmg); axis tight;
xlabel('M'); ylabel('Relative Difference'); title('d\sigma/dM');
subplot(2,2,4)
semilogx(mg,(mfg-mft)./mfg); axis tight;
xlabel('M'); ylabel('Relative Difference'); title('dn/dM');

subplot(2,2,1)
loglog(mg,smg,mt,smt);
ylabel('\sigma')
subplot(2,2,2)
loglog(mg,fsg,mt,fst)
ylabel('f(\sigma)')
subplot(2,2,3)
loglog(mg,dsdmg,mt,dsdmt)
ylabel('d\sigma/dM')
subplot(2,2,4)
loglog(mg,mfg,mt,mft);

%%%%%%%%%%%%%%%%%%

[obs,par,opt] = read_mcmc_conf('greer/mf/test.in');
opt.MMIN=0.1;
cosm=get_default_cosm;
[mg,mfg]=compute_halo_mf(0.25,cosm,opt);
area = 0.6; zmin=0.2; zmax=0.3;
vol = (area/4.0/pi) * co_vol([zmin,zmax],cosm,opt,'Mpc'); %in Mpc/h
%vol = vol/(cosm.h)^4;

mg=mg;
mfg = mfg*vol;

loglog(mg,mfg)
