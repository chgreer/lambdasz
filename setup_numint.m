function opt = setup_numint(z,opt,par);
a = 1./(z+1);
%setup integration weights
mmin = log(opt.MMIN);
mmax = log(opt.MMAX);
ymin = log(opt.YMIN);
ymax = log(opt.YMAX);
lam_min = log(opt.LMIN/60);
lam_max = log(opt.LMAX/60);

n    = opt.NNODES;

hm = (mmax - mmin)/(n-1);
log_m = (mmin:hm:mmax).';
wm = ones(1,n); wm(2:2:n-1) = 4; wm(3:2:n-2) = 2; wm=wm*hm/3; %weights

hl = (lam_max - lam_min)/(n-1);
log_l = (lam_min:hl:lam_max).';
wl = ones(1,n); wl(2:2:n-1) = 4; wl(3:2:n-2) = 2; wl=wl*hl/3;

hy = (ymax - ymin)/(n-1);
log_y = (ymin:hy:ymax).';
wy = ones(1,n); wy(2:2:n-1) = 4; wy(3:2:n-2) = 2; wy=wy*hy/3;

[log_mm,log_yy] = meshgrid(log_m,log_y); 
log_mm = log_mm(:); log_yy = log_yy(:);
wmy = wy.'*wm; wmy = wmy(:).';

opt.log_m = log_m;
opt.wm = wm;
opt.log_y = log_y;
opt.wy = wy;
opt.log_m2dy = log_mm;
opt.log_y2dm = log_yy;
opt.wmy = wmy;

[log_mm,log_ll] = meshgrid(log_m,log_l); 
log_mm = log_mm(:); log_ll = log_ll(:);
wml = wl.'*wm; wml = wml(:).';

opt.log_l = log_l;
opt.wl = wl;
opt.log_l2dm = log_ll;
opt.log_m2dl = log_mm;
opt.wml = wml;

opt.mass = exp(log_m);
opt.mass2d = exp(opt.log_m2dy);
opt.mf = 0.*opt.mass;
opt.mf2d = 0.*opt.mass2d;

if(exist(opt.mf_fn,'file'))
  opt.mf = load(opt.mf_fn,'-ascii','mf');
else
  for ii=1:length(opt.mass)
    opt.mf(ii) = mass_function(opt.mass(ii),a,par.cosm,opt);
  end
  opt.mf = opt.mf(:);
  mf = opt.mf;
  save(opt.mf_fn,'-ascii','mf');
end
for ii=1:length(opt.mass)
  opt.mf2d(opt.mass2d==opt.mass(ii))=opt.mf(ii);
end

return