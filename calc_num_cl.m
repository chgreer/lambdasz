function ncl = calc_num_cl(mmin,zmin,zmax,area,cosm,opt)
%
%ncl = calc_num_cl(mmin,zmin,zmax,area,cosm,opt)
%
%Calculates the expected number of galaxy clusters
%in a given volume given a cosmology.
%Volume bounded by zmin,zmax, and area (in steradian)

n = opt.NNODES;
hx = (log(opt.MMAX)-log(mmin))/(n-1);
x = log(mmin):hx:log(opt.MMAX);
wx = ones(1,n); wx(2:2:n-1) = 4; wx(3:2:n-2) = 2; wx=wx*hx/3;
hy = (zmax-zmin)/(n-1);
y = zmin:hy:zmax;
wy = ones(1,n); wy(2:2:n-1) = 4; wy(3:2:n-2) = 2; wy=wy*hy/3;

[x2d,y2d] = meshgrid(x,y); x2d=x2d(:); y2d=y2d(:);
w2d = wy.'*wx; w2d=w2d(:).';

ncl = (area/4/pi) * w2d * ncl_integrand(x2d,y2d,cosm,opt);

%vol = (obs.area/4/pi).*co_vol([obs.zmin,obs.zmax],par.cosm,opt,'Mpc');
%ncl = vol * wx * mass_integrand(x,mean([zmin zmax]),cosm,opt);


%zmid = mean([zmin,zmax]);
%vol = (area/4.0/pi) * co_vol([zmin,zmax],cosm,opt,'Mpc'); %in Mpc/h
%ncl = wx * mass_integrand(x,zmid,cosm,opt)' * vol;
return

function ncl_int = ncl_integrand(m,z,cosm,opt)
  m = exp(m);
  a = 1./(1+z);
  ncl_int = zeros(size(m));
  
  for i=1:length(m)
     ncl_int(i) = mass_function(m(i),a(i),cosm,opt);
  end
  uniq_z = unique(z);
  uniq_codVdz = co_dVdz(uniq_z',cosm);
  for i=1:length(uniq_z)
     codVdz(z==uniq_z(i))=uniq_codVdz(i);
  end

  ncl_int = m.*ncl_int.*codVdz';
return

function m_int = mass_integrand(m,z,cosm,opt)
   m = exp(m);
   a = 1./(1+z);
   for i=1:length(m)
        mf(i) = mass_function(m(i),a,cosm,opt);
   end
   m_int = m.*mf;
return
