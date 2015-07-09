function dm = Delta_m(k,a,cosm,opt)

gf_scale = growth_function(a,cosm,opt)/cosm.gf0;

dm =  cosm.psAmp*(k.^(3.0+cosm.tilt))/2.0/pi/pi;
dm = dm.*trans_fun(k,a,cosm,opt).^2.*gf_scale.^2;

return


