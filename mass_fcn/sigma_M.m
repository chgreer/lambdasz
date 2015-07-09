function sm = sigma_M(m,a,cosm,opt)
%
%calculates the sqrt(var) of the matter density field

rc = rho_crit(1.0,cosm,'tinker','expan');

r = ( (3.0.*m)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);
%fprintf('mass: %e \t R: %f \n',m*1e14,r);
sm = sigma_R(r,a,cosm,opt);

return

