function [log_y,var_log_y] =log_y_give_lam(lam,mean_rln,z,par,opt)
%
%pYgiveLam = calcpYgiveLam
%
%p(Y|lambda) = P(Y,lambda)/P(lambda)
%P(Y,lambda) = \int dM P(Y,lambda|M) P(M)
%P(M) = const * dN/dM;
%
%so
%
%P(Y|lambda) = [\int dM P(Y,lambda|M) dN/dM] / [\int dM P(lambda|M) dN/dM]
%
a = 1./(z+1);
log_l_in = log(lam/60);

%setup integration weights
opt = setup_numint(z,opt,par);

for ii=1:length(log_l_in)
	nrm = opt.wm * lambda_give_mass_ntgrand(log_l_in(ii),mean_rln,par,opt);
  ntegrand =  joint_give_mass_ntgrand(log_l_in(ii),mean_rln,par,opt);
	log_y(ii) = opt.wmy * (opt.log_y2dm .* ntegrand);
	log_y(ii) = log_y(ii)./nrm;
  log_ysq(ii) = opt.wmy * (opt.log_y2dm .* opt.log_y2dm .* ntegrand);
  log_ysq(ii) = log_ysq(ii)./nrm;
end
var_log_y = log_ysq - log_y.^2;

return

function prob = prob_ylam_give_m(log_l,mean_rln,par,opt)
	y_slope = mean_rln(1);	  lam_slope = mean_rln(4);
  y_norml = mean_rln(2);	  lam_norml = mean_rln(5);
  y_scatt = mean_rln(3);	  lam_scatt = mean_rln(6);
  rho =  mean_rln(7);
  log_mstar = log(par.mass_pivot);

  log_lm = lam_norml + lam_slope.*(opt.log_m2dy-log_mstar);
  log_ym = y_norml + y_slope.*(opt.log_m2dy-log_mstar);
  log_x = [ (log_l - log_lm), (opt.log_y2dm-log_ym) ];

	C = [lam_scatt.^2, rho.*lam_scatt.*y_scatt; ...
        rho.*lam_scatt.*y_scatt, y_scatt.^2];
  detC = det(C);

  t = 1-(rho^2);
  q = -(0.5/t)*( (log_x(:,1)/lam_scatt).^2 - ...
    2*rho*log_x(:,1).*log_x(:,2)/lam_scatt/y_scatt + ...
    (log_x(:,2)/y_scatt).^2 );
	prob = exp(q)./(2*pi*sqrt(detC)); %./exp(log_y)./exp(log_lam);
return

function prob = prob_lam_give_m(log_l,mean_rln,par,opt)
  lam_slope = mean_rln(4);
  lam_norml = mean_rln(5);
  lam_scatt = mean_rln(6);
  log_mstar = log(par.mass_pivot);

  log_lm = lam_norml + lam_slope.*(opt.log_m-log_mstar);
  q = -(0.5/lam_scatt/lam_scatt) * (log_l-log_lm).^2;
  prob = exp(q)/sqrt(2*pi)/lam_scatt; %/exp(log_lam);
return

function ntgrand = lambda_give_mass_ntgrand(log_l,mean_rln,par,opt)
	prob = prob_lam_give_m(log_l,mean_rln,par,opt);
  ntgrand = prob.*opt.mass.*opt.mf;
return

function ntgrand = joint_give_mass_ntgrand(log_l,mean_rln,par,opt)
  prob = prob_ylam_give_m(log_l,mean_rln,par,opt);
  ntgrand = prob.*opt.mass2d.*opt.mf2d;
return
