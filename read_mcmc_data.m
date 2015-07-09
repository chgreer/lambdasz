function obs = read_mcmc_data(obs,par,opt)

[obs.name, obs.lambda, obs.lam_err, obs.zobs, ...
	r500, r500p, r500m, y500, y500p, y500m, ...
  r2500, r2500p, r2500m, y2500, y2500p, y2500m, ...
	y1Mpc,y1Mpcp,y1Mpcm] = textread(obs.data_fn,...
  	'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
  	'commentstyle','shell');

ii = find(y1Mpc<=0);
obs.name(ii)
obs.missing_y = ii;

switch(obs.rad)
	case 500
  	r = r500;
  	rp = r500p;
  	rm = r500m;
  	y = y500;
  	yp = y500p;
  	ym = y500m;
	case 2500
  	r = r2500;
  	rp = r2500p;
  	rm = r2500m;
  	y =  y2500;
  	yp = y2500p;
  	ym = y2500m;
	case '1Mpc'
  	r = ones(size(y1Mpc));
  	rp = zeros(size(r));
  	rm = zeros(size(r));
  	y = y1Mpc
  	yp = y1Mpcp;
  	ym = y1Mpcm;
	otherwise
  	error('yvol_vs_lam.m :: radius not supported');
end

obs.y = y;
obs.yp = yp;
obs.ym = ym;
obs.logn = log(obs.lambda/60);
obs.logy = log(y);
obs.logn_err = obs.lam_err./obs.lambda;
obs.logy_err = (abs(yp)+abs(ym))./(2.*y);
obs.ncl = length(y);

obs.aobs = 1./(1+obs.zobs);

return