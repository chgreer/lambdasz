function [obs,par,opt] = read_mcmc_conf(fn,skipncl);
%
%[obs,par,opt] = read_mcmc_conf(fn,skipncl)
%
%reads in the parameter file fn and returns
%the three structures that run the mcmc routine
%
%skipncl is an option flag to skip the lengthy
%calculation for ncl_expect if it's not needed
%default: off
%

if nargin<2
  skipncl=0;
end

%initialize structures 
obs.name = 'obs';
opt.name = 'opt';
par.name = 'par';


opt_keys = {'options', '', 'nstep', 'i', 10; ...
            'options', '', 'burnin', 'i', 3;, ...
            'options', '', 'burnsteps', 'd', 100;, ...
			    	'options', '', 'method', '', 'full'; ...
	  			  'options', '', 'fn', '', '/tmp/junk.out'; ...
            'options', '', 'MMIN', 'd', 1; ...
            'options', '', 'MMAX', 'd', 200; ...
            'options', '', 'LMIN', 'i', 30; ...
            'options', '', 'LMAX', 'i', 1000; ...
            'options', '', 'YMIN', 'd', 0.1; ...
	    			'options', '', 'YMAX', 'd', 1000; ...
            'options', '', 'NUMPTS', 'i', 2000; ...
            'options', '', 'CHOICE_OF_MF', '', 'tinker'; ...
            'options', '', 'CHOICE_OF_TF', '', 'eisenhu97'; ...
            'options', '', 'INTEGRATOR', '', 'composite'; ...
            'options', '', 'NNODES', 'i', 161};

par_keys = {'params', '', 'init', 'd', [1.666,0.0,0.166,1.0,0.0,0.25,0.0]; ...
            'params', '', 'k', 'd', 0.8; ...
            'params', '', 'var', 'd', 0.03; ...
            'params', '', 'covar_fn', '', '';... 
            'params', '', 'fit_par', 'd', [1,1,1,1,1,1,0]; ...
            'params', '','prior_cent_val','d',[1.666,0,0.166,1.0,0,0.25,nan]; ...
	    			'params', '', 'prior_widths', 'd', ones(1,7); ...
            'params', '', 'mass_pivot', 'd', ones(1,1); ...
            'params', '', 'prior_pivot', 'd', ones(1,2); ...
            'params', '', 'flat_var_priors', 'd', zeros(1,2); ...
            'params', '', 'cosm', '', 'lcdm' };

obs_keys = {'obs', '', 'zmin', 'd', 0.2; ...
            'obs', '', 'zmax', 'd', 0.3; ...
            'obs', '', 'area', 'd', 7398; ...
            'obs', '', 'ncl', 'i', 28; ...
            'obs', '', 'type', '', 'ysz'; ...
            'obs', '', 'issim', 'i', 1; ...
						'obs', '', 'rad', 'i', 500; ...
						'obs', '', 'simcosm', '', ''; ...
						'obs', '', 'simpivot', 'd', []; ...
            'obs', '', 'sim_yerr', 'd', 0.1; ...
            'obs', '', 'missing_y', 'd', []; ...
            'obs', '', 'data_fn', '', 'greer/mf/maxbcg_data.txt'};

%read the options
keys = opt_keys(:,3);
conf = inifile(fn,'read',opt_keys);

for i=1:length(keys)
   opt = setfield(opt,keys{i},conf{strcmp(keys,keys{i})});
end

%do param struct
keys = par_keys(:,3);
conf = inifile(fn,'read',par_keys);

for i=1:length(keys)
   par = setfield(par,keys{i},conf{strcmp(keys,keys{i})});
end

%do obs struct
keys = obs_keys(:,3);
conf = inifile(fn,'read',obs_keys);

for i=1:length(keys)
   obs = setfield(obs,keys{i},conf{strcmp(keys,keys{i})});
end

%temporary filename that allows run-once generation of the mass function file
opt.mf_fn = tempname;
opt.mf2d_fn = tempname;

%fn for a .mat save file of run data
last_dot = strfind(opt.fn,'.');
last_dot = last_dot(end);
opt.mat_fn = opt.fn(1:last_dot); %strtok(opt.fn,'.');
opt.mat_fn = strcat(opt.mat_fn,'mat');

%don't fit some parameters
par.fit_par = logical(par.fit_par);
par.par_name = {'slope_y()','norm_y()','sig_y()','slope_n()','norm_n()',...
	'sig_n()','rho()'};

if(obs.rad == 1)
  obs.rad = '1Mpc';
end

%filename for a .mat save file of parameter covariances
if(exist(par.covar_fn,'file'))
   %fprintf('Loading parameter covariances from %s\n',par.covar_fn);
   %par.covar = load(par.covar_fn,'-ascii');
   [data] = read_markov_output(par.covar_fn);
   [data] = unwrap_markov_data(data);
   par.covar = cov(data);
else
   %fprintf('No parameter covariance file. Using diagional matrix.\n');
   %par.covar = diag(par.prior_widths);
   %par.covar(isnan(par.covar)) = 0; %don't step in parameters with nan
   par.covar=par.var*diag(par.fit_par);
   %this breaks if init params are 0. no steps taken
   %par.covar = diag((par.var.*par.init).^2);
end
[par.rotn,canon] = eig(par.covar);
par.eigval = diag(canon);

%do some calculations
par.cosm = get_default_cosm(par.cosm,opt);
obs.area = obs.area/3282.8; %into str

%set up the priors
if length(par.var)<length(par.init)
   par.var = par.var(1) .* par.init;
end

par.raw_prior_cent_val = par.prior_cent_val;
par.raw_prior_widths = par.prior_widths;

%convert the input prior to the mass_pivot point we're using

%fprintf('Adjusting priors to data mass pivot.\n');
%fprintf('Prior pivot: %.1f, %.1f\n',par.prior_pivot(1),par.prior_pivot(2));
%fprintf('Mass pivot: %.1f, %.1f\n\n',par.mass_pivot,par.mass_pivot);

%npar=7
%fprintf('Before shift:\n ');
%for i=1:length(par.par_name);
%  fprintf('%s: %.2f+/-%.2f\n ',par.par_name{i},par.prior_cent_val(i),par.prior_widths(i));
%end
%fprintf('\n');

%ynorm_cent_val = par.prior_cent_val(2) + ...
	%par.prior_cent_val(1)*log(par.mass_pivot/par.prior_pivot(1));

%find prior cov matrix
%C = diag(par.prior_widths.^2);
%prior_cov_mat = C;
%bignum = 1e8;
%ind = find(isnan(C));
%C(ind) = bignum;
%
%C22 =C(1:2,1:2);
%%M = [1,0;log(par.mass_pivot/par.prior_pivot(1)), 1];
%C22 = M*C22*M';
%C22y = C22;
%ynorm_prior_width = sqrt(C22(2,2));
%
%lnorm_cent_val = par.prior_cent_val(5) + ...
  %par.prior_cent_val(4)*log(par.mass_pivot/par.prior_pivot(2));
%
%C22 =C(4:5,4:5);
%M = [1,0;log(par.mass_pivot/par.prior_pivot(2)), 1];
%C22 = M*C22*M';
%C22l = C22;
%lnorm_prior_width = sqrt(C22(2,2));
%
%
%if(isnan(par.prior_cent_val(2))==0)
	%par.prior_cent_val(2) = ynorm_cent_val;
	%par.prior_widths(2) = ynorm_prior_width;
	%prior_cov_mat(1:2,1:2) = C22y;
%end
%if(isnan(par.prior_cent_val(5))==0)
	%par.prior_cent_val(5) = lnorm_cent_val;
  %par.prior_widths(5) = lnorm_prior_width;
  %prior_cov_mat(4:5,4:5) = C22l;
%end
%
%fprintf('After shift:\n ');
%for i=1:npar
  %fprintf('%s: %.2f+/-%.2f\n ',par.par_name{i},par.prior_cent_val(i),par.prior_widths(i));
%end
%ii = find(isnan(par.prior_cent_val)==0);
%par.prior_cov_mat = prior_cov_mat(ii,ii)



prior_cov_mat = diag(par.prior_widths.^2);
ii = find(isnan(par.prior_cent_val)==0);
par.prior_cov_mat = prior_cov_mat(ii,ii);

%if(isreal(ynorm_cent_val))
%  prior_cov_mat(1:2,1:2) = C22y;
%end
%if(isreal(lnorm_cent_val))
%	prior_cov_mat(4:5,4:5) = C22l;
%end

%if parameter is free from gaussian prior, set it to nan in conf file
%prior_cent_val = par.prior_cent_val(~isnan(par.prior_cent_val));
%prior_widths = par.prior_widths(~isnan(par.prior_cent_val));

%if length(par.prior_widths)<length(prior_cent_val)
   %par.prior_widths = par.prior_widths(1) .* ones(size(prior_cent_val));
%end

%sigma = prior_widths; %prior_cent_val .* prior_widths;
%par.prior_cov_mat = diag(sigma.*sigma);

%if(isreal(ynorm_cent_val) && isreal(par.prior_cent_val(1)))
  %par.prior_cov_mat(1:2,1:2) = C22y;
%end
%if(isreal(lnorm_cent_val) && isreal(par.prior_cent_val(4)))
	%par.prior_cov_matI

par.yAlpha = par.init(1); par.yBeta  = par.init(2); par.ySigma = par.init(3);
par.nAlpha = par.init(4); par.nBeta  = par.init(5); par.nSigma = par.init(6);
par.rhoYN = par.init(7);

%read in data or generate the fake stuff
%if (~skipncl)
%   obs.ncl_expect = round(calc_num_cl(opt.MMIN,obs.zmin,obs.zmax,obs.area,par.cosm,opt));
%   %fprintf('Expecting %d clusters in this cosmology, etc.\n',obs.ncl_expect);
%else
%  obs.ncl_expect=1000;
%end

obs.vol = (obs.area/4/pi).*co_vol([obs.zmin,obs.zmax],par.cosm,opt,'Mpc');

if (obs.issim)
	if (isempty(obs.simcosm))
	  obs.simcosm = par.cosm;
	else
		obs.simcosm = get_default_cosm(obs.simcosm,opt);
	end
	if (isempty(obs.simpivot))
		obs.simpivot = par.mass_pivot;
	end
end
  %disp('Generating fake data.');
%  obs=gen_fake_obs(obs,par,opt);
%else
  %%disp('Reading data.')
  %[obs.name, obs.zobs, lambda, lam_err, y, yerr_p, yerr_m] = ...
		%textread(obs.data_fn,'%s %f %f %f %f %f %f','commentstyle','shell');
  %ii=find(y>0);
  %obs.missing_y = find(y<=0); 
  %obs.logn = log(lambda/60);
  %obs.logy = log(y/1e-5);
  %obs.logn_err = lam_err./lambda;
  %obs.logy_err = (yerr_p-yerr_m)./(2*y);
  %obs.ncl = length(y);
%end
%obs.aobs = 1./(1+obs.zobs);

 

return
