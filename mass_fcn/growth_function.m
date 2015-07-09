function gf = growth_function(a,cosm,opt)
%
%gf = growth_function(a,cosm,opt)
%
%computes the growth function at a given 
%value of the expansion parameter.
%
%inputs:
%
%a -> expansion parameter for calculation
%cosm -> cosmology structure (help update_cosm)
%opt  -> options structure
%        Requres: opt.TFCHOICE
%                 opt.INTEGRATOR
%           at the very least.
%



if nargin<2
   cosm = get_default_cosm();
end

switch opt.INTEGRATOR
  case 'quad'
     gf = quad(@(x)growth_function_integrand(x,cosm),1.0e-8,a);
  case 'composite'
     n = opt.NNODES; %number of nodes for integration
     h = (a-1.0e-8)/(n-1); %spacing for nodes
     x = (1.0e-8:h:a).'; %node locations
     w = ones(1,n); w(2:2:n-1) = 4; w(3:2:n-2) = 2; w=w*h/3; %weights
     gf = w * growth_function_integrand(x,cosm);
  otherwise
     error('growth_function.m: Integrator not specified');
end

gf = gf .* 2.50 .* cosm.Omega_m .* hz_evol(a,cosm,'expan') .* (1.0 + cosm.zeq);
%gf = gf .* hz_evol(a,cosm,'expan');
%gf = gf .* (1.0 + cosm.zeq);

return

function gfi = growth_function_integrand(a,cosm)

	gfi = a.*hz_evol(a,cosm,'expan');
	gfi = gfi .* gfi .* gfi;
	gfi = 1.0./gfi;

return
