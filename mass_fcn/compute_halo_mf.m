function [m mf] = compute_halo_mf(zobs,cosm,opt)
%
%function [m mf] = compute_halo_mf(zobs,cosm,opt)
%
%computes the halo mass function dndm as a function
%of mass, at a given redshift, zobs (optional), for 
%a given cosmology, cosm (optional), for a given
%set of options, opt (optional). Units are in
%terms of 10^12 Msun, physical, not comoving.
%
%inputs:
%
%z -> redshift (default z=0);
%cosm -> cosmology structure (see help update_cosm)
%opt -> options structure to select mass function
%        implementation, integration method to use
%        etc. 
%        
%        Required fields:
%           opt.NUMPTS -> number of points at which
%                          to evaluate mf.
%           opt.MMIN -> minimum mass in 10^12 Msun 
%           opt.MMAX -> maximum mass in 10^12 Msun
%           opt.CHOICE_OF_MF -> mass function implementation
%                                only 'tinker' implemented now
%           opt.TFCHOICE -> transfer function implementation
%                             only 'eisenhu97' implemented now
%           opt.INTEGRATOR -> integration method to use
%                             choices: 
%                                'quad' (slow, accurate)
%                                'composite' (faster, accurate)
%           opt.NNODES -> number of nodes to use for 'composite'
%                          more nodes is slower, but more accurate
%
%


if nargin<3
   opt.NUMPTS = 100;
   opt.MMIN = 1.00;
   opt.MMAX = 1.0e2;
   opt.CHOICE_OF_MF = 'tinker';
   opt.CHOICE_OF_TF  = 'eisenhu97';
   opt.INTEGRATOR = 'composite';
   %opt.INTEGRATOR = 'quad';
   opt.NNODES = 199;
end

if nargin<2
   cosm = get_default_cosm('lcdm',opt);
end

if nargin<1
   zobs = 0;
end

mmax=log10(opt.MMAX);
mmin=log10(opt.MMIN);


a = 1.0/(1.0+zobs);
m = logspace(mmin,mmax,opt.NUMPTS);
mf = 0*m;

for i=1:opt.NUMPTS
   mf(i) = mass_function(m(i),a,cosm,opt);
end
return

