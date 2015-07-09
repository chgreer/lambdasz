function rnd = corr_rand(n,covm,U)
%generates correlated, normally distributed random numbers
%see: http://www.sitmo.com/doc/Generating_Correlated_Random_Numbers
%
%Takes a covariance matrix cov and generates n observations with size
%equal to rank cov.
%
%E.g. >> 
%cov = [0.5 sqrt(0.5*0.1)*0.2; sqrt(0.5*0.1)*0.2 0.1]
%rnd = corr_rand(1000000,cov);
%
%rnd will now have 1000000 pairs of variables if cov is 2x2.

if(nargin<3)
  %get the eigenvector matrix V and the eigen value matrix D.
  [V,D] = eigs(covm);
  
  %generate the transformation matrix.
  U = (V*sqrt(D))';
end

%transform independent observations into correlated ones
rnd = randn(n,length(U))*U;

return
