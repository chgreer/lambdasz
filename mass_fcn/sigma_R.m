function sR = sigma_R(R,a,cosm,opt)

%fprintf('Calling sigma_R.m.\n');
%fprintf('R: %f     a: %f\n',R,a);

switch opt.INTEGRATOR
 case 'quad'
    fprintf('Calling quad to integrate.\n');
    int_val =quad(@(x)sigma_R_integrand(x,R,a,cosm,opt),-0.95,0.95);
 case 'composite'
    n = opt.NNODES; %number of nodes for integration
    h = (0.95+0.95)/(n-1); %spacing for nodes
    x = (-0.95:h:0.95).'; %node locations
    w = ones(1,n); w(2:2:n-1) = 4; w(3:2:n-2) = 2; w=w*h/3; %weights
    int_val = w * sigma_R_integrand(x,R,a,cosm,opt);
 otherwise
     error('sigma_R -- integration method not specified');
end

sR = sqrt(int_val);

return

function sRint = sigma_R_integrand(x,R,a,cosm,opt)

k = exp( x./(1.0-x.*x) );
w = tophat_window(k.*R);

%fprintf('R: %f   a: %f \n',R,a);

val = Delta_m(k,a,cosm,opt) .* w .* w;
val = val .* (1.0+x.*x)./(1.0-x.*x)./(1.0-x.*x);

sRint = val;

return

