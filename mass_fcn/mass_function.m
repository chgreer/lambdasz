function dndm = mass_function(m,a,cosm,opt)
%
%function dndm = mass_function(m,a,cosm,opt)
%
%returns the differential mass function dndm
%at a single value of mass, m,
%as a function of expansion parameter, a,
%cosmology, cosm, with options, opt.
%
%dndm in halos/10^14 Msun/Mpc^3 relative to 
%M_{500m}.
%
%m -> mass in 10^14 Msun M500c -- CRITICAL!!!
%a -> expansion parameter 
%cosm -> cosmology structure (see help update_cosm)
%opt -> options structure to select mass function
%        implementation, integration method to use
%        etc. see help compute_halo_mf for example 

%internally uses 1e14 suns

z = 1/a-1;

%works internally on mean, so convert:
m = crit2mean(m,500,z,cosm);

rhoc = rho_crit(1.0,cosm,'tinker','expan');
   
switch opt.CHOICE_OF_MF
   case 'tinker'
      %tinker mass function
      %see Astrophys.J.688:709-728,2008

      overD = 500;
			
      switch overD
        case 500
          tinkerA=0.2150;
					tinkera=1.5850;
          tinkerb=1.9600;
          tinkerc=1.3950;
        otherwise
		      Delta_grid = [200 300 400 600 800 1200 1600 2400 3200];
      		A_grid = [0.186 0.200 0.212 0.218 0.248 0.255 0.260 0.260 0.260];
      		a_grid = [1.47 1.52 1.56 1.61 1.87 2.13 2.30 2.53 2.66];
					b_grid = [2.57 2.25 2.05 1.87 1.59 1.51 1.46 1.44 1.41];
					c_grid = [1.19 1.27 1.34 1.45 1.58 1.80 1.97 2.24 2.24];
		
      		tinkerA = interp1(Delta_grid,A_grid,overD);
      		tinkera = interp1(Delta_grid,a_grid,overD);
      		tinkerb = interp1(Delta_grid,b_grid,overD);
      		tinkerc = interp1(Delta_grid,c_grid,overD);
      end
      log_alpha = -(0.75/log10(overD/75))^(1.2);
      alpha = 10^(log_alpha);

%      tinkerA = 0.185866 * a^(0.14);
%      tinkera = 1.466904 * a^(0.06); 
%      tinkerb = 2.571104 * a^alpha;
%      tinkerc = 1.193958;

      tinkerA = tinkerA * a^(0.14);
      tinkera = tinkera * a^(0.06); 
      tinkerb = tinkerb * a^alpha;
      tinkerc = tinkerc;

      sm = sigma_M(m,a,cosm,opt); %sqrt(var) of density field
      ds2dm = dsigma2_dM(m,a,cosm,opt); %direvative of sqrt(var)

      f_sigma = tinkerA * ( ((sm/tinkerb)^(-tinkera)) +1 ) ...
				.* exp(-tinkerc./(sm.*sm));

      dndm = - f_sigma * cosm.Omega_m.*rhoc.*ds2dm./sm./sm./m/2;

   otherwise
      disp('MF Not implimented.');
end

%output at 500c, not 500m
%z = 1/a - 1;
%scal = m/conv_Mdelta(m,z,500,500,'mean','crit');
%dndm = dndm/scal;

return


function sm = sigma_m(m,a,cosm,opt)
%
%calculates the sqrt(var) of the matter density field
 
rc = rho_crit(1.0,cosm,'tinker','expan');

r = ( (3.0.*m)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);
%fprintf('mass: %e \t R: %f \n',m*1e14,r);
sm = sigma_R(r,a,cosm,opt);

return

function dsdm = dsigma_dm(m,a,cosm,opt)
%
%calculate the derivative of the sqrt(var) of the
%matter density field

rc = rho_crit(1.0,cosm,'tinker','expan');

%temp=(3.0*m)/(4.0*pi*cosm.Omega_m*rc); 
%temp=(temp^(1.0/3.0));

r = ( (3.0.*m)./(4.0.*pi.*cosm.Omega_m.*rc) ).^(1/3);


switch opt.INTEGRATOR
   case 'quad'
      dsdm = quad(@(x)dsigma_dR_integrand(x,a,r,cosm,opt),-0.95,0.95,1e-6);
   case 'composite'
      n = opt.NNODES; %number of nodes for integration
      h = (0.95+0.95)/(n-1); %spacing for nodes
      x = (-0.95:h:0.95).'; %node locations
      w = ones(1,n); w(2:2:n-1) = 4; w(3:2:n-2)= 2; w=w*h/3; %weights
      dsdm = w * dsigma_dR_integrand(x,a,r,cosm,opt); %integral
   otherwise
     error('mass_function.m: function dsdm -- integration method not specified'); 
end

dsdm = dsdm * 1.0/3.0*(r./m);

return


function dsdRi = dsigma_dR_integrand(x,a,R,cosm,opt)

k=exp( x./(1.0-x.*x) );
dm=Delta_m(k,a,cosm,opt);

kR=k.*R;

w = tophat_window(kR);

dwdr = 3.*k.*((kR.^2-3)*sin(kR)+3*kR*cos(kR))/(kR.^4);
dw2dr = 2.*w.*dwdr;

dsdRi = dm * dw2dr * (1.0+x.*x)./(1.0-x.*x)./(1.0-x.*x);

%temp1=tophat_window(kR);
%temp2=3.0./kR.*( sin(kR)./kR - temp1 );% /* [window function] */
%val=val.*2.0.*temp1.*temp2.*k;

%val=val.*(1.0+x.*x)./(1.0-x.*x)./(1.0-x.*x); %/* Measure for integration variable. */
%dsdRi = val;

return;






