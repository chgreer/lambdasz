function tfn = trans_fun(k,a,cosm,opt)
%
%function tfn = trans_fun(k,a,cosm,opt)
%
%computes the transfer function that takes the
%primordial power spectrum into the matter power
%spectrum.
%
%inputs:
%
%k -> fourier mode (Mpc^-1).
%a -> expansion parameter for calculation
%cosm -> cosmology structure (help update_cosm)
%opt  -> options structure
%        Requres: opt.TFCHOICE
%                 opt.INTEGRATOR
%           at the very least.
%



switch opt.CHOICE_OF_TF
   case 'eisenhu97'
      %See astro-ph/9710252v1 section 3.3
      tfn = low_baryon_matter_tf(k,a,cosm,opt);
   otherwise
      error('trans_fun.m: CHOICE_OF_TF not valid');
end

return

function lbmtf = low_baryon_matter_tf(k,a,cosm,opt)
  lbmtf = master_tf(k,cosm); 
  if(cosm.fnu>1e-6)
     %corrections only necessary if fnu>0
     lbmtf = lbmtf.*growth_cbnu(k,a,cosm,opt)/growth_function(a,cosm,opt);
  end
return

function mtf = master_tf(k,cosm)
   q = k./cosm.keq;

   temp=0.43.*k.*cosm.shor;
   temp=temp.*temp.*temp.*temp;
   Geff=cosm.om_m.*( sqrt(cosm.alpha_nu)+(1.0-sqrt(cosm.alpha_nu))./(1.0+temp) );

   qeff=k.*cosm.Tcmb.*cosm.Tcmb./Geff;

   L=log( exp(1.0)+1.84.*qeff.*sqrt(cosm.alpha_nu)./(1.0-0.949.*cosm.fnub) );
   C=14.4+325.0./(1.0+60.5.*(qeff.^(1.11)));
   T=L./(L+C.*qeff.*qeff);

   temp=3.92.*q.*sqrt(cosm.nNu/cosm.fnu);
   Bcorr=1.2.*(cosm.fnu.^(0.64)).*(cosm.nNu.^(0.3+0.6.*cosm.fnu));
   Bcorr=Bcorr ./ ( (temp.^(-1.6))+(temp.^(0.8)) );
   Bcorr=1.0+Bcorr;
   mtf=T.*Bcorr;
return

function gcbnu = growth_cbnu(k,a,cosm,opt)

   q=k/cosm.keq;
   yfs=(cosm.nNu*q/cosm.fnu);
   yfs=yfs.*yfs;
   yfs=yfs.*17.2*cosm.fnu*( 1.0+0.488*(cosm.fnu^(-7.0/6.0)) );

   temp = growth_function(a,cosm,opt);

   val=((temp./(1.0+yfs)).^(0.7));
   val = val + (cosm.fcb^(0.7/cosm.pcb));
   val=(val.^(cosm.pcb/0.7));
   val=val.*(temp.^(1.0-cosm.pcb));
   gcbnu = val;

return


