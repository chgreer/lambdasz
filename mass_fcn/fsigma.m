function fs = fsigma(m,a,cosm,opt)

      overD = 200;

      log_alpha = -(0.75/log10(overD/75))^(1.2);
      alpha = 10^(log_alpha);

      tinkerA = 0.185866 * a^(0.14);
      tinkera = 1.466904 * a^(0.06);
      tinkerb = 2.571104 * a^alpha;
      tinkerc = 1.193958;
%      fprintf('a1: %f  a2: %f  a3: %f  a4: %f\n',tinkerA,tinkera,tinkerb,tinkerc);
 
      sm = sigma_M(m,a,cosm,opt); %sqrt(var) of density field

      f_sigma = tinkerA * ( ((sm/tinkerb)^(-tinkera))+1 ) .* exp(-tinkerc./(sm.*sm));

      fs=f_sigma;

return
