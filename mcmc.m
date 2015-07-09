function [par,obs,opt] = mcmc(fn,verbose)

%if(nargin<2)
%  verbose=1;
%else
%verbose=0; %keep it quiet at first, but then reset once burnin is done
%end

disp('Reading input file.');
[obs,par,opt] = read_mcmc_conf(fn,1);

if(obs.issim)
  obs = gen_fake_obs(obs,par,opt);
else
  obs = read_mcmc_data(obs,par,opt);
end

%open up output file
fid=fopen(opt.fn,'w');

%print parameter names to be fit in output file
for i=1:length(par.fit_par)
  if (par.fit_par(i))
    fprintf(fid,'%s ',par.par_name{i});
  end
end
fprintf(fid,'\n');

%necessary housekeeping
n_accept = 0; %num of accepted mcmc steps
n_rep = 1; %n_rep is the number of times a step is repeated
step_num = 1; %step through on the actual mcmc
nstep = opt.nstep+(opt.burnin)*opt.burnsteps; %total num of loop iter

if(opt.nstep <= opt.burnsteps)
  error('Burnin longer than chain. Stupid.');
end

%set up the chain. first step is initial guess.
par.steps(1,:) = par.init;
par.log_like(1) = calc_like(obs,par,opt); 

%start main loop
tic
fprintf('Optimizing steps for MCMC using %d steps %d times.\n', ...
   opt.burnsteps,opt.burnin);

%this loop indexes over i, but does something funny after
%(opt.burnin)*opt.burnsteps steps, it resets. thus all the arrays
%are indexed by step_num because matlab won't allow the index
%to be modified in the loop
for(i=1:nstep)

 %every burnsteps steps reevaluate the parameter steps
 %do this opt.burnin times then reset
 if(opt.burnin)
   if(mod(step_num,opt.burnsteps)==0)
     fprintf('Retuning after %d steps: %d more times to go. ',...
        opt.burnsteps,opt.burnin-1);
     fprintf('Num_accepted steps: %d (%.2f%%).\n',n_accept, ...
	100*n_accept/opt.burnsteps);
     if(n_accept<0.01*opt.burnsteps)
        fprintf('Something weird happening with this trial..., pausing...');
        %return
     end
     par.covar = cov(par.steps);
     [par.rotn,canon] = eig(par.covar);
     par.eigval = diag(canon); 
     opt.burnin = opt.burnin - 1;
     step_num=1;
     n_accept=0;
     n_rep=1;
     par.steps(1,:) = par.init;
     par.yAlpha = par.init(1); par.yBeta  = par.init(2); 
     par.ySigma = par.init(3);
     par.nAlpha = par.init(4); par.nBeta  = par.init(5); 
     par.nSigma = par.init(6); par.rhoYN = par.init(7);
     par.log_like(1) = calc_like(obs,par,opt);
   end
   if(opt.burnin==0)
      fprintf('Last reop step. Starting MCMC chain now.\n');
      continue;
   end
 end

 %take a step
 step = par.rotn*(par.k.*sqrt(par.eigval).*randn(size(par.init'))).*par.fit_par';;
 par.test = step'+par.steps(step_num,:);

 %set the test parameters
 par.yAlpha = par.test(1); par.yBeta  = par.test(2); par.ySigma = par.test(3);
 par.nAlpha = par.test(4); par.nBeta  = par.test(5); par.nSigma = par.test(6);
 par.rhoYN = par.test(7);  

 %calc the likelihood
 par.log_like_test = calc_like(obs,par,opt);
 
 %metropolis tests against exisiting probability 
 pr = exp(par.log_like_test-par.log_like(step_num));

 testval = rand(1);
 
 %metropolis testing
 if( (pr>testval) )

   %test accepted; go to new val
   par.steps(step_num+1,:) = par.test;
   par.log_like(step_num+1) = par.log_like_test;
   n_accept = n_accept+1;

   %write out old value to file
   if(~opt.burnin)
      for(jj=1:length(par.fit_par))
        if(par.fit_par(jj))
          fprintf(fid,'%f ',par.steps(step_num,jj));
        end
      end
      fprintf(fid,'%d\n',n_rep);
   end
   % reset n_rep
   n_rep=1;
 else
   %test failed stay at old spot
   par.steps(step_num+1,:) = par.steps(step_num,:);
   par.log_like(step_num+1) = par.log_like(step_num);

   %increment n_rep
   n_rep = n_rep+1;
 end

 %print values to screen every so often
 if(mod(step_num,1000)==0)
   fprintf('\n*********** Step %d of %d **************\n',step_num,opt.nstep);
   fprintf('%.4f %.4f %.4f %.4f %.4f %.4f %.4f %d\n',...
      par.steps(step_num,1),par.steps(step_num,2),par.steps(step_num,3),...
      par.steps(step_num,4),par.steps(step_num,5),par.steps(step_num,6),...
      par.steps(step_num,7),n_rep);
   fprintf('log(L_test): %.3f   pr: %.3f   testval: %.3f\n',...
      par.log_like_test,pr,testval);
   fprintf('%d of %d accepted (%.2f percent)\n',n_accept,step_num,...
      n_accept*100/step_num);
   t=toc;
   fprintf('Time run: %.2f min  Time left: %.2f min  Iter: %.2f sec\n',...
      t/60,(nstep-i)*t/i/60,t/i);
 end

 step_num = step_num+1;
 
 if(isnan(par.log_like_test))
    error('Bad likelihood, a nan!')
 end
end

%write out the last step if it hasn't been
if(n_rep>1)
   for(i=1:length(par.fit_par))
     if(par.fit_par(i))
       fprintf(fid,'%f ',par.steps(step_num,i));
     end
   end
   fprintf(fid,'%d\n',n_rep);
end
save(opt.mat_fn);
fclose(fid);

