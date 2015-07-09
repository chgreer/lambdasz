function pivot = find_pivot_point(nreal,ncl,zmin,zmax,area,cosm)

fn_in=tempname;
fn_ot=sprintf('/home/greer/sza_analysis/greer/mf/pivots/%s_%dncl_%.1f-%.1fzrange_%ddeg.pivot',cosm,ncl,zmin,zmax,area);
fid=fopen(fn_in,'w');
fprintf(fid,'[options]\n');

fprintf(fid,'[params]\n');
fprintf(fid,'cosm = %s \n',cosm);

fprintf(fid,'[obs]\n');
fprintf(fid,'issim = 1 \n');
fprintf(fid,'zmin = %f \n',zmin);
fprintf(fid,'zmax = %f \n',zmax);
fprintf(fid,'area = %f \n',area);
fprintf(fid,'ncl = %f \n',ncl);

fclose(fid);

[obs,par,opt]=read_mcmc_conf(fn_in,1);

fprintf('Iteration:');
for i = 1:nreal
  %if(mod(i,25)==0)
     fprintf('.%d',i);
  %end
  obs = gen_fake_obs(obs,par,opt);
  med(i) = median(exp(obs.logm));
end
fprintf('\n');
fid = fopen(fn_ot,'w');
fprintf(fid,'%d',round(median(med)));
fclose(fid);

%for i = 1:length(files)
  %fprintf('Finding pivot for %s\n',files{i});
  %tic
  %[obs,par,opt] = read_mcmc_conf(files{i});
  %for j = 1:nreal
     %if(~mod(j,20))
       %fprintf('realization %d/%d \n',j,nreal);
     %end
     %obs = gen_fake_obs(obs,par,opt);
     %[~,ii] = sort(obs.logn,1,'descend');
     %%exp(obs.logm(ii(1:13)))'
     %m13(j) = obs.logm(ii(13));
  %end
  %pivot(i) = exp(median(m13));
  %piv=round(pivot(i));
  %%cmd=sprintf('sed -n ''H;${x;s/^\\n//;s/prior_pivot = 0, 0\\n/prior_pivot = %.1f,%.1f\\nmass_pivot = %.1f\\n/;p;}'' %s > %s',piv,piv,piv,files{i},files{i});
  %cmd=sprintf('sed ''H;${x;s/^\\n//;s/prior_pivot = 0, 0\\n/prior_pivot = %.1f,%.1f\\nmass_pivot = %.1f\\n/;p;}'' %s > %s',piv,piv,piv,files{i},files{i});
  %system('set +o noclobber');
  %tmp=tempname;
  %cmd=sprintf('sed -n ''H;${x;s/^\\n//;s/prior_pivot = 0, 0\\n/prior_pivot = %.1f,%.1f\\nmass_pivot = %.1f\\n/;p;}'' %s > %s',piv,piv,piv,files{i},tmp);
  %[st,res] = system(cmd);
  %[st,res] = system(sprintf('cp %s %s',tmp,files{i}));
  %fprintf('Pivot point: %.1f x 10^14 Msun, rounded to %.1f x 10^14 Msun\n',pivot(i),round(pivot(i)));
  %toc
  %system('set -o noclobber');
  %fprintf('\n');
  %%disp(cmd);
%end
   

%m13 = zeros(nreal,1);

%for i=1:nreal
%   obs = gen_fake_obs(obs,par,opt);
%   [~,ii] = sort(obs.logn,1,'descend');
%   ii13 = ii(13);
%   m13(i) = obs.logm(ii13);
%end

%pivot = median(m13);
%fprintf('sed -n ''H;${x;s/^\\n//;s/prior_pivot\\n/prior_pivot = %.1f,%.1f\\n&/;p;}'' tmp.in',1.0,1.0);
%end
