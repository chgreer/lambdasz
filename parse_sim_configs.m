function parse_sim_configs()

d='.';
list = dir([d,sprintf('/*out')]);
files = {list.name};
fprintf('Found %d files\n',length(files));

for i = 1:length(files)
  file = files{i};
  fprintf('Parsing %s \n',file);
  ind = strfind(file,'.');
  ind = max(ind);
  fn = file(1:ind-1);
  fnin = strcat(file(1:ind),'in');
  fnout = strcat(file(1:ind),'out');
  fnpdf = strcat(file(1:ind),'pdf');
  fnpdfdat=strcat(file(1:ind-1),'_data.pdf'); 
  fprintf('Data pdf: %s\n', fnpdfdat);
	if strfind(fn,'rozo')
		datafn = '/home/greer/sza_analysis/greer/maxbcg/maxbcg_ys_rozo.txt';
	elseif (strfind(fn,'arnaud'))
		datafn = '/home/greer/sza_analysis/greer/maxbcg/maxbcg_ys_arnaud.txt';
  else
		datafn = '';
	end

	if strfind(fn,'R500')
		rad=500;
	elseif strfind(fn,'R2500')
		rad=2500;
	elseif strfind(fn,'1Mpc')
		rad='1Mpc';
	else
		rad=500;
	end

  data = read_markov_output(fnout);
  n_uniq = length(data);
  data = unwrap_markov_data(data);
  n_step = length(data);
  accept_frac(i) = n_uniq/n_step;


  ind = strfind(file,'_');
  ind = max(ind);
  fnrt{i} = file(1:ind-1);
  %don't plot janky pdfs
	if ( accept_frac > 0.04 ) 
    if(~exist(fnpdf,'file'))
	  	h = plot_maxbcg_mcmc(fn);
  	  close(h);
		end
		
		%if(exist(datafn,'file'))
		%	fprintf('Plotting %s\n',fnpdfdat);
		%	if(~exist(fnpdfdat,'file'))
		%		h = yvol_vs_lam(datafn,rad,fnin,fnpdfdat);
		%		close(h);
		%	else
		%		fprintf('%s already exists!!\n',fnpdfdat);
		%	end
		%end

	end
end

[fnrt,~,ind] = unique(fnrt);
simmed=[];
simmedvar=[];
simpos=[];
simneg=[];

fid=fopen('sim_out_data.txt','w');
for i = 1:length(fnrt)
  ind_i = find(ind==i);
  list = dir([d,sprintf('/%s*out',fnrt{i})]);
  reals = {list.name};
  medval = [];
  poserr = [];
  negerr = [];
  for j = 1:length(reals)
    if(accept_frac(ind_i(j))>0.1)
	    [data,labels,xx,yy,zz] = read_markov_output(reals{j},0,0);
  	  medval = [medval; xx];
    	poserr = [poserr; yy];
    	negerr = [negerr; zz];
    end
  end
  nreal = size(medval,1);
	%stdmedval = std(medval);
  ii = randperm(nreal);
  ii = ii(end);
	medval = medval(ii,:);
	poserr = poserr(ii,:);
	negerr = negerr(ii,:);

	if length(medval)==6
		medval(end+1)=0;
		%stdmedval(end+1)=0;
		poserr(end+1)=0;
		negerr(end+1)=0;
	end		
  %keyboard
  simmed = [simmed; medval];
  %simmedvar = [simmedvar; stdmedval];
  simpos = [simpos; poserr];
  simneg = [simneg; negerr];
  fprintf(fid,'%s',fnrt{i});
	%keyboard
  for j = 1:length(simmed(i,:))
    fprintf(fid,' %.3f %.3f %.3f',simmed(i,j),simpos(i,j),...
			simneg(i,j));
  end
  fprintf(fid,'\n');
end

fclose(fid);   

end
