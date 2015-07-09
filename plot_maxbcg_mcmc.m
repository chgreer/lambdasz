function h1 = plot_maxbcg_mcmc(fn)

in_fn = strcat(fn,'.in');
out_fn = strcat(fn,'.out');
pdf_fn = strcat(fn,'.pdf');

%read in input file to get priors
[obs,par,opt]=read_mcmc_conf(in_fn,1);

if(nargin<4)
  mass_pivot=[par.mass_pivot,par.mass_pivot];
% mass_pivot=[4.4,4.4];
end

if(nargin<3)
  plot_fn='';
end

%read in data
[data,labels,medd,posd,negd]=read_markov_output(out_fn);
[data] = unwrap_markov_data(data);
npar = size(data,2);

%find prior cov matrix
C = diag(par.prior_widths.^2);
bignum = 1e8;
ind = find(isnan(C));
C(ind) = bignum;

%data(:,2) = data(:,2) + data(:,1)*log(mass_pivot(1)/par.prior_pivot(1));
%data(:,5) = data(:,5) + data(:,4)*log(mass_pivot(2)/par.prior_pivot(2));

fprintf('Adjusting priors to data mass pivot.\n');
fprintf('Prior pivot: %.1f, %.1f\n',par.prior_pivot(1),par.prior_pivot(2));
fprintf('Mass pivot: %.1f, %.1f\n\n',mass_pivot(1),mass_pivot(2));

%get original ones, so we can rotate prior matrix properly
par.prior_cent_val = par.raw_prior_cent_val;
par.prior_widths = par.raw_prior_widths;

fprintf('Before shift:\n ');
for i=1:npar
  fprintf('%s: %.2f+/-%.2f\n ',labels(i,:),par.prior_cent_val(i),par.prior_widths(i));
end
fprintf('\n');

par.prior_cent_val(2) = par.prior_cent_val(2) + ...
	par.prior_cent_val(1)*log(mass_pivot(1)/par.prior_pivot(1));
%par.init(2) = par.init(2) + par.init(1)*log(mass_pivot(1)/par.prior_pivot(1));
%keyboard
%par.prior_widths(2) = par.prior_widths(2) + ...
%	par.prior_cent_val(1)*log(mass_pivot(1)/par.prior_pivot(1));

C22 =C(1:2,1:2);
M = [1,0;log(mass_pivot(1)/par.prior_pivot(1)), 1];
C22 = M*C22*M';
par.prior_width(2) = sqrt(C22(2,2));
C(1:2,1:2) = C22;

par.prior_cent_val(5) = par.prior_cent_val(5) + ...
	par.prior_cent_val(4)*log(mass_pivot(2)/par.prior_pivot(2));
%par.init(5) = par.init(5) + par.init(4)*log(mass_pivot(2)/par.prior_pivot(2));
%par.prior_widths(5) = par.prior_widths(5) + ...
        %par.prior_cent_val(4)*log(mass_pivot(2)/par.prior_pivot(2));

C22 =C(4:5,4:5);
M = [1,0;log(mass_pivot(2)/par.prior_pivot(2)), 1];
C22 = M*C22*M';
par.prior_width(5) = sqrt(C22(2,2));
C(4:5,4:5) = C22;

fprintf('After shift:\n ');
for i=1:npar
  fprintf('%s: %.2f+/-%.2f\n ',labels(i,:),par.prior_cent_val(i),par.prior_widths(i));
  isrho(i) = any(strfind(labels(i,:),'rho'));
end

if(any(isrho))
  npar=npar-1;
end

nsteps=23;
ngauss=101;
h1 = figure('Units','pixels','Position',[100 100 1000 750],...
   'menubar','none','toolbar','none');

ha = tight_subplot(npar,npar,[0.0,0.0],[0.085,0.03],[0.085,0.03]);
ii = 1;

slabels = {strcat('B_{Y}'), strcat('A_{Y}'),...
    strcat('\sigma_{Y}'),...
   'B_{\lambda}','A_{\lambda}','\sigma_{\lambda}','\rho'};


for m = 1:npar 
  for n = 1:npar
    % Plot only on and to left of diagonal
    axes(ha(ii));
    if (m==n | n<m)
      if (n==m)
        %plot histograms
        [counts,xout] = hist(data(:,m),nsteps);
        hbar = bar(xout,counts/max(counts),'hist');
				axis tight;
        xl = xlim;
        xmin(m) = xl(1);
        xmax(m) = xl(2);
        hold on;
        
        %plot gaussian for priors
        h = (xl(2)-xl(1))/(ngauss-1);
        xx= (xl(1)-1):h:(xl(2)+1);
       
        if(isnan(par.prior_cent_val(m)))
          %disp('No prior on this param!'); %continue; %no prior!
        else
          x0 = par.prior_cent_val(m);
          sig2 = par.prior_widths(m)^2;
          gauss = exp( -0.5 * ((xx-x0).^2)./sig2 )./sqrt(2*pi*sig2);
          hl = line([x0 x0],[0 1]);
          hg = plot(xx,gauss/max(gauss));
        end
				set(hbar,'EdgeColor','none','LineWidth',0.2);
        hp=findobj(hbar,'type','patch');
        hfill=hatchfill(hp,'cross',45,3);
        set(hfill,'Color',[0 0 0.8]);
			  hold off;
      else

        %plot contours
        d2=[data(:,m) data(:,n)];
        N = hist3(d2,[nsteps nsteps]);

        %find confidence levels
        for i=1:max(N(:))
          ind = find(N>i);
          enc(i) = sum(sum(N(ind)))/sum(sum(N));
        end
        [~,V(3)] = min(abs(enc-0.68));
        [~,V(2)] = min(abs(enc-0.95));
        [~,V(1)] = min(abs(enc-0.99));
        ninterp=2;
        N = interp2(N,ninterp);

        % Generate grid for 2-D projected view of intensities 
        xb = linspace( min(d2(:,1)) , max(d2(:,1)) , size(N,1));
        yb = linspace( min(d2(:,2)) , max(d2(:,2)) , size(N,1));
        
        % Make a filled contour plot on this grid 
        cmap = [0 0 1.0; 0 0 0.4; 0 0 0.7; 0 0 0];
        h = contourfcmap(yb, xb,N,V,cmap(2:3,:),cmap(4,:),cmap(1,:));

        %find max point to plot and input point to plot
        [~,ind] = max(N(:));
        [xi,yi] = ind2sub(size(N),ind);
        hold on;
        hmax = plot(yb(yi),xb(xi),'wx');

        if(obs.issim)
          hin = plot(par.init(n),par.init(m),'wd');
        end

				%if(m==3)
				%	keyboard
				%end
				%if both params have gaussian priors, draw ellipse
				%keyboard
				%C22 = [C(n,n), C(n,m); C(m,n), C(m,m)];
        %C22 = [C(m,m), C(m,n); C(n,m), C(n,n)];

				%j = 0; k = 0;
				%if( (m==2)  )
					%j = log(par.prior_pivot(1)/mass_pivot(1));
			%		j = log(mass_pivot(1)/par.prior_pivot(1));
				%end
				%if( m==5 )
					%j = log(par.prior_pivot(2)/mass_pivot(2));
				%	j = log(mass_pivot(2)/par.prior_pivot(2));
				%end
				%if( (n==2) )
					%k = log(par.prior_pivot(1)/mass_pivot(1));
			%		k = log(mass_pivot(1)/par.prior_pivot(1));
				%end
				%if (n==5)
					%k = log(par.prior_pivot(2)/mass_pivot(2));
			%		k = log(mass_pivot(2)/par.prior_pivot(2));
				%end
				%M = [1 j; k 1];
				%M = [1 k; j 1];
%keyboard
				%C22 = M*C22*M';
				mu0 = [par.prior_cent_val(n),par.prior_cent_val(m)];
				h0 = plot(mu0(1),mu0(2),'rx','MarkerSize',16,'LineWidth',2);
				if(all(eig(C([n m],[n m]))>0))
					h_pr = error_ellipse(C([n m],[n m]),mu0,0.68);
					set(h_pr,'LineWidth',2,'LineStyle','--','Color','red');
				end
	
				

        hold off;
      end

			%set same xlimits on everything
      xlim([xmin(n) xmax(n)]);

      %plotting done -- now format
      if(exist('hbar','var'))
         set(hbar,'LineWidth',2);
      end
      if(exist('hl','var'))
        set(hl,'Color','r','LineWidth',2,'LineStyle','--');
      end
      if(exist('hg','var'))
        set(hg,'Color','r','LineWidth',2,'LineStyle','--');
      end
      if(exist('hmax','var'))
        set(hmax,'MarkerSize',12,'LineWidth',1.5);
      end
      if(exist('hin','var'))
        set(hin,'MarkerSize',8,'MarkerFaceColor',[1 1 1], ...
					'MarkerEdgeColor','none'); 
      end

      %format ticks so they don't overlap
      %xtck = get(gca,'XTick');
      %ntick = 6;
      %if(length(xtck)<ntick)
        %xl = xlim;
        %h = (xl(2)-xl(1))/(ntick-1);
        %if(h < 0.1)
         %fact = 10;
         %form='%.1f';
        %else
         %fact = 10;
         %form='%.1f';
        %end
        %h = round(h*fact)/fact;
        %xl = sign(xl).*ceil(abs(xl)*fact)/fact;
        %%xlim(xl);
        %xtck = xl(1):h:xl(2);
        %set(gca,'XTick',xtck);
        %xtcklab = num2str(xtck',form);
      %else
        %xtcklab = get(gca,'XTickLabel');
      %end
      %if(length(xtck)>=ntick)
      	%xtcklab(1,:)=' ';
      	%xtcklab(end,:)=' ';
      %end
      %set(gca,'XTickLabel',xtcklab);
%
      %ytck = get(gca,'YTick');
      %if(length(ytck)<ntick)
        %yl = ylim;
        %h = (yl(2)-yl(1))/(ntick-1);
        %if(h < 0.1)
          %fact = 10;
					%form='%.1f';
        %else
          %fact = 10;
					%form='%.1f';
        %end
        %h = round(h*fact)/fact;
        %yl = sign(yl).*ceil(abs(yl)*fact)/fact;
        %%ylim(yl);
        %%ytck = yl(1):h:yl(2);
        %set(gca,'YTick',ytck);
        %ytcklab = num2str(ytck',form);
      %else
        %ytcklab = get(gca,'YTickLabel');
				%ytcklab(1,:)=' ';
				%ytcklab(end,:)=' ';
      %end
      %set(gca,'YTickLabel',ytcklab);

      %y label only the left side     
      if(n==1)
        hYLabel=ylabel(slabels{m});
        if(m==1)
          hYLabel=ylabel('');
					set(gca,'YTickLabel',[]);
					set(gca,'YTick',[]);
        end
      else
        hYLabel=ylabel('');
        set(gca,'YTickLabel',[]);
      end
 
     %xlabel only the bottom
     if(m==npar)
        hXLabel=xlabel(slabels{n});
        set(hXLabel,'Position',get(hXLabel,'Position')-0.02);
      else
        hXLabel=xlabel('');
        set(gca,'XTickLabel',[]);
      end

      %set general format properties
      set([hXLabel, hYLabel],'FontName', 'Helvetica','FontSize',14);
      set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.04 .04] , ...
        'LineWidth'   , 1         );
      %set(gcf, 'Color', 'none'); % Sets figure background
      set(gca, 'Color', 'none'); % Sets axes background
    else
       axis off;
    end
    ii=ii+1;
  end
end

info{1}=strcat('Run Info:');
if(obs.issim == 1)
  info{end+1} = 'Simulation';
else
  info{end+1} = 'Real Data!';
end
info{end+1}=strcat('Obs: ',obs.type);
info{end+1}=strcat('N_{cl}:',num2str(obs.ncl));
info{end+1}=strcat('Cosm: ',par.cosm.name);
info{end+1}=strcat('Area: ',num2str(obs.area*3282.8),'sq deg');
if(obs.issim ==1)
	info{end+1}=strcat('\sigma_{obs} : ',num2str(obs.sim_yerr));
end
info{end+1}=strcat('Obs M_*: ',num2str(par.mass_pivot),' \times 10^{14} Msun');
info{end+1}=strcat('Prior M_*: ',num2str(par.prior_pivot(1)),' \times 10^{14} Msun');
%info{end+1}=strcat('A_{\lambda} Prior Width:',num2str(par.prior_widths(5)));
%annotation('textbox',[0.80 0.85 0.1 0.1],'String',info,'FontSize',13,...
%	'FitBoxToText','on','LineStyle','none');

res{1}=strcat('Constraints:');
res{end+1}='';
res{end+1}=strcat('B_{',obs.type,'} = ',num2str(medd(1),'%.2f'),...
	'^{+',num2str(abs(posd(1)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(1)),'%.2f'),'}');
res{end+1}='';
res{end+1}=strcat('A_{',obs.type,'} = ',num2str(medd(2),'%.2f'),...
        '^{+',num2str(abs(posd(2)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(2)),'%.2f'),'}');
res{end+1}='';
res{end+1}=strcat('\sigma_{',obs.type,'} = ',num2str(medd(3),'%.2f'),...
        '^{+',num2str(abs(posd(3)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(3)),'%.2f'),'}');
res{end+1}='';
res{end+1}=strcat('B_{\lambda} = ',num2str(medd(4),'%.2f'),...
        '^{+',num2str(abs(posd(4)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(4)),'%.2f'),'}');
res{end+1}='';
res{end+1}=strcat('A_{\lambda} = ',num2str(medd(5),'%.2f'),...
        '^{+',num2str(abs(posd(5)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(5)),'%.2f'),'}');
res{end+1}='';
res{end+1}=strcat('\sigma_{\lambda} = ',num2str(medd(6),'%.2f'),...
        '^{+',num2str(abs(posd(6)),'%.2f'),'}',...
        '_{-',num2str(abs(negd(6)),'%.2f'),'}');

%annotation('textbox',[0.65 0.85 0.1 0.1],'String',res,'FontSize',13,...
%	'FitBoxToText','on','LineStyle','none');

eqn{1} = 'Mean Relations:';
eqn{end+1} = '<\lambda|\lambda_*> = A_{\lambda} + B_{\lambda} log(M/M_*)';
eqn{end+1} = '<Obs|Obs_*> = A_{obs} + B_{obs} log(M/M_*)';
if(obs.issim)
	eqn{end+1} = 'White diamond is sim input.';
end
eqn{end+1} = 'White X is max of 2d histogram.';
%annotation('textbox',[0.35 0.85 0.1 0.1],'String',eqn,'FontSize',13,'FitBoxToText','on','LineStyle','none');


%save to plot
if(pdf_fn)
  export_fig(pdf_fn,'-pdf','-eps','-transparent');
end

return
