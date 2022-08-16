% plot_phase_retrieval_betas.m -- This function plots the reconstruction
% quality trends for the proposed method as a function of regularization
% parameter (beta). This code loads results files created by
% run_phase_retrieval_betas.m. You should run that first. This figure is
% included in Fig. 3 of the paper.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% files to plot
filepath = fullfile('..','..','results','phase_retrieval');
fileprefix = 'betas_';
filesuffix = '.mat';

[filenames,filepath] = uigetfile({[fileprefix '*' filesuffix];['*' filesuffix]},'Plot phase retrieval betas...',filepath,'MultiSelect','on');
if isempty(filenames) || isequal(filenames,0) || isequal(filepath,0), return; end
if ~iscell(filenames), filenames = {filenames}; end

Nplots = length(filenames);

outfilepath = uigetdir('Save betas results to...');
if isequal(outfilepath,0), return; end

%% do plotting of median, deciles (10-90%)
boxes = true; %boxes = false; % set to true to plot quartile boxes; otherwise, just medians are plotted

min_error_use = cell(Nplots,1);
min_error_betas = cell(Nplots,1);

for iplot = 1:Nplots
    Sdata = load(fullfile(filepath,filenames{iplot}),'errors','Ms','Ks','N','betas','Ntrials','betas_grid');
    
    if min(diff(Sdata.betas))*(log10(Sdata.betas(end))-log10(Sdata.betas(1))) < min(diff(log10(Sdata.betas)))*(Sdata.betas(end)-Sdata.betas(1))
        scale = 'log';
    else
        scale = 'linear';
    end
    
    errors_sorted = sort(Sdata.errors.^2,1);
    errors_quartiles = reshape(interp1((0:Sdata.Ntrials-1).'+0.5,reshape(errors_sorted,Sdata.Ntrials,[]),[0.25,0.5,0.75]*Sdata.Ntrials,'linear','extrap'),[3,length(Sdata.Ks),length(Sdata.Ms),length(Sdata.betas)]);
    errors_quartiles(errors_quartiles < 0) = 0;
    errors_use = reshape(errors_quartiles(2,:,:,:),[length(Sdata.Ks),length(Sdata.Ms),length(Sdata.betas)]);
    [min_error_use{iplot},min_error_betas{iplot}] = min(errors_use,[],3);
    min_error_betas{iplot} = Sdata.betas(min_error_betas{iplot});
    
    if length(Sdata.Ks) < length(Sdata.Ms)
        inplot = Sdata.Ms;
        outplot = Sdata.Ks;
        errors_use = permute(errors_use,[2,1,3]);
        errors_quartiles = permute(errors_quartiles,[1,3,2,4]);
        Sdata.betas_grid = permute(Sdata.betas_grid,[2,1,3]);
        namestr = '%s: K = %d';
        titlestr = 'M = %d';
        filenamestr = '%s_K%d';
    else
        inplot = Sdata.Ks;
        outplot = Sdata.Ms;
        namestr = '%s: M = %d';
        titlestr = 'K = %d';
        filenamestr = '%s_M%d';
    end
    
    incols = ceil(sqrt(length(inplot))); inrows = ceil(length(inplot)/incols);
    
    for iout = 1:length(outplot)
        hfig = figure('Name',sprintf(namestr,filenames{iplot},outplot(iout)),'PaperPositionMode','auto');
        u = get(hfig,'Units'); set(hfig,'Units','inches'); pos = get(hfig,'Position'); c = pos(1)+pos(3)/2; t = pos(2)+pos(4); set(hfig,'Position',[c-1.6,t-2.5,3.2,2.5]); set(hfig,'Units',u);
        haxs = NaN(length(inplot),1);
        for iin = 1:length(inplot)
            haxs(iin) = subplot(inrows,incols,iin); 
            betas1 = col(Sdata.betas_grid(iin,iout,:)).';
            if isequal(scale,'linear')
                boxwidth = 2/3*min(diff(betas1));
            else
                boxwidth = 10^(2/3*min(diff(log10(betas1))));
            end
            errors_range = [floor(log10(min(col(errors_quartiles(:,iin,iout,:))))),ceil(log10(max(col(errors_quartiles(:,iin,iout,:)))))];
            errors_ticks = min(2,errors_range(2)-errors_range(1));
            errors_ticks = errors_range(1)+(0:errors_ticks).*floor((errors_range(2)-errors_range(1))/errors_ticks);
            hs = semilogy(betas1,reshape(errors_use(iin,iout,:),1,length(betas1)),'k');
            set(haxs(iin),'FontSize',8,'YLim',10.^errors_range,'YTick',10.^errors_ticks);
            if iin == incols*(inrows-1)+ceil(incols/2), xlabel(haxs(iin),'Regularization parameter (\beta)'); end
            if iin == 1, ylabel(haxs(iin),'Squared error'); end
            title(haxs(iin),sprintf(titlestr,inplot(iin)));
            if isequal(scale,'linear')
                xlim([min([betas1,boxwidth])-boxwidth,max(betas1)+boxwidth]);
                if boxes
                    line(bsxfun(@plus,betas1,boxwidth/2.*[-1;-1;1;1;-1]),reshape(errors_quartiles([1;end;end;1;1],iin,iout,:),5,length(betas1)),zeros(5,length(betas1)),'LineStyle','-','Color','k','Marker','none');
                end
                set(haxs(iin),'Box','off','XTick',unique(round(betas1),'sorted'));
            else
                set(haxs(iin),'XScale',scale,'XLim',[min(betas1)/boxwidth,max(betas1)*boxwidth]);
                if boxes
                    wp = sqrt(boxwidth); wn = 1/wp;
                    line(bsxfun(@times,betas1,[wn;wn;wp;wp;wn]),reshape(errors_quartiles([1;end;end;1;1],iin,iout,:),5,length(betas1)),zeros(5,length(betas1)),'LineStyle','-','Color','k','Marker','none');
                end
                set(haxs(iin),'Box','off','XTick',10.^unique(round(log10(betas1)),'sorted'));
            end
        end
        drawnow;
        tis = get(haxs,'TightInset'); tis = cat(1,tis{:});
        indscapt = [1;(inrows-1)*incols+ceil(incols/2)];
        ticapt = [tis(indscapt(1),1),tis(indscapt(2),2)];
        tis = [max([tis([1:indscapt(1)-1,indscapt(1)+1:end],1),tis([1:indscapt(2)-1,indscapt(2)+1:end],2)],[],1),max(tis(:,[3,4]),[],1)];
        width = (1-(incols-1)*tis(1)-ticapt(1))/incols-tis(3);
        widths = (width+tis(3)) + [ticapt(1),repmat(tis(1),1,incols-1)];
        height = (1-(inrows-1)*tis(2)-ticapt(2))/inrows-tis(4);
        heights = (height+tis(4)) + [repmat(tis(2),inrows-1,1);ticapt(2)];
        lefts = [ticapt(1),tis(1)+cumsum(widths(1:end-1))];
        bottoms = [tis(2)+flipud(cumsum(heights(end:-1:2)));ticapt(2)];
        for iin = 1:length(inplot)
            newpos = [lefts(mod(iin-1,incols)+1),bottoms(floor((iin-1)/incols)+1),width,height];
            set(haxs(iin),'Position',newpos);
        end
        [~,outfilename] = fileparts(filenames{iplot});
        saveas(hfig,fullfile(outfilepath,[sprintf(filenamestr,outfilename,outplot(iout)),'.fig']));
        print(hfig,fullfile(outfilepath,[sprintf(filenamestr,outfilename,outplot(iout)),'.eps']),'-deps2','-r0');
    end

end

