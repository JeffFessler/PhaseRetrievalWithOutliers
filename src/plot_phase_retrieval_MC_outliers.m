% plot_phase_retrieval_MC_outliers.m -- This function plots the
% reconstruction quality trends for the proposed method and/or L1-modified
% sparse Fienup method as a function of outliers and measurements (M). This
% code loads results files created by run_phase_retrieval_MC_noise.m or
% run_sparse_fienup_MC_noise.m. You should run those first. These figures
% are included in Fig. 6 of the paper (median PSERs) and Figs. 5 (median),
% 7, 8 (mean) of the supplement.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

filepath = fullfile('..','..','results','phase_retrieval');
fileprefix = 'MC_';
filesuffix = '.mat';
algoproposed = 'Proposed'; % highlight this one

[filenames,filepath] = uigetfile({[fileprefix '*' filesuffix];['*' filesuffix]},'Plot phase retrieval MC...',filepath,'MultiSelect','on');
if isempty(filenames) || isequal(filenames,0) || isequal(filepath,0), return; end
if ~iscell(filenames), filenames = {filenames}; end

Nplots = length(filenames);

if Nplots == 0, return; end
load(fullfile(filepath,filenames{1}),'Ks','AWGN_SNRs','AWLN_SNRs','outlier_variances');
AWGN_SNRs = unique(AWGN_SNRs,'stable');
AWLN_SNRs = unique(AWLN_SNRs,'stable');
outlier_variances = unique(outlier_variances,'stable');
Nfigs = length(Ks)*length(AWGN_SNRs)*length(AWLN_SNRs)*length(outlier_variances);
[Ks_grid,AWGN_SNRs_grid,AWLN_SNRs_grid,outlier_variances_grid] = ndgrid(Ks,AWGN_SNRs,AWLN_SNRs,outlier_variances);

%% 

hfigs = zeros(Nfigs,2);
for ifig = 1:Nfigs
    hfigs(ifig,1) = figure('Name',sprintf('Median: K = %d, AWGN = %g, AWLN = %g, Outliers Range = [%g,%g]',Ks_grid(ifig),AWGN_SNRs_grid(ifig),AWLN_SNRs_grid(ifig),1,1+sqrt(12*outlier_variances_grid(ifig))),'PaperPositionMode','auto');
    u = get(hfigs(ifig,1),'Units'); set(hfigs(ifig,1),'Units','inches'); pos = get(hfigs(ifig,1),'Position'); c = pos(1)+pos(3)/2; t = pos(2)+pos(4); set(hfigs(ifig,1),'Position',[c-3.25,t-1.6,6.5,1.6]); set(hfigs(ifig,1),'Units',u);
    hfigs(ifig,2) = figure('Name',sprintf('Mean: K = %d, AWGN = %g, AWLN = %g, Outliers Range = [%g,%g]',Ks_grid(ifig),AWGN_SNRs_grid(ifig),AWLN_SNRs_grid(ifig),1,1+sqrt(12*outlier_variances_grid(ifig))),'PaperPositionMode','auto');
    u = get(hfigs(ifig,2),'Units'); set(hfigs(ifig,2),'Units','inches'); pos = get(hfigs(ifig,2),'Position'); c = pos(1)+pos(3)/2; t = pos(2)+pos(4); set(hfigs(ifig,2),'Position',[c-3.25,t-1.6,6.5,1.6]); set(hfigs(ifig,2),'Units',u);
end
minPSNRs = Inf(1,Nfigs,2); maxPSNRs = -Inf(1,Nfigs,2);
haxs = NaN(Nplots+1,Nfigs,2);
for iplot = 1:Nplots
    Sdata = load(fullfile(filepath,filenames{iplot}),'PSNRs_mean','errors','Ms','Ks','N','outliers','algotitle');
    Sdata.outliers = unique(Sdata.outliers,'stable');
    PSNRs_use = cat(4,-10.*log10(reshape(median(Sdata.errors.^2,1),size(Sdata.PSNRs_mean))),Sdata.PSNRs_mean);
    PSNRs_use = reshape(PSNRs_use,[length(Sdata.Ks),length(Sdata.Ms),length(AWGN_SNRs),length(AWLN_SNRs),length(Sdata.outliers),length(outlier_variances),2]);
    PSNRs_use = reshape(permute(PSNRs_use,[2,5,1,3,4,6,7]),[length(Sdata.Ms),length(Sdata.outliers),Nfigs,2]);
    minPSNRs = min([reshape(PSNRs_use,[],Nfigs,2);minPSNRs],[],1);
    maxPSNRs = max([reshape(PSNRs_use,[],Nfigs,2);maxPSNRs],[],1);
    for ifig = 1:2*Nfigs
        set(0,'CurrentFigure',hfigs(ifig));
        PSNRs_use1 = PSNRs_use(:,:,ifig);

        haxs(iplot,ifig) = subplot(1,Nplots,iplot);
        imagesc(log2(Sdata.outliers),log2(Sdata.Ms),PSNRs_use1); colormap(gray);
        set(haxs(iplot,ifig),'FontSize',8,'Box','on','YDir','reverse');
        set(haxs(iplot,ifig),'XTick',log2(Sdata.outliers),'XTickLabel',arrayfun(@(outlier) sprintf('%d',outlier),Sdata.outliers,'UniformOutput',false));
        if ~isempty(Sdata.algotitle)
            ht = title(haxs(iplot,ifig),Sdata.algotitle); 
            if strcmpi(Sdata.algotitle,algoproposed), set(ht,'FontWeight','bold'); end
        end
        if iplot == floor(Nplots/2)+1
            xlabel(haxs(iplot,ifig),'Outliers'); drawnow;
        end
        if iplot > 1
            set(haxs(iplot,ifig),'YTick',[]);
        else
            ylabel(haxs(iplot,ifig),'Measurement fraction (M/N)'); 
            set(haxs(iplot,ifig),'YTick',log2(fliplr(Sdata.Ms(end:-2:1))),'YTickLabel',arrayfun(@(M) sprintf('%g',M/Sdata.N),fliplr(Sdata.Ms(end:-2:1)),'UniformOutput',false));
            ti = get(haxs(iplot,ifig),'TightInset');
            lspace = ceil((ti(1))/0.01)*0.01;
        end
    end
end

for ifig = 1:Nfigs*2
    for iplot2 = 1:Nplots
        set(haxs(iplot2,ifig),'CLim',[floor(minPSNRs(ifig)/10)*10,round(maxPSNRs(ifig)/10)*10]);
    end
    
    set(0,'CurrentFigure',hfigs(ifig));
    haxs(Nplots+1,ifig) = colorbar('FontSize',8); drawnow;
    set(haxs(Nplots+1,ifig),'YTickMode','manual','YTickLabelMode','manual'); drawnow;
    tlabels = get(haxs(Nplots+1,ifig),'YTickLabel');
    if ~iscell(tlabels), tlabels = num2cell(tlabels,2); end
    tlabels{end} = [tlabels{end},' dB'];
    set(haxs(Nplots+1,ifig),'YTickLabel',tlabels);
    pos = get(haxs(Nplots+1,ifig),'Position'); 
    ispace = ceil(pos(3)/0.01)*0.01;
    try
        ti = get(haxs(Nplots+1,ifig),'TightInset');
    catch
        ti = [0,0,0.05,0];
    end
    rspace = ceil((ispace+ti(3))/0.01)*0.01;
    
    pos = get(haxs(floor(Nplots/2)+1,ifig),'Position');
    bspace = ceil(pos(2)/0.01)*0.01; height = floor(pos(4)/0.01)*0.01;

    width = (1-lspace-rspace)/Nplots-ispace;
    for iplot2 = 1:Nplots
        set(haxs(iplot2,ifig),'Position',[lspace+(iplot2-1)*(width+ispace),bspace,width,height]);
    end
    set(haxs(Nplots+1,ifig),'Position',[lspace+Nplots*(width+ispace),bspace,ispace,height]);
end

%% save
[outfilename,outfilepath] = uiputfile('*','Save MC outliers results to...');
if isequal(outfilename,0) || isequal(outfilepath,0), return; end
[~,outfilename] = fileparts(outfilename);

for ifig = 1:Nfigs
    if length(Ks) > 1, identstr = sprintf('_K%d',Ks_grid(ifig)); else identstr = ''; end
    if length(AWGN_SNRs) > 1, identstr = [identstr,sprintf('_AWGN%2g',AWGN_SNRs_grid(ifig))]; end
    if length(AWLN_SNRs) > 1, identstr = [identstr,sprintf('_AWLN%2g',AWLN_SNRs_grid(ifig))]; end
    if length(outlier_variances) > 1, identstr = [identstr,sprintf('_outmax%.2g',1+sqrt(12*outlier_variances_grid(ifig)))]; end
    saveas(hfigs(ifig,1),fullfile(outfilepath,[outfilename,identstr,'.fig']));
    print(hfigs(ifig,1),fullfile(outfilepath,[outfilename,identstr,'.eps']),'-deps2','-r0');
    saveas(hfigs(ifig,2),fullfile(outfilepath,[outfilename,identstr,'_means.fig']));
    print(hfigs(ifig,2),fullfile(outfilepath,[outfilename,identstr,'_means.eps']),'-deps2','-r0');
end

