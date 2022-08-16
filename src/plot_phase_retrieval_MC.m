% plot_phase_retrieval_MC.m -- This function plots the reconstruction
% quality trends for the proposed method and/or L1-modified sparse Fienup
% method as a function of sparsity (K) and measurements (M). This code
% loads results files created by run_phase_retrieval_MC.m or
% run_sparse_fienup_MC.m. You should run those first. These figures are
% included in Figs. 4, 5 of the paper (median PSERs) and Figs. 3, 4 of the
% supplement (mean PSERs).
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

%%
hfigs(1) = figure('Name','Median','PaperPositionMode','auto');
hfigs(2) = figure('Name','Mean','PaperPositionMode','auto');
for ifig = 1:2
    u = get(hfigs(ifig),'Units'); set(hfigs(ifig),'Units','inches'); pos = get(hfigs(ifig),'Position'); c = pos(1)+pos(3)/2; t = pos(2)+pos(4); set(hfigs(ifig),'Position',[c-3.25,t-1.6,6.5,1.6]); set(hfigs(ifig),'Units',u);
end
drawnow;
minPSNRs = Inf(1,2); maxPSNRs = -Inf(1,2);
haxs = NaN(Nplots+1,2);
for iplot = 1:Nplots
    Sdata = load(fullfile(filepath,filenames{iplot}),'PSNRs_mean','errors','Ms','Ks','N','algotitle');
    PSNRs_use = cat(3,-10.*log10(reshape(median(Sdata.errors.^2,1),size(Sdata.PSNRs_mean))),Sdata.PSNRs_mean);
    minPSNRs = min([reshape(PSNRs_use,[],2);minPSNRs],[],1);
    maxPSNRs = max([reshape(PSNRs_use,[],2);maxPSNRs],[],1);
    for ifig = 1:2
        set(0,'CurrentFigure',hfigs(ifig));
        haxs(iplot,ifig) = subplot(1,Nplots,iplot); imagesc(log2(Sdata.Ms),Sdata.Ks,PSNRs_use(:,:,ifig)); colormap(gray); 
        set(haxs(iplot,ifig),'FontSize',8,'Box','on','YDir','normal','XTick',log2(fliplr(Sdata.Ms(end:-2:1))),'XTickLabel',arrayfun(@(M) sprintf('%g',M/Sdata.N),fliplr(Sdata.Ms(end:-2:1)),'UniformOutput',false));
        if ~isempty(Sdata.algotitle)
            ht = title(haxs(iplot,ifig),Sdata.algotitle); 
            if strcmpi(Sdata.algotitle,algoproposed), set(ht,'FontWeight','bold'); end
        end
        if (iplot == floor(Nplots/2)+1)
            xlabel(haxs(iplot,ifig),'Measurement fraction (M/N)'); drawnow;
            pos = get(haxs(iplot,ifig),'Position');
            bspace = ceil(pos(2)/0.01)*0.01; height = floor(pos(4)/0.01)*0.01;
        end
        if iplot > 1
            set(haxs(iplot,ifig),'YTick',[]);
        else
            ylabel(haxs(iplot,ifig),'Sparsity fraction (K/N)'); 
            set(haxs(iplot,ifig),'YTick',Sdata.Ks,'YTickLabel',arrayfun(@(K) sprintf('%d/%d',K,Sdata.N),Sdata.Ks,'UniformOutput',false));
            ti = get(haxs(iplot,ifig),'TightInset');
            lspace = ceil((ti(1))/0.01)*0.01;
        end
    end
end

%%
for ifig = 1:2
    set(0,'CurrentFigure',hfigs(ifig));
    for iplot2 = 1:Nplots
        set(haxs(iplot2,ifig),'CLim',[floor(minPSNRs(ifig)/10)*10,round(maxPSNRs(ifig)/10)*10]);
    end
    set(hfigs(ifig),'CurrentAxes',haxs(Nplots,ifig));
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

    width = (1-lspace-rspace)/Nplots-ispace;
    for iplot2 = 1:Nplots
        set(haxs(iplot2,ifig),'Position',[lspace+(iplot2-1)*(width+ispace),bspace,width,height]);
    end
    set(haxs(Nplots+1,ifig),'Position',[lspace+Nplots*(width+ispace),bspace,ispace,height]);
end

%% save
[outfilename,outfilepath] = uiputfile('*','Save MC results to...');
if isequal(outfilename,0) || isequal(outfilepath,0), return; end
[~,outfilename] = fileparts(outfilename);

saveas(hfigs(1),fullfile(outfilepath,[outfilename,'.fig']));
print(hfigs(1),fullfile(outfilepath,[outfilename,'.eps']),'-deps2','-r0');
saveas(hfigs(2),fullfile(outfilepath,[outfilename,'_means.fig']));
print(hfigs(2),fullfile(outfilepath,[outfilename,'_means.eps']),'-deps2','-r0');
