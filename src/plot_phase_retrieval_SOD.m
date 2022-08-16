% plot_phase_retrieval_SOD.m -- This function plots the images
% reconstructed by the proposed method and/or competing methods and reports
% the PSER values (in dB) for those reconstructions. This function loads
% results files created by run_phase_retrieval_SOD.m or
% run_sparse_fienup_SOD.m. You should run those first. These figures are
% included in Fig. 8 of the paper and Fig. 11 of the supplement.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

setup_IRT; % needs Image Reconstruction Toolbox

% files to plot
filepath = fullfile('..','..','results','phase_retrieval');
fileprefix = 'SOD_';
filesuffix = '.mat';
algoproposed = 'Proposed'; % highlight this one
plotoriginal = false; plotoriginal = true; % if true, include plot of original image
square = false; square = true; % makes figure square if true, otherwise makes in a row

[filenames,filepath] = uigetfile({[fileprefix '*' filesuffix];['*' filesuffix]},'Plot phase retrieval images...',filepath,'MultiSelect','on');
if isempty(filenames) || isequal(filenames,0) || isequal(filepath,0), return; end
if ~iscell(filenames), filenames = {filenames}; end

Nplots = length(filenames);

if Nplots == 0, return; end
load(fullfile(filepath,filenames{1}),'Ms','Nimage');
Nfigs = length(Ms);
if square
    nrows = ceil(sqrt(Nplots+plotoriginal)); ncols = ceil((Nplots+plotoriginal)/nrows);
    figwidth = 3.2;
else
    nrows = 1; ncols = Nplots+plotoriginal;
    figwidth = 6.5;
end

%% 

hfigs = zeros(Nfigs,1);
for ifig = 1:Nfigs
    figheight = nrows*figwidth/ncols;
    hfigs(ifig) = figure('Name',sprintf('M/N = %g',Ms(ifig)/Nimage),'PaperPositionMode','auto');
    u = get(hfigs(ifig),'Units'); set(hfigs(ifig),'Units','inches'); pos = get(hfigs(ifig),'Position'); c = pos(1)+pos(3)/2; t = pos(2)+pos(4); set(hfigs(ifig),'Position',[c-figwidth/2,t-figheight,figwidth,figheight]); set(hfigs(ifig),'Units',u);
end
minPSNRs = Inf(1,Nfigs); maxPSNRs = -Inf(1,Nfigs);
haxs = NaN(Nplots+plotoriginal,Nfigs);
for iplot = 1:Nplots
    Sdata = load(fullfile(filepath,filenames{iplot}));
    
    for ifig = 1:Nfigs
        set(0,'CurrentFigure',hfigs(ifig));
        [error_best,ibest] = min(Sdata.errors(ifig,:));
        if isfield(Sdata,'betas')
            str_best = {sprintf('\\beta = %g',Sdata.betas(ibest)),sprintf('PSER = %.3g dB',20*log10(max(abs(Sdata.image(:)))/error_best))};
        elseif isfield(Sdata,'Ks')
            str_best = {sprintf('K = %d',Sdata.Ks(ibest)),sprintf('PSER = %.3g dB',20*log10(max(abs(Sdata.image(:)))/error_best))};
        else
            str_best = sprintf('PSER = %.3g dB',20*log10(max(abs(Sdata.image(:)))/error_best));
        end
        image_best = reshape(Sdata.images_best(:,ifig,ibest),size(Sdata.image));

        if iplot == 1 && plotoriginal
            haxs(1,ifig) = subplot(nrows,ncols,1);
            im('notick',abs(Sdata.image),[0,max(abs(Sdata.image(:)))]); colormap(gray);
            set(haxs(1,ifig),'FontSize',8,'Box','off','XTick',[],'YTick',[],'Color','none','Visible','off');
            title(haxs(1,ifig),'True Image','Visible','on');
        end
        haxs(iplot+plotoriginal,ifig) = subplot(nrows,ncols,iplot+plotoriginal);
        im('notick',abs(image_best),[0,max(abs(Sdata.image(:)))]); colormap(gray);
        set(haxs(iplot+plotoriginal,ifig),'FontSize',8,'Box','off','XTick',[],'YTick',[],'Color','none','Visible','off');
        if ~isempty(Sdata.algotitle)
            ht = title(haxs(iplot+plotoriginal,ifig),Sdata.algotitle,'Visible','on'); 
            if strcmpi(Sdata.algotitle,algoproposed), set(ht,'FontWeight','bold'); end
        end
        xlabel(haxs(iplot+plotoriginal,ifig),str_best,'Visible','on');
    end
end
drawnow;

%%
for ifig = 1:Nfigs
    % tighten position around image
    u = get(haxs(:,ifig),'Units'); set(haxs(:,ifig),'Units','pixels');
    pos = get(haxs(:,ifig),'Position'); pos = cat(1,pos{:});
    pbas = get(haxs(:,ifig),'PlotBoxAspectRatio'); pbas = cat(1,pbas{:});
    ratios = [pos(:,3)./pbas(:,1),pos(:,4)./pbas(:,2)];
    [~,imin] = min(ratios,[],2);
    fixwidths = imin == 2;
    if any(fixwidths)
        newwidth = pbas(fixwidths,1).*ratios(fixwidths,2);
        pos(fixwidths,1) = pos(fixwidths,1) + (pos(fixwidths,3)-newwidth)./2;
        pos(fixwidths,3) = newwidth;
    end
    if any(~fixwidths)
        newheight = pbas(~fixwidths,2).*ratios(~fixwidths,1);
        pos(~fixwidths,2) = pos(~fixwidths,2) + (pos(~fixwidths,4)-newheight)./2;
        pos(~fixwidths,4) = newheight;
    end
    set(haxs(:,ifig),{'Position'},num2cell(pos,2)); set(haxs(:,ifig),{'Units'},u);
    
    tis = get(haxs(:,ifig),'TightInset'); tis = cat(1,tis{:});
    tis(end+1:nrows*ncols,:) = 0; tis = reshape(tis,[ncols,nrows,4]);
    tis_lr = reshape(max(tis(:,:,[1,3]),[],2),[],2);
    tis_bt = reshape(max(tis(:,:,[2,4]),[],1),[],2);
    
    width = (1 - sum(tis_lr(:)))/ncols;
    height = (1-sum(tis_bt(:)))/nrows;
    lefts = [tis_lr(1,1);tis_lr(2:end,1)+cumsum(tis_lr(1:end-1,1)+tis_lr(1:end-1,2)+width)];
    bottoms = [tis_bt(1:end-1,1)+flipud(cumsum(tis_bt(end:-1:2,1)+tis_bt(end:-1:2,2)+height));tis_bt(end,1)];
    for iplot = 1:Nplots+plotoriginal
        set(haxs(iplot,ifig),'Position',[lefts(mod(iplot-1,ncols)+1),bottoms(floor((iplot-1)/ncols)+1),width,height]);
    end
end

%% save
[outfilename,outfilepath] = uiputfile('*','Save image results to...');
if isequal(outfilename,0) || isequal(outfilepath,0), return; end
[~,outfilename] = fileparts(outfilename);

for ifig = 1:Nfigs
    if length(Ms) > 1, identstr = sprintf('_M%d',Ms(ifig)); else identstr = ''; end
    saveas(hfigs(ifig),fullfile(outfilepath,[outfilename,identstr,'.fig']));
    print(hfigs(ifig),fullfile(outfilepath,[outfilename,identstr,'.eps']),'-deps2','-r0');
end

