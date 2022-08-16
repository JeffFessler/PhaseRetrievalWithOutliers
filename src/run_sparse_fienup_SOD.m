% run_sparse_fienup_SOD.m -- Run the L1-modified sparse Fienup
% algorithm on the star of David phantom for Fig. 8 in the paper. See
% run_sparse_fienup_MC.m for more information. See "TODO"s for places that
% require attention before running this code.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

setup_IRT; % code uses Image Reconstruction Toolbox

algotitle = 'L1-Fienup'; % name of algorithm

rng('default');

%% create SOD phantom
[image,dictionary,x,specs] = generate_pattern('SOD.txt',make_disc(10.5,1));
betas = sum(abs(x(:))); % essentially true value

Ninit = 50; % number of initial conditions per trial
I_mults = 10000; % default; will match # from running proposed method below
Nx = numel(x); Nimage = numel(image);
Ms = ceil(Nimage/2);
M_sampling = 'random_plus_zero';
thresh_sf = 1e-4; % stopping criterion
I_ap = 50; % # of alternating projection iterations per initialization

% where to store temporary and final results
% TODO: folder must exist
outfolder = fullfile('..','..','results','phase_retrieval');

% TODO: specify parameters for simulation (examples for figures in paper
% follow)

% AWGN + 1% outliers (Fig. 8)
AWGN_SNR = 60; % dB
outliers_pct = 1; outliers_range = [1,2];
AWLN_SNR = Inf; % dB
outfile = 'SOD_AWGN1pctout_sfl1'; % name of final results file
multsfiles = {'SOD_AWGN1pctout_p1q2.mat'};

% AWGN + 0.1% outliers (Fig. 11(c) in supplement)
AWGN_SNR = 60; % dB
outliers_pct = 0.1; outliers_range = [1,2];
AWLN_SNR = Inf; % dB
outfile = 'SOD_AWGN0p1pctout_sfl1'; % name of final results file
multsfiles = {'SOD_AWGN0p1pctout_p1q2.mat'};

% AWGN + 0.01% outliers (Fig. 11(b) in supplement)
AWGN_SNR = 60; % dB
outliers_pct = 0.01; outliers_range = [1,2];
AWLN_SNR = Inf; % dB
outfile = 'SOD_AWGN0p01pctout_sfl1'; % name of final results file
multsfiles = {'SOD_AWGN0p01pctout_p1q2.mat'};

% AWGN + 0.001% outliers (Fig. 11(a) in supplement)
AWGN_SNR = 60; % dB
outliers_pct = 0.001; outliers_range = [1,2];
AWLN_SNR = Inf; % dB
outfile = 'SOD_AWGN0p001pctout_sfl1'; % name of final results file
multsfiles = {'SOD_AWGN0p001pctout_p1q2.mat'};

%% do experiments
Ms_grid = Ms(:);
[~,betas_grid] = ndgrid(Ms,betas); % set of all (M,beta) pairs to test

samps = false(Nimage,numel(Ms_grid)); % measurement indices (for M < N)
ys = cell(1,numel(Ms_grid)); % measurements for each trial
errors = zeros(numel(Ms_grid),numel(betas)); % error of global minimizer
objectives_true = zeros(numel(Ms_grid),numel(betas)); % obj. function value for true signal
objectives = zeros(numel(Ms_grid),numel(betas)); % obj. function values for global minimizer
xs_prs = zeros(Nx,numel(Ms_grid),numel(betas)); % global minimizer signals
images_best = zeros(Nimage,numel(Ms_grid),numel(betas)); % global minimizer images
pr_opts = cell(numel(Ms_grid),numel(betas)); % reconstruction options
infos = cell(numel(Ms_grid),numel(betas)); % reconstruction info

%% set # of multiplies to match proposed method
I_mults_grid = I_mults.*ones(size(Ms_grid));
if exist('multsfiles','var') && ~isempty(multsfiles)
    for ifile = 1:length(multsfiles)
        try
            Smultsfile = load(fullfile(outfolder,multsfiles{ifile}),'Ms_grid','infos');
            [Smultsfile.Ms_unique,~,Smultsfile.Ms_inds] = unique(Smultsfile.Ms_grid(:),'stable');
            [Smultsfile.my_Ms_unique,~,Smultsfile.my_Ms_inds] = unique(Ms_grid(:),'stable');
            [~,ia,ib] = intersect(Smultsfile.Ms_unique,Smultsfile.my_Ms_unique);
            for ii = 1:length(ia)
                I_mults_grid(Smultsfile.my_Ms_inds == ib(ii)) = max(max(col(I_mults_grid(Smultsfile.my_Ms_inds == ib(ii)))),max(col(arrayfun(@(info) info.nAmults, Smultsfile.infos(:,Smultsfile.Ms_inds == ia(ii))))));
            end
        catch 
            fprintf(2,'Warning: unable to get # of multiplies from file \"%s\".\n',multsfiles{ifile});
        end
        clear Smultsfile ia ib ii;
    end
end

%% do simulation
for i1 = 1:numel(Ms_grid)
    %% setup parameters
    M = Ms_grid(i1);
    outliers = ceil(outliers_pct*M/100);

    %% set up system matrix
    if M <= Nimage
        samp = false(size(image));
        switch lower(M_sampling)
            case 'random_plus_zero'
                samp(1) = true;
                samp(randperm(Nimage-1,M-1)+1) = true;
                samp = fftshift(samp);
            case 'random'
                samp(randperm(Nimage,M)) = true;
            otherwise
                fail('Invalid M_sampling');
        end
        samps(:,i1) = samp(:);
        
        specs_use = specs;
        conv_fft_forw = @(arg,in) (1/sqrt(Nimage)).*specs_use.pattern_fft.*in;
        conv_fft_back = @(arg,in) (1/sqrt(Nimage)).*conj(specs_use.pattern_fft).*in;
        conv_fft_fatrix = fatrix2('imask',samp,'omask',samp,'forw',conv_fft_forw,'back',conv_fft_back);
        A = conv_fft_fatrix*Gdft('samp',samp,'fftshift',1,'ifftshift',1);
        A_curvature = max(abs(specs_use.pattern_fft(:))).^2;
    else
        fail('M > N not supported in 2-D case');
    end

    %% generate signal, samples
    y = sample_gen(A,x(:),'AWGN_SNR',AWGN_SNR,'outliers',outliers,'outliers_range',outliers_range,'AWLN_SNR',AWLN_SNR,'q',2);
    norm_factor = 1; 
    ys{i1} = y;

    for i2 = 1:numel(betas)
        beta = betas_grid(sub2ind([numel(Ms_grid),numel(betas)],i1,i2));
        objectives_true(i1,i2) = sum(col(abs(y-abs(A*x(:)).^2).^2))./norm_factor;
        
        %% perform reconstruction
        % options for phase retrieval algorithm
        pr_ap_opts_use = {'I_ap',I_ap,'thresh_sf',thresh_sf};
        pr_opts_use = {'I_mults_min',I_mults_grid(i1),'Ninit_max',Inf,'Ninit_min',Ninit,'I_mults_max',Inf,'A_unitary',false,'ap_opts',pr_ap_opts_use,'norm_factor',norm_factor,'pcg_opts',{'niter',5}};
        
        pr_opts{i1,i2} = pr_opts_use;
        
        [x_pr,objectives(i1,i2),infos{i1,i2}] = phase_retrieval_sparse_fienup(A,y,beta,pr_opts_use{:});
        x_pr = reshape(x_pr,size(x));
        
        %% evaluate
        xs_prs(:,i1,i2) = x_pr(:);
        [image_best,errors(i1,i2)] = find_best_match_complex(dictionary*x_pr,image);
        images_best(:,i1,i2) = image_best(:);
    end
end

samps = reshape(samps,[Nimage,size(Ms_grid)]); % measurement indices (for M < N)
ys = reshape(ys,size(Ms_grid)); % measurements for each trial
errors = reshape(errors,size(betas_grid)); % error of global minimizer
objectives_true = reshape(objectives_true,size(betas_grid)); % obj. function value for true signal
objectives = reshape(objectives,size(betas_grid)); % obj. function values for global minimizer
xs_prs = reshape(xs_prs,[Nx,size(betas_grid)]); % global minimizer signals
images_best = reshape(images_best,[Nimage,size(betas_grid)]); % global minimizer images
pr_opts = reshape(pr_opts,size(betas_grid)); % reconstruction options
infos = reshape(cat(1,infos{:}),size(betas_grid));

%% save
if ~exist(fullfile(outfolder,sprintf('%s.mat',outfile)),'file')
    save(fullfile(outfolder,sprintf('%s.mat',outfile)),'-v7.3');
else
    ii = 2;
    while exist(fullfile(outfolder,sprintf('%sv%d.mat',outfile,ii)),'file') > 0
        ii = ii + 1;
    end
    save(fullfile(outfolder,sprintf('%sv%d.mat',outfile,ii)),'-v7.3');
end
return;

%% plot
figure; im(abs(reshape(images_best,[size(image),size(betas_grid)])));
