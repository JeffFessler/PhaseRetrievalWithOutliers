% run_phase_retrieval_SOD.m -- Run the proposed method on the star of David
% phantom for Fig. 8 in the paper. See run_phase_retrieval_MC.m for more
% information. See "TODO"s for places that require attention before running
% this code.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

setup_IRT; % code uses Image Reconstruction Toolbox

algotitle = 'Proposed'; % name of algorithm

rng('default');

%% create SOD
[image,dictionary,x,specs] = generate_pattern('SOD.txt',make_disc(10.5,1));

Ninit = 50; % number of initial conditions per trial
Nx = numel(x); Nimage = numel(image); % size of coefficients, image
Ms = ceil(Nimage/2); % measurements
M_sampling = 'random_plus_zero';

% where to store temporary and final results
% TODO: folder must exist
outfolder = fullfile('..','..','results','phase_retrieval');

% TODO: specify parameters for simulation (examples for figures in paper
% follow)

% AWGN + 1% outliers (Fig. 8)
betas = linspace(0.1,0.4,7); % regularization parameters to try
p = 1; q = 2;
AWGN_SNR = 60;
outliers_pct = 1; outliers_range = [1,2];
AWLN_SNR = Inf;
outfile = 'SOD_AWGN1pctout_p1q2'; % name of final results file

% AWGN + 0.1% outliers (Fig. 11(c) in supplement)
betas = 0.3; % regularization parameters to try
p = 1; q = 2;
AWGN_SNR = 60;
outliers_pct = 0.1; outliers_range = [1,2];
AWLN_SNR = Inf;
outfile = 'SOD_AWGN0p1pctout_p1q2'; % name of final results file

% AWGN + 0.01% outliers (Fig. 11(b) in supplement)
betas = 0.3; % regularization parameters to try
p = 1; q = 2;
AWGN_SNR = 60;
outliers_pct = 0.01; outliers_range = [1,2];
AWLN_SNR = Inf;
outfile = 'SOD_AWGN0p01pctout_p1q2'; % name of final results file

% AWGN + 0.001% outliers (Fig. 11(a) in supplement)
betas = 0.3; % regularization parameters to try
p = 1; q = 2;
AWGN_SNR = 60;
outliers_pct = 0.001; outliers_range = [1,2];
AWLN_SNR = Inf;
outfile = 'SOD_AWGN0p001pctout_p1q2'; % name of final results file

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
    y = sample_gen(A,x(:),'AWGN_SNR',AWGN_SNR,'outliers',outliers,'outliers_range',outliers_range,'AWLN_SNR',AWLN_SNR,'q',q);
    norm_factor = 1;
    ys{i1} = y;

    for i2 = 1:numel(betas)
        beta = betas_grid(sub2ind([numel(Ms_grid),numel(betas)],i1,i2));
        objectives_true(i1,i2) = sum(col(abs(y-abs(A*x(:)).^q).^p))./norm_factor + beta*norm(col(x),1);
        
        %% perform reconstruction
        % options for phase retrieval algorithm
        if p == 1 && q == 2
            pr_ADMM_opts_use = {'mu_u',1,'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',100,'A_curvature',A_curvature,'A_unitary',false};
            pr_ox_opts_use = {'eps_ox',1e-10,'I_ox',20,'ADMM_opts',pr_ADMM_opts_use};
            pr_opts_use = {'Ninit_min',Ninit,'ox_opts',pr_ox_opts_use,'norm_factor',norm_factor};
        elseif p == 2 && q == 2
            pr_ADMM_opts_use = {'mu_u',1,'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',100,'A_curvature',A_curvature,'A_unitary',false};
            pr_ox_opts_use = {'eps_ox',1e-10,'I_ox',20,'ADMM_opts',pr_ADMM_opts_use};
            pr_opts_use = {'Ninit_min',Ninit,'ox_opts',pr_ox_opts_use,'norm_factor',norm_factor};
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
            pr_opts_use = [];
        end
        pr_opts{i1,i2} = pr_opts_use;
        
        if p == 1 && q == 2
            [x_pr,objectives(i1,i2),infos{i1,i2}] = phase_retrieval_p1_q2([numel(x),1],A,y,beta,pr_opts_use{:});
            x_pr = reshape(x_pr,size(x));
        elseif p == 2 && q == 2
            [x_pr,objectives(i1,i2),infos{i1,i2}] = phase_retrieval_p2_q2([numel(x),1],A,y,beta,pr_opts_use{:});
            x_pr = reshape(x_pr,size(x));
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
            x_pr = []; % should not get here
        end
        
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
figure; im(abs(reshape(images_best,[size(image),size(betas_grid)])),[0,1]);
