% run_phase_retrieval_MC_noise.m -- Run the proposed method for a range of
% outliers and/or noise levels. See run_phase_retrieval_MC.m for more
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
Ntrials = 50; % # of trials
Ninit = 50; % # of initializations
N = 128; % signal length
Ms = ceil((2.^linspace(-2,0,5))*N); % measurements
M_sampling = 'random_plus_zero';
Ks = [3,5]; % # of sparse points
sparse_min = 0; sparse_max = 1; % magnitude range for nonzeros in signal

% where to store temporary and final results
% TODO: folder must exist
outfolder = fullfile('..','..','results','phase_retrieval');

% TODO: specify parameters for simulation (examples for figures in paper
% follow)

% vary range and # of outliers (Fig. 6)
AWGN_SNRs = 40; % dB
outliers = 2.^(1:4); % absolute #
outlier_variances = [1/24,1/12]; % outliers relative to max measurement
AWLN_SNRs = Inf; % dB
p = 1;
q = 2;
betas_k = 10.^[-0.9,-0.9]; % from run_phase_retrieval_MC.m for these K's
outfile = 'MC_AWGNoutliers_p1q2'; % name of final results file

% varying noise level (Fig. 7)
AWGN_SNRs = 20:10:60; % dB
outliers = 5; % absolute #
outlier_variances = 1/12; % outliers relative to max measurement
AWLN_SNRs = Inf; % dB
p = 1;
q = 2;
betas_k = 10.^[-0.9,-0.9]; % from run_phase_retrieval_MC.m for these K's
outfile = 'MC_AWGNnoise_p1q2'; % name of final results file

%% do experiments

[AWGN_SNRs,AWLN_SNRs,outliers,outlier_variances] = ndgrid(AWGN_SNRs,AWLN_SNRs,outliers,outlier_variances);
AWGN_SNRs = AWGN_SNRs(:); AWLN_SNRs = AWLN_SNRs(:);
outliers = outliers(:); outlier_variances = outlier_variances(:);
[Ks_grid,Ms_grid,AWGN_SNRs_grid] = ndgrid(Ks,Ms,AWGN_SNRs); % set of all (K,M) pairs to test
[betas_grid,~,AWLN_SNRs_grid] = ndgrid(betas_k,Ms,AWLN_SNRs);
[~,~,outliers_grid] = ndgrid(Ks,Ms,outliers);
[~,~,outlier_variances_grid] = ndgrid(Ks,Ms,outlier_variances);

% storage
errors = zeros([Ntrials,size(Ks_grid)]); % error of global minimizer for each trial
objectives_true = zeros([Ntrials,size(Ks_grid)]); % obj. function value for true signal
objectives = zeros([Ntrials,size(Ks_grid)]); % obj. function values for global minimizer
xs_true = zeros([N,Ntrials,size(Ks_grid)]); % true signals for each trial
samps = false([N,Ntrials,size(Ks_grid)]); % measurement indices for each trial (for M < N)
ys = cell(size(Ks_grid)); % measurements for each trial
xs_best = zeros([N,Ntrials,size(Ks_grid)]); % global minimizer signals for each trial
infos = cell([Ntrials,size(Ks_grid)]); % reconstruction info

parpool; % for parallel processing

for i12 = 1:numel(Ks_grid)
    try
        % load existing temporary results if re-run after crash (note: no
        % validation on loaded results)
        load(fullfile(outfolder,sprintf('%s_part%03d.mat',outfile,i12)));
        loaded = true;
    catch eload
        loaded = false;
    end

    if ~loaded
        %% setup parameters
        K = Ks_grid(i12);
        M = Ms_grid(i12);
        AWGN_SNR = AWGN_SNRs_grid(i12);
        AWLN_SNR = AWLN_SNRs_grid(i12);
        outlier = outliers_grid(i12);
        outlier_variance = outlier_variances_grid(i12);
        outlier_range = 1+[0,sqrt(12*outlier_variance)];
        beta = betas_grid(i12);        

        if p == 1 && q == 2
            pr_ADMM_opts = {'mu_u',1,'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',100,'A_curvature',1,'A_unitary',(M >= N)};
            pr_ox_opts = {'eps_ox',1e-10,'I_ox',10,'ADMM_opts',pr_ADMM_opts};
            pr_opts = {'Ninit_min',Ninit*(1 + (M >= N)),'ox_opts',pr_ox_opts};
        elseif p == 2 && q == 2
            pr_ADMM_opts = {'mu_u',1,'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',100,'A_curvature',1,'A_unitary',(M >= N)};
            pr_ox_opts = {'eps_ox',1e-10,'I_ox',10,'ADMM_opts',pr_ADMM_opts};
            pr_opts = {'Ninit_min',Ninit*(1 + (M >= N)),'ox_opts',pr_ox_opts};
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
        end

        % storage for this (K,M) pair
        errors1 = zeros(Ntrials,1);
        objectives_true1 = zeros(Ntrials,1);
        objectives1 = zeros(Ntrials,1);
        xs_true1 = zeros(N,Ntrials);
        samps1 = false(N,Ntrials);
        ys1 = zeros(M,Ntrials);
        xs_best1 = zeros(N,Ntrials);
        infos1 = cell(Ntrials,1);

        parfor itrial = 1:Ntrials
            %% set up system matrix
            if M < N
                % undersampled case
                samp = false(N,1);
                switch lower(M_sampling)
                    case 'random_plus_zero'
                        samp(1) = true;
                        samp(randperm(N-1,M-1)+1) = true;
                    case 'random'
                        samp(randperm(N,M)) = true;
                    case 'low-freq'
                        samp(1:ceil(M/2)) = true;
                        samp(end-floor(M/2)+1:end) = true;
                    otherwise
                        fail('Invalid M_sampling');
                end
                samps1(:,itrial) = samp;
                A = (1/sqrt(N))*myDFT(N,N,samp);
            else
                % critical/oversampled case
                A = (1/sqrt(M))*myDFT(N,M);
            end

            %% generate signal, samples
            x = signal_gen(N,K,sparse_min,sparse_max);
            xs_true1(:,itrial) = x;
            y = sample_gen(A,x,'AWGN_SNR',AWGN_SNR,'outliers',outlier,'outliers_range',outlier_range,'AWLN_SNR',AWLN_SNR,'q',q);
            ys1(:,itrial) = y;
            norm_factor = compute_beta_normalization(y,outlier,outlier_range,AWGN_SNR,AWLN_SNR,p);
            objectives_true1(itrial) = sum(col(abs(y-abs(A*x).^q).^p))./norm_factor + beta*norm(col(x),1);

            %% perform reconstruction
            pr_opts_use = pr_opts;
            pr_opts_use = [pr_opts_use,{'norm_factor',norm_factor}];
            if p == 1 && q == 2
                [x_pr,objectives1(itrial),infos1{itrial}] = phase_retrieval_p1_q2(size(x),A,y,beta,pr_opts_use{:});
            elseif p == 2 && q == 2
                [x_pr,objectives1(itrial),infos1{itrial}] = phase_retrieval_p2_q2(size(x),A,y,beta,pr_opts_use{:});
            else
                fail('Unknown p = %d, q = %d combination!',p,q);
                x_pr = []; % should not get here
            end

            %% evaluate
            [x_pr,errors1(itrial)] = find_best_match_complex(x_pr,x);
            xs_best1(:,itrial) = x_pr;
        end
        
        %% write temporary results
        try
            save(fullfile(outfolder,sprintf('%s_part%03d.mat',outfile,i12)),'-v7.3','errors1','objectives_true1','objectives1','xs_true1','samps1','ys1','xs_best1','K','M','beta','pr_opts','infos1','AWGN_SNR','AWLN_SNR','outlier','outlier_variance','outlier_range');
        catch esave
            fprintf(2,'Error saving ''%s_part%03d.mat'': %s\n',outfile,i12,esave.message);
        end
    end

    %% update storage
    errors(:,i12) = errors1;
    objectives_true(:,i12) = objectives_true1;
    objectives(:,i12) = objectives1;
    xs_true(:,:,i12) = xs_true1;
    samps(:,:,i12) = samps1;
    ys{i12} = ys1;
    xs_best(:,:,i12) = xs_best1;
    infos(:,i12) = infos1;

    %% report findings
    fprintf(1,'[K = %d, M = %d, out = (%g,%g,%d,%g)]: Error: %g\n',K,M,AWGN_SNR,AWLN_SNR,outlier,outlier_variance,sqrt(mean(errors1.^2)));

    clear errors1 objectives_true1 objectives1 infos1 xs_true1 samps1 ys1 xs_best1;
end

% matlabpool('close');
delete(gcp('nocreate'));

infos = reshape(cat(1,infos{:}),[Ntrials,size(Ks_grid)]);

%% final reporting
MSEs_mean = reshape(mean(errors.^2,1),size(Ks_grid));
PSNRs_mean = -10.*log10(MSEs_mean);

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
delete(fullfile(outfolder,sprintf('%s_part*.mat',outfile)));
