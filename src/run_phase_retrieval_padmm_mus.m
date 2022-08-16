% run_phase_retrieval_padmm_mus.m -- Run the proposed method for a range of
% measurements or sparsities, and ADMM penalty parameters (beta). This code
% tests convergence of one run of ADMM/PADMM (not the complete phase
% retrieval) adaptive penalty parameters against fixed values. See "TODO"s
% for places that require attention before running this code.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

setup_IRT; % code uses Image Reconstruction Toolbox

algotitle = 'Proposed'; % name of algorithm

rng('default');
N = 128; % signal length
Ms = ceil([0.5,1]*N); % undersampled and critically sampled cases
M_sampling = 'random_plus_zero';
Ks = [6,8]; % # of sparse points
sparse_min = 0; sparse_max = 1; % magnitude range for nonzeros in signal
ADMMiters = 100; % iterations of ADMM
ADMMoveriterfactor = 4; % run more to get a "true" solution
mus_range = logspace(-2,2,5); % range of fixed mus to test
adaptive_mus_range = [1,0.1]; % range of initial mus for the adaptive method to test

% where to store temporary and final results
% TODO: folder must exist
outfolder = fullfile('..','..','results','phase_retrieval');

% TODO: specify parameters for simulation (examples for figures in paper
% follow)

% AWGN + 5 outliers (Fig. 2)
AWGN_SNR = 40; % dB
outliers = 5; % absolute #
outliers_range = [1,2]; % outliers relative to max measurement
AWLN_SNR = Inf; % dB
p = 1;
q = 2;
betas_k = [0.2,0.25];
outfile = 'mus_AWGN5out_padmm'; % name of final results file

%% do experiments

[Ks_grid,Ms_grid] = ndgrid(Ks,Ms); % set of all (K,M) pairs to test
betas_grid = ndgrid(betas_k,Ms);

% storage
Nmus = length(mus_range)+length(adaptive_mus_range);
errors = zeros([Nmus,size(Ks_grid)]); % error of global minimizer for each trial
objectives_true = zeros(size(Ks_grid)); % obj. function value for true signal
objectives = zeros([Nmus,size(Ks_grid)]); % obj. function values for all initializers
xs_true = zeros([N,size(Ks_grid)]); % true signals for each trial
samps = false([N,size(Ks_grid)]); % measurement indices for each trial (for M < N)
ys = cell(size(Ks_grid)); % measurements for each trial
xs_best = zeros([N,Nmus,size(Ks_grid)]); % global minimizer signals for each trial
info_prs = cell([Nmus,size(Ks_grid)]);

for i12 = 1:numel(Ks_grid)
    %% setup parameters
    K = Ks_grid(i12);
    M = Ms_grid(i12);
    beta = betas_grid(i12);
    
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
        samps(:,i12) = samp;
        A = (1/sqrt(N))*myDFT(N,N,samp);
    else
        % critical/oversampled case
        A = (1/sqrt(M))*myDFT(N,M);
    end
    
    %% generate signal, samples
    x = signal_gen(N,K,sparse_min,sparse_max);
    xs_true(:,i12) = x;
    y = sample_gen(A,x,'AWGN_SNR',AWGN_SNR,'outliers',outliers,'outliers_range',outliers_range,'AWLN_SNR',AWLN_SNR,'q',q);
    ys{i12} = y;
    objectives_true(i12) = sum(col(abs(y-abs(A*x).^q).^p)) + beta*norm(col(x),1);
    
    % do reconstructions (all start from same s value)
    ss_init_use = {Wirt_init(y,A)};

    for imu = 1:length(mus_range)
        if p == 1 && q == 2
            pr_ADMM_opts = {'mu_u',mus_range(imu),'mu_I_adapt',0,'eps_ADMM',1e-10,'I_ADMM',ADMMoveriterfactor*ADMMiters,'A_curvature',1,'A_unitary',(M >= N),'I_ADMM_min',Inf,'verbose',true};
            pr_ox_opts = {'eps_ox',1e-10,'I_ox',1,'ADMM_opts',pr_ADMM_opts};
            pr_opts = {'Ninit_min',1,'ox_opts',pr_ox_opts};
        elseif p == 2 && q == 2
            pr_ADMM_opts = {'mu_u',mus_range(imu),'mu_I_adapt',0,'eps_ADMM',1e-10,'I_ADMM',ADMMoveriterfactor*ADMMiters,'A_curvature',1,'A_unitary',(M >= N),'I_ADMM_min',Inf,'verbose',true};
            pr_ox_opts = {'eps_ox',1e-10,'I_ox',1,'ADMM_opts',pr_ADMM_opts};
            pr_opts = {'Ninit_min',1,'ox_opts',pr_ox_opts};
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
        end
        
        %% perform reconstruction
        pr_opts_use = [pr_opts,{'ss_init',ss_init_use}];
        if p == 1 && q == 2
            [x_pr,objectives(imu,i12),info_prs{imu,i12}] = phase_retrieval_p1_q2(size(x),A,y,beta,pr_opts_use{:});
        elseif p == 2 && q == 2
            [x_pr,objectives(imu,i12),info_prs{imu,i12}] = phase_retrieval_p2_q2(size(x),A,y,beta,pr_opts_use{:});
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
            x_pr = []; % should not get here
        end
        
        %% evaluate
        [x_pr,errors(imu,i12)] = find_best_match_complex(x_pr,x);
        xs_best(:,imu,i12) = x_pr;
        
    end
    
    for imu = 1:length(adaptive_mus_range)
        if p == 1 && q == 2
            pr_ADMM_opts_adapt = {'mu_u',adaptive_mus_range(imu),'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',ADMMoveriterfactor*ADMMiters,'A_curvature',1,'A_unitary',(M >= N),'I_ADMM_min',Inf,'verbose',true};
            pr_ox_opts_adapt = {'eps_ox',1e-10,'I_ox',1,'ADMM_opts',pr_ADMM_opts_adapt};
            pr_opts_adapt = {'Ninit_min',1,'ox_opts',pr_ox_opts_adapt};
        elseif p == 2 && q == 2
            pr_ADMM_opts_adapt = {'mu_u',adaptive_mus_range(imu),'mu_I_adapt',10,'eps_ADMM',1e-10,'I_ADMM',ADMMoveriterfactor*ADMMiters,'A_curvature',1,'A_unitary',(M >= N),'I_ADMM_min',Inf,'verbose',true};
            pr_ox_opts_adapt = {'eps_ox',1e-10,'I_ox',1,'ADMM_opts',pr_ADMM_opts_adapt};
            pr_opts_adapt = {'Ninit_min',1,'ox_opts',pr_ox_opts_adapt};
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
        end
        
        %% perform reconstruction
        pr_opts_use = [pr_opts_adapt,{'ss_init',ss_init_use}];
        if p == 1 && q == 2
            [x_pr,objectives(length(mus_range)+imu,i12),info_prs{length(mus_range)+imu,i12}] = phase_retrieval_p1_q2(size(x),A,y,beta,pr_opts_use{:});
        elseif p == 2 && q == 2
            [x_pr,objectives(length(mus_range)+imu,i12),info_prs{length(mus_range)+imu,i12}] = phase_retrieval_p2_q2(size(x),A,y,beta,pr_opts_use{:});
        else
            fail('Unknown p = %d, q = %d combination!',p,q);
            x_pr = []; % should not get here
        end
        
        %% evaluate
        [x_pr,errors(length(mus_range)+imu,i12)] = find_best_match_complex(x_pr,x);
        xs_best(:,length(mus_range)+imu,i12) = x_pr;
        
    end
    
    %% report findings
    fprintf(1,'[K = %d, M = %d]: Error: %g\n',K,M,min(errors(:,i12)));
end

info_prs = reshape(cat(1,info_prs{:}),[Nmus,size(Ks_grid)]);

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

%% postprocessing to create figures
for i12 = 1:numel(Ks_grid)
    figure;
    niters_uses = arrayfun(@(info_pr) info_pr.ox_info(1).ADMM_info(1).niters/ADMMoveriterfactor,info_prs(:,i12));
    
    [Phi_bests] = arrayfun(@(info_pr) min(info_pr.ox_info(1).ADMM_info(1).Phis(1:info_pr.ox_info(1).ADMM_info(1).niters+1)),info_prs(:,i12));
    Phi_best = min(Phi_bests);

    hs = zeros(Nmus,1);
    hs(1) = semilogy((0:niters_uses(1)).',info_prs(1,i12).ox_info(1).ADMM_info(1).Phis(1:niters_uses(1)+1)-Phi_best); hold all;
    hs(2:Nmus) = arrayfun(@(niters_use,info_pr) semilogy((0:niters_use).',info_pr.ox_info(1).ADMM_info(1).Phis(1:niters_use+1)-Phi_best),niters_uses(2:end),info_prs(2:end,i12));

    Phi_labels = arrayfun(@(niters_use,info_pr) min(info_pr.ox_info(1).ADMM_info(1).Phis(1:niters_use+1)),niters_uses,info_prs(:,i12));
    [~,icurveorder] = sort(Phi_labels,'descend');
    mus_strings = [arrayfun(@(mu) sprintf('\\mu = %g',mu),mus_range,'UniformOutput',false),arrayfun(@(mu) sprintf('Adaptive, \\mu = %g',mu),adaptive_mus_range,'UniformOutput',false)].';
    legend(hs(icurveorder),mus_strings(icurveorder),'Location','Best');
    
    xlabel('iteration (i)'); ylabel('\Phi({\bf{x}}^{i};{\bf{s}})-\Phi^*'); title('Fixed versus Adaptive PADMM');
    set(gca,'Box','off','XLim',[0,max(niters_uses)]);
end

