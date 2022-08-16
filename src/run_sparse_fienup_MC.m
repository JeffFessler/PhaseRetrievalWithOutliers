% run_sparse_fienup_MC.m -- Run the L1-modified sparse Fienup
% algorithm for a range of sparsities and measurements. For each generated
% signal, samples are generated with the desired noise/outliers, and the
% phase_retrieval_sparse_fienup.m function is called to run the L1-modified
% sparse Fienup algorithm for multiple initializations. The best signal
% reconstructed by the algorithm is then passed to
% find_best_match_complex.m to match it with the true signal and compute
% its RMS error. The Monte Carlo experiments average the squared errors
% over all the trials to get PSER (called PSNR here) values for the figures
% in the paper. See "TODO"s for places that require attention before
% running this code.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

setup_IRT; % code uses Image Reconstruction Toolbox

algotitle = 'L1-Fienup'; % name of algorithm

rng('default');

Ntrials = 50; % # of trials
Ninit = 50; % constrain by # of multiplies
N = 128; % signal length
Ms = ceil((2.^linspace(-3,0,7))*N); % measurements
M_sampling = 'random_plus_zero';
Ks = 3:8; % # of sparse points
sparse_min = 0; sparse_max = 1; % magnitude range for nonzeros in signal
I_mults = 10000; % default; will match # from running proposed method below
thresh_sf = 1e-4; % stopping criterion
I_ap = 50; % # of alternating projection iterations per initialization

% where to store temporary and final results
% TODO: folder must exist
outfolder = fullfile('..','..','results','phase_retrieval');

% TODO: specify parameters for simulation (examples for figures in paper
% follow)

% AWGN+5 outliers (Fig. 4)
AWGN_SNR = 40; % dB
outliers = 5; outliers_range = [1,2];
AWLN_SNR = Inf; % dB
betas_k = 0.6*Ks; % based on experiment varying betas for l1proj = true, for sparse_min = 0
outfile = 'MC_AWGN5out_sfl1'; % name of final results file
multsfiles = {'MC_AWGN5out_p1q2.mat'};

% AWLN+5 outliers (Fig. 5)
AWGN_SNR = Inf; % dB
outliers = 5; outliers_range = [1,2];
AWLN_SNR = 40; % dB
betas_k = 0.6*Ks; % based on experiment varying betas for l1proj = true, for sparse_min = 0
outfile = 'MC_AWLN5out_sfl1'; % name of final results file
multsfiles = {'MC_AWLN5out_p1q2.mat'};

%% do experiments

[Ks_grid,Ms_grid] = ndgrid(Ks,Ms); % set of all (K,M) pairs to test
betas_grid = betas_k(:)*ones(size(Ms));

% storage
errors = zeros([Ntrials,size(Ks_grid)]); % error of global minimizer for each trial
objectives_true = zeros([Ntrials,size(Ks_grid)]); % obj. function value for true signal
objectives = zeros([Ntrials,size(Ks_grid)]); % obj. function values for all initializers
xs_true = zeros([N,Ntrials,size(Ks_grid)]); % true signals for each trial
samps = false([N,Ntrials,size(Ks_grid)]); % measurement indices for each trial (for M < N)
ys = cell(size(Ks_grid)); % measurements for each trial
xs_best = zeros([N,Ntrials,size(Ks_grid)]); % global minimizer signals for each trial
infos = cell([Ntrials,size(Ks_grid)]); % recon info for each trial

%% set # of multiplies to match proposed method
I_mults_grid = I_mults.*ones(size(Ks_grid));
if exist('multsfiles','var') && ~isempty(multsfiles)
    for ifile = 1:length(multsfiles)
        try
            Smultsfile = load(fullfile(outfolder,multsfiles{ifile}),'Ks_grid','Ms_grid','infos');
            [Smultsfile.KMs_unique,~,Smultsfile.KMs_inds] = unique([Smultsfile.Ks_grid(:),Smultsfile.Ms_grid(:)],'rows','stable');
            [Smultsfile.my_KMs_unique,~,Smultsfile.my_KMs_inds] = unique([Ks_grid(:),Ms_grid(:)],'rows','stable');
            [~,ia,ib] = intersect(Smultsfile.KMs_unique,Smultsfile.my_KMs_unique,'rows');
            for ii = 1:length(ia)
                I_mults_grid(Smultsfile.my_KMs_inds == ib(ii)) = max(max(col(I_mults_grid(Smultsfile.my_KMs_inds == ib(ii)))),max(col(arrayfun(@(info) info.nAmults, Smultsfile.infos(:,Smultsfile.KMs_inds == ia(ii))))));
            end
        catch 
            fprintf(2,'Warning: unable to get # of multiplies from file \"%s\".\n',multsfiles{ifile});
        end
        clear Smultsfile ia ib ii;
    end
end

%% iterate experiments
parpool;

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
        beta = betas_grid(i12);        

        pr_ap_opts = {'I_ap',I_ap,'thresh_sf',thresh_sf};
        pr_opts = {'I_mults_min',I_mults_grid(i12),'Ninit_max',Inf,'Ninit_min',Ninit*(1 + (M >= N)),'I_mults_max',Inf,'A_unitary',true,'ap_opts',pr_ap_opts};

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
            y = sample_gen(A,x,'AWGN_SNR',AWGN_SNR,'outliers',outliers,'outliers_range',outliers_range,'AWLN_SNR',AWLN_SNR,'q',2);
            ys1(:,itrial) = y;
            norm_factor = compute_beta_normalization(y,outliers,outliers_range,AWGN_SNR,AWLN_SNR,2); % not really needed by sparse Fienup
            objectives_true1(itrial) = sum(col(abs(y-abs(A*x).^2).^2))./norm_factor;

            %% perform reconstruction
            pr_opts_use = [pr_opts,{'norm_factor',norm_factor}];
            beta = sum(abs(x(:)));
            [x_sf,objectives1(itrial),infos1{itrial}] = phase_retrieval_sparse_fienup(A,y,beta,pr_opts_use{:});

            %% evaluate
            [x_sf,errors1(itrial)] = find_best_match_complex(x_sf,x);
            xs_best1(:,itrial) = x_sf;
        end
        
        %% write temporary results
        try
            save(fullfile(outfolder,sprintf('%s_part%03d.mat',outfile,i12)),'-v7.3','errors1','objectives_true1','objectives1','infos1','xs_true1','samps1','ys1','xs_best1','K','M','beta','pr_opts');
        catch esave
            fprintf(2,'Error saving ''%s_part%03d.mat'': %s\n',outfile,i12,esave.message);
        end
    end

    %% update storage
    errors(:,i12) = errors1;
    objectives_true(:,i12) = objectives_true1;
    objectives(:,i12) = objectives1;
    infos(:,i12) = infos1;
    xs_true(:,:,i12) = xs_true1;
    samps(:,:,i12) = samps1;
    ys{i12} = ys1;
    xs_best(:,:,i12) = xs_best1;

    %% report findings
    fprintf(1,'[K = %d, M = %d]: Error: %g\n',K,M,sqrt(mean(errors1.^2)));

    clear errors1 objectives_true1 objectives1 infos1 xs_true1 samps1 ys1 xs_best1;
end

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
