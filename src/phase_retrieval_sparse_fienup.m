function [x_best,f_best,info] = phase_retrieval_sparse_fienup(A,y,beta,varargin)
% [x_best,f_best,info] = phase_retrieval_sparse_fienup(A,y,beta,key1,value1,key2,value2,...)
%
% Inputs:
%  A - measurement matrix
%  y - measurements
%  beta - l1-ball radius constraint
%  keys, values - see "opts" structure below
%
% Outputs:
%  x_best - best reconstruction (over all initializations)
%  f_best - data discrepancy for x_best
%  info - structure containing recon information (optional)
%
% This function reconstructs a compressible signal x, from measurements y,
% using the measurement matrix A and a l1-sparsity constraint ||x||_1 <=
% beta. The bulk of the sparse Fienup algorithm is in sparse_fienup.m, so a
% detailed description of the algorithm can be found there.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

opts.Ninit_max = Inf; % max # of number of initial guesses
opts.Ninit_min = 1; % minimum number of initial guesses to test
opts.ap_opts = {}; % options for alternating projections (see sparse_fienup.m)
opts.A_unitary = false; % set to true to use A' instead of A'(AA')^{-1}
opts.pcg_opts = {}; % options for solving data projection step
opts.xs_init = {}; % leave empty to use ss_init
opts.ss_init = {}; % leave empty to use random
opts.I_mults_max = Inf; % max # of multiplications by A
opts.I_mults_min = 1; % min # of multiplications by A
opts.norm_factor = 1; % normalization for objective function (not really used)

opts = vararg_pair(opts,varargin);

% storage
x_best = [];
f_best = Inf;

if nargout >= 3
    info.fs = [];
    tm = tic();
end

% loop
ii = 0; nAmults = 0;
while (nAmults < opts.I_mults_min && ii < opts.Ninit_max) || (ii < opts.Ninit_min && nAmults < opts.I_mults_max)
    ii = ii + 1;
    nAmults_prev = nAmults; nAmults = 0;
    if ~isempty(opts.xs_init) && numel(opts.xs_init) >= ii
        x = opts.xs_init{ii};
    else
        if ~isempty(opts.ss_init) && numel(opts.ss_init) >= ii
            x = opts.ss_init{ii};
        else
            x_phase = (2*pi).*rand(size(y),class(y));
            x = sqrt(max(y,0)).*exp(complex(0,x_phase));
        end
        if ~opts.A_unitary
            [x,pcg_info] = my_pcg(x,A*A',x,opts.pcg_opts{:});
            x = reshape(x,size(y));
            nAmults = nAmults + 2*(1+size(pcg_info,1));
        end
        x = A'*x; % projection
        nAmults = nAmults + 1;
    end

    if exist('info','var')
        [x,f,nAmults,ap_info] = sparse_fienup(x,A,y,beta,nAmults,'A_unitary',opts.A_unitary,'I_mults',opts.I_mults_max-nAmults-nAmults_prev,'pcg_opts',opts.pcg_opts,opts.ap_opts{:},'norm_factor',opts.norm_factor);
        info.fs(ii) = f;
    else
        [x,f,nAmults] = sparse_fienup(x,A,y,beta,nAmults,'A_unitary',opts.A_unitary,'I_mults',opts.I_mults_max-nAmults-nAmults_prev,'pcg_opts',opts.pcg_opts,opts.ap_opts{:},'norm_factor',opts.norm_factor);
    end
    nAmults = nAmults_prev + nAmults;
    
    if f < f_best || isempty(x_best)
        x_best = x;
        f_best = f;
        if exist('info','var'), info.ap_info = ap_info; end
    end
end

if exist('info','var'), info.time = toc(tm); info.ninit = ii; info.nAmults = nAmults; end

end
