function [x_best,f_best,info] = phase_retrieval_p1_q2(sz_x,A,y,beta,varargin)
% [x_best,f_best,info] = phase_retrieval_p1_q2(sz_x,A,y,beta,key1,value1,key2,value2,...)
%
% Inputs:
%  sz_x - size of signal/image
%  A - measurement matrix
%  y - measurements
%  beta - regularization parameter
%  keys, values - see "opts" structure below
%
% Outputs:
%  x_best - best reconstruction (over all initializations)
%  f_best - objective value for x_best
%  info - structure containing recon information (optional)
%
% This function reconstructs a compressible signal x, from measurements y,
% using the measurement matrix A and a l1 regularizer weighted by beta. The
% bulk of the proposed algorithm is in ox_fessler_p1_q2.m and in
% PADMM_1split_p1_q2.m, so a detailed description of the algorithm can be 
% found in those places. This code implements the 1-norm data fit term for
% squared-magnitude (q = 2) measurements.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

opts.Ninit_max = Inf; % number of initial guesses
opts.Ninit_min = 1; % minimum number of initial guesses to test
opts.I_mults_max = Inf; % maximum number of multiplies by A
opts.I_mults_min = 1; % minimum number of multiplies by A
opts.ox_opts = {}; % options for optimization transfer (see ox_fessler_p1_q2.m)
opts.xs_init = {}; % leave empty to use default
opts.ss_init = {}; % leave empty to randomize
opts.norm_factor = 1; % normalization for data term

opts = vararg_pair(opts,varargin);

% storage
x_best = [];
f_best = Inf;

if nargout >= 3
    info.fs = [];
    info.ox_info = [];
    tm = tic();
end

% loop
ii = 0; nAmults = 0;
while (ii < opts.Ninit_min && nAmults < opts.I_mults_max) || (nAmults < opts.I_mults_min && ii < opts.Ninit_max)
    ii = ii + 1;
    
    if ~isempty(opts.ss_init) && numel(opts.ss_init) >= ii
        s = opts.ss_init{ii};
    else
        % choose random initial majorizer point
        s_phase = (2*pi).*rand(size(y),class(y));
        s = sqrt(max(y,0)).*exp(complex(0,s_phase));
    end
    
    if ~isempty(opts.xs_init) && numel(opts.xs_init) >= ii
        x = opts.xs_init{ii};
        Ax = []; % don't precompute
    else
        x = zeros(sz_x,class(y));
        Ax = zeros(size(y),class(y)); % A*0
    end

    nAmults_prev = nAmults;
    % run optimization transfer
    if exist('info','var')
        [x,f,nAmults,ox_info] = ox_fessler_p1_q2(x,s,A,y,beta,opts.ox_opts{:},'Ax0',Ax,'norm_factor',opts.norm_factor);
        info.fs(ii) = f;
    else
        [x,f,nAmults] = ox_fessler_p1_q2(x,s,A,y,beta,opts.ox_opts{:},'Ax0',Ax,'norm_factor',opts.norm_factor);
    end
    nAmults = nAmults_prev + nAmults;
    
    if f < f_best || isempty(x_best)
        x_best = x;
        f_best = f;
        if exist('info','var')
            info.ox_info = ox_info;
        end
    end
end

if exist('info','var'), info.time = toc(tm); info.ninit = ii; info.nAmults = nAmults; end

end
