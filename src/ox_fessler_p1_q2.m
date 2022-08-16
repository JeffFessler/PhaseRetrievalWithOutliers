function [x,f,nAmults,info] = ox_fessler_p1_q2(x,s,A,y,beta,varargin)
% [x,f,nAmults,info] = ox_fessler_p1_q2(x,s,A,y,beta,key1,value1,key2,value2,...)
%
% Inputs:
%  x - initial guess (generally zero)
%  s - initial majorization vector
%  A - measurement matrix
%  y - measurements
%  beta - regularization parameter
%  keys, values - see "opts" structure below
%
% Outputs:
%  x - best reconstruction (over all iterations)
%  f - objective value for x
%  nAmults - number of multiplications
%  info - structure containing recon information (optional)
%
% This function reconstructs a compressible signal x, from measurements y,
% using the measurement matrix A and a l1 regularizer weighted by beta.
% This function contains the majorize-minimize portion of the proposed
% method. Given an initial guess and majorization vector, it calls
% PADMM_1split_p1_q2.m to find a new x that minimizes the majorizer, and
% then it updates the majorization vector based on that new x and repeats
% the process until convergence or the maximum # of iterations is reached.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

tm = tic();

opts.eps_ox = 0.01; % convergence in true objective function threshold
opts.I_ox = 10; % max # of iterations
opts.ADMM_opts = {}; % options for ADMM/PADMM (see PADMM_1split_p1_q2.m)
opts.Ax0 = []; % (optional) precomputed A*x
opts.norm_factor = 1;

opts = vararg_pair(opts,varargin);

Psifun = @(x,Ax) norm(col(y-abs(Ax).^2),1)./opts.norm_factor + beta*sum(abs(col(x))); % true objective

% initialize
nAmults = 0;
if ~isempty(opts.Ax0), Ax = opts.Ax0; else Ax = A*x; nAmults = 1; end
if isempty(s), s = Ax; end
weights = 1;
f = Inf;

if nargout >= 4
    info.fs = [f;zeros(opts.I_ox,1)];
    info.Phidiffs = zeros(opts.I_ox,1);
    info.xdiffs = zeros(opts.I_ox,1);
    info.sdiffs = zeros(opts.I_ox,1);
    info.ADMM_info = cell(opts.I_ox,1);
end

% loop
iter = 0; u = []; lm = []; mus_last = [];
while iter < opts.I_ox
    % store previous
    x_prev = x; 
    s_prev = s;
    f_prev = f;
    
    iter = iter + 1; nAmults_prev = nAmults;
    
    % compute random weights for this iteration
    Phi_prev = surrogate_p1_q2(x,s,Ax,y,weights./opts.norm_factor,beta);

    opts_extra = {'Ax0',Ax,'norm_factor',opts.norm_factor};
    if ~isempty(u), opts_extra = [opts_extra,{'u0',u}]; end %#ok<AGROW>
    if ~isempty(lm), opts_extra = [opts_extra,{'lm_u0',lm}]; end %#ok<AGROW>
    if ~isempty(mus_last), opts_extra = [opts_extra,{'mu_u',mus_last}]; end %#ok<AGROW>
    
    % update x using ADMM
    if exist('info','var')
        [x,Phi,Ax,u,lm,nAmults,mus_last,info.ADMM_info{iter}] = PADMM_1split_p1_q2(x,s,A,y,weights,beta,opts.ADMM_opts{:},opts_extra{:});
    else
        [x,Phi,Ax,u,lm,nAmults,mus_last] = PADMM_1split_p1_q2(x,s,A,y,weights,beta,opts.ADMM_opts{:},opts_extra{:});
    end
    if ~isscalar(Phi), [Phi,ibest] = min(Phi(:)); x = x{ibest}; Ax = Ax{ibest}; u = u{ibest}; lm = lm{ibest}; nAmults = nAmults(ibest); end % case for multiple mu's
    if exist('info','var')
        info.Phidiffs(iter) = Phi_prev - Phi; 
        info.xdiffs(iter) = norm(col(x-x_prev)); 
    end
    nAmults = nAmults_prev + nAmults;
    
    % update s
    s = Ax;
    if exist('info','var'), info.sdiffs(iter) = norm(col(s-s_prev)); end
    
    % update and check objective
    f = Psifun(x,Ax);
    if exist('info','var'), info.fs(iter+1) = f; end
    if abs(f_prev - f) < opts.eps_ox, break; end
end

if exist('info','var')
    info.time = toc(tm); 
    info.niters = iter; 
    info.ADMM_info = cat(1,info.ADMM_info{1:iter});
end

end

