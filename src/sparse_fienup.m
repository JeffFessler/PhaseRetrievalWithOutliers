function [x_best,f_best,nAmults,info] = sparse_fienup(x,A,y,beta,nAmults,varargin)
% [x_best,f_best,nAmults,info] = sparse_fienup(x,A,y,beta,nAmults,key1,value1,key2,value2,...)
%
% Inputs:
%  x - initial guess
%  A - measurement matrix
%  y - measurements
%  beta - l1-ball radius constraint
%  nAmults - # of multiplies already
%  keys, values - see "opts" structure below
%
% Outputs:
%  x_best - best reconstruction (during alternating projection iterations)
%  f_best - data discrepancy for x_best
%  nAmults - updated # of multiplies
%  info - structure containing recon information (optional)
%
% This function reconstructs a compressible signal x, from measurements y,
% using an initial guess x, the measurement matrix A, and a l1-sparsity
% constraint ||x||_1 <= beta. This algorithm is essentially the sparse
% Fienup algorithm described in
%
% S. Mukherjee and C. S. Seelamantula, "An iterative algorithm for phase
% retrieval with sparsity constraints: application to frequency domain
% optical coherence tomography," in Proc. IEEE Conf. Acoust. Speech Sig.
% Proc., 2012, pp. 553-6.
%
% The algorithm alternates between enforcing the l1-norm constraint and the
% measurements. The l1-norm constraint is enforced using a simple sorting
% algorithm (see l1proj.m), while the measurement constraint is enforced by
% transforming the signal, keeping the phase while setting the magnitudes
% to match the (square root) measurements, and projecting the signal in the
% signal domain to the subspace corresponding to this transform domain
% constraint. In the DFT case, this is the same as transforming the DFT
% measurements back to the image domain.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

opts.I_ap = 10; % max # of iterations
opts.thresh_sf = 1e-4; % stopping criterion for xdiff
opts.A_unitary = false; % set to true to use A' instead of A'(AA')^{-1}
opts.pcg_opts{:} = {};
opts.I_mults = Inf; % maximum # of multiplications by A
opts.verbose = false; % info contains verbose details for all iterations
opts.norm_factor = 1; % not really used

opts = vararg_pair(opts,varargin);

% determine what sub-algorithm to use for x-update
X_ALG_SHRINK = 0;
x_alg = X_ALG_SHRINK;

% initialize
% sz_x = size(x);
Ax = A*x; nAmults = nAmults + 1;

if nargout >= 4
    info = struct('alg',x_alg);
    if opts.verbose
        info.xdiffs = zeros(opts.I_ap,1);
        info.fs = zeros(opts.I_ap+1,1);
    end
    info.nAmults = 0;
    tm = tic();
end

% loop
iter = 0;
f = norm(col(abs(Ax).^2-y))^2/opts.norm_factor;
if sum(abs(x(:)),1) <= beta
    f_best = f; % initial function value; f(x_best) must be as good as this
else
    f_best = Inf;
end
x_prev = x;
x_best = x;
    
if opts.verbose && exist('info','var'), info.fs(1) = f; end
while iter < opts.I_ap && nAmults < opts.I_mults
    % next iteration
    iter = iter + 1;
    
    % update x using sparsity penalty
    % update x
    switch x_alg
        case X_ALG_SHRINK
            % x <- argmin_x 1/2*||x-xhat||^2 s.t. ||x||_1 <= beta
            x = l1proj(x,beta);
    end
    Ax = A*x;
    nAmults = nAmults + 1;
    
    f = norm(col(abs(Ax).^2-y))^2/opts.norm_factor;
    xdiff = norm(col(x-x_prev));
    if opts.verbose && exist('info','var'), info.fs(iter+1) = f; info.xdiffs(iter) = xdiff; end
    
    if f <= f_best
        f_best = f;
        x_best = x;
    end
    
    if xdiff <= opts.thresh_sf*norm(col(x_prev))
        break;
    end
    
    % store previous
    x_prev = x;
    
    % update x using projection to data
    x = sqrt(max(y,0)).*exp(complex(0,angle(Ax))) - Ax;
    if ~opts.A_unitary
        [x,pcg_info] = my_pcg(x,A*A',x,opts.pcg_opts{:});
        x = reshape(x,size(y));
        nAmults = nAmults + 2*(1+size(pcg_info,1));
    end
    x = x_prev + A'*x;
    nAmults = nAmults + 1;
    
end

if exist('info','var'), info.time = toc(tm); info.niters = iter; info.nAmults = nAmults; end

end
