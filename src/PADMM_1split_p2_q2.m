function [x_best,Phi_best,Ax_best,u_best,lm_best,nAmults,mus_last,info] = PADMM_1split_p2_q2(x0,s,A,y,weights,beta,varargin)
% [x_best,Phi_best,Ax_best,u_best,lm_best,nAmults,mus_last,info] = PADMM_1split_p2_q2(x0,s,A,y,weights,beta,key1,value1,key2,value2,...)
%
% Inputs:
%  x0 - initial guess
%  s - majorization vector
%  A - measurement matrix
%  y - measurements
%  weights - weights for data fit term
%  beta - regularization parameter
%  keys, values - see "opts" structure below
%
% Outputs:
%  x_best - best reconstruction (over all iterations)
%  Phi_best - majorizer value for x_best
%  Ax_best - A*x_best
%  u_best - u for best x
%  lm_best - lagrange multiplier for best x
%  nAmults - number of multiplications
%  mus_last - last value of mu for adaptive method
%  info - structure containing recon information (optional)
%
% This function performs the ADMM or preconditioned variant for the
% proposed method with the quadratic data fit term. It iterates the
% x-update (using either shrinkage if A is unitary or FISTA-like shrinkage
% if A is not), the u-update (using the calculations in the supplement),
% and the lagrange multipliers update. The algorithm also includes
% (optional) automatic tuning of mu. 
% 
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

tm = tic();

% default option values
opts.mu_u = 1; % initial value of mu (NEW: can use vector to repeat for various mu's)
opts.mu_I_adapt = 0; % set to > 0 to adapt every mu_I_adapt iterations
opts.mu_thresh = 2; % factor of primal residual being larger or smaller than dual residual
opts.mu_scale = 2; % factor to multiply or divide mu during adaptive step
opts.eps_ADMM = 0.01; % threshold for convergence (set to zero to disable)
opts.I_ADMM_min = 1; % minimum # of iterations for ADMM before termination
opts.I_ADMM = 10; % maximum # of iterations for ADMM before early termination
opts.A_unitary = false; % if true, use shrinkage directly; if false, use FISTA for x-update
opts.A_curvature = []; % curvature for FISTA; will compute if needed
opts.Ax0 = []; % (optional) precomputed A*x0
opts.u0 = []; % can resume from previous u
opts.lm_u0 = []; % can resume from previous d_u = lm_u/mu
opts.verbose = false; % if true, include iteration details in info structure
opts.norm_factor = 1;

opts = vararg_pair(opts,varargin);

% determine what sub-algorithm to use for x-update
X_ALG_SHRINK = 0;
X_ALG_PADMM = 2;
if opts.A_unitary % use shrinkage for x-update
    x_alg = X_ALG_SHRINK;
else % use FISTA-like x-update
    x_alg = X_ALG_PADMM;
end

% initialize
% sz_x = size(x0); 
nAmults = 0;
if isempty(opts.Ax0), opts.Ax0 = A*x0; nAmults = 1; end
if isempty(opts.u0), opts.u0 = opts.Ax0; end
if isempty(opts.lm_u0), opts.lm_u0 = zeros(size(opts.u0),class(opts.u0)); end
if isempty(s), s = opts.u0; end
if x_alg == X_ALG_PADMM, t0 = 1; end

% for tighter majorizer when |s|^q > y > 0
use_y_s = (y > 0 & abs(s).^2 > y);
s(use_y_s) = sqrt(y(use_y_s)).*exp(complex(0,angle(s(use_y_s))));

Phifun = @(x,Ax) surrogate_p2_q2(x,s,Ax,y,weights./opts.norm_factor,beta);

Phi0 = Phifun(x0,opts.Ax0); % initial function value; Phi(x_best) must be as good as this

    function out = power_iterate(x)
        out = col(A'*(A*x));
        nAmults = nAmults + 2*size(x,2);
    end

if isempty(opts.A_curvature) && x_alg == X_ALG_PADMM
    % calculate curvature if needed & not provided
    if isnumeric(A)
        opts.A_curvature = normest(A).^2;
    else
        % power iteration
        opts.A_curvature = eigs(@power_iterate,numel(x0),1,'lm',struct('issym',1,'isreal',0));
    end
end

% iterate over mu's
nmus = numel(opts.mu_u);
Phi_best = zeros(1,nmus);
x_best = cell(1,nmus);
Ax_best = cell(1,nmus);
u_best = cell(1,nmus);
lm_best = cell(1,nmus);
nAmults = nAmults*ones(1,nmus);

if nargout >= 8
    info.alg = x_alg;
    info.niters = zeros(1,nmus);
    if opts.verbose
        info.Phis = [Phi0;zeros(opts.I_ADMM,1)]*ones(1,nmus);
        info.rps = zeros(opts.I_ADMM,nmus);
        info.xdiffs = zeros(opts.I_ADMM,nmus);
        info.udiffs = zeros(opts.I_ADMM,nmus);
        if opts.mu_I_adapt > 0
            info.mus = ones(opts.I_ADMM+1,1)*(opts.mu_u(:).'); 
            info.rds = zeros(opts.I_ADMM,nmus);
        end
        switch x_alg
            case X_ALG_PADMM
                info.ts = zeros(opts.I_ADMM,nmus);
        end
    end
    info.times = toc(tm).*ones(1,nmus);
end

for imu = 1:nmus
    tm = tic();
    x = x0; Ax = opts.Ax0;
    u = opts.u0; d_u = opts.lm_u0./opts.mu_u(imu);
    if x_alg == X_ALG_PADMM, z = x0; Az = Ax; t = t0; end
    Phi = Inf; % avoid convergence in first iteration
    Phi_best(imu) = Phi0; % track best majorizer objective function value
    x_best{imu} = x; % track best x (NEW: for each mu)
    Ax_best{imu} = Ax; % A*x for best x
    u_best{imu} = u;
    lm_best{imu} = opts.lm_u0;
    
    % loop
    iter = 0;
    while iter < opts.I_ADMM
        % store previous
        x_prev = x; Ax_prev = Ax;
        u_prev = u;
        Phi_prev = Phi;
        if x_alg == X_ALG_PADMM, t_prev = t; end
        
        iter = iter + 1;
        
        % update x
        switch x_alg
            case X_ALG_SHRINK
                % x <- argmin_x mu_u/2*||x-A'(u-d_u)||^2 + beta*||x||_1
                x = shrink(A'*(u-d_u),beta/opts.mu_u(imu));
                Ax = A*x;
                nAmults(imu) = nAmults(imu) + 2;
                
            case X_ALG_PADMM
                % x <- argmin_x mu_u*c/2*||x-(z-A'(Az-u+d_u)/c)||^2 + beta*||x||_1
                x = shrink(z-(A'*(Az-u+d_u))./opts.A_curvature,beta/opts.mu_u(imu)/opts.A_curvature);
                Ax = A*x;
                nAmults(imu) = nAmults(imu) + 2;
                if (x-z)'*(x-x_prev) < 0
                    t = 1;
                    z = x; Az = Ax;
                else
                    t = 0.5 + sqrt(0.25 + t.^2);
                    z = x + ((t_prev - 1)/t).*(x-x_prev);
                    Az = Ax + ((t_prev - 1)/t).*(Ax-Ax_prev);
                end
                if opts.verbose && exist('info','var')
                    info.ts(iter,imu) = t;
                end
        end
        
        if opts.verbose && exist('info','var'), info.xdiffs(iter,imu) = norm(col(x-x_prev)); end
        
        % check for convergence
        Phi = Phifun(x,Ax);
        rp = norm(col(Ax-u),2);
        if Phi <= Phi_best(imu) % update best solution
            Phi_best(imu) = Phi;
            x_best{imu} = x;
            Ax_best{imu} = Ax;
            u_best{imu} = u;
            lm_best{imu} = opts.mu_u(imu).*d_u;
        end
        if opts.verbose && exist('info','var'), info.Phis(iter+1,imu) = Phi; info.rps(iter,imu) = rp; end
        
        % check convergence
        if iter >= opts.I_ADMM_min && abs(Phi_prev - Phi) < opts.eps_ADMM && rp < opts.eps_ADMM, break; end
        if iter >= opts.I_ADMM_min && iter >= opts.I_ADMM && Phi_best(imu) < Phi0, break; end
        
        % update u
        etas = (opts.norm_factor*opts.mu_u(imu)/2)./weights;
        d = Ax+d_u;
        
        u = update_us_p2_q2(d,y,s,etas);
        
        if opts.verbose && exist('info','var'), info.udiffs(iter,imu) = norm(col(u-u_prev)); end
        
        % update dual variable
        d_u = d_u + (Ax-u);
        
        if opts.mu_I_adapt > 0 && mod(iter,opts.mu_I_adapt) == 0
            % adapt mu according to heuristic from boyd:10:doa
            epsp = max(norm(col(Ax),2),norm(col(u),2));
            epsd = opts.mu_u(imu)*norm(col(A'*d_u),2);
            rd = opts.mu_u(imu)*norm(col(A'*(u-u_prev)),2); % dual residual
            nAmults(imu) = nAmults(imu) + 2;
            if rd*epsp > opts.mu_thresh*rp*epsd % dual >> primal => decrease mu
                opts.mu_u(imu) = opts.mu_u(imu)/opts.mu_scale;
                d_u = d_u.*opts.mu_scale;
            elseif rp*epsd > opts.mu_thresh*rd*epsp % dual << primal => increase mu
                opts.mu_u(imu) = opts.mu_u(imu)*opts.mu_scale;
                d_u = d_u./opts.mu_scale;
            end
            if opts.verbose && exist('info','var')
                info.rds(iter,imu) = rd;
                info.mus(iter+1:end,imu) = opts.mu_u(imu); 
            end
        end
    end
    
    if exist('info','var'), info.niters(imu) = iter; info.times(imu) = info.times(imu) + toc(tm); end

end

if nmus == 1, x_best = x_best{1}; Ax_best = Ax_best{1}; u_best = u_best{1}; lm_best = lm_best{1}; end
mus_last = opts.mu_u;

end
