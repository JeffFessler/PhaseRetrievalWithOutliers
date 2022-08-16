function [y,varargout] = sample_gen(A,x,varargin)
% [y,y_true] = sample_gen(A,x,key1,value1,key2,value2,...)
%
% Inputs:
%  A - measurement matrix
%  x - true signal
%  keys, values - see "opts" structure below
%
% Outputs:
%  y - squared-magnitude measurements
%  y_true - true value of measurements (no noise or outliers) (optional)
%
% This function generates samples from the signal x, using the system
% matrix A, and adding (optional) noise and/or outliers to the signal.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

opts.AWGN_SNR = Inf; % dB scale; if finite, add AWGN to samples
opts.outliers_oldstyle = false; % if true, use ICIP version (set to 2*max amplitude)
opts.outliers = 0; % if > 0, add outliers to samples (uniform around true value)
opts.outliers_range = [1,1]; % if outliers > 0, sets range of outliers; relative to max(|y|)^q
opts.AWLN_SNR = Inf; % dB scale; if finite, add additive white Laplace noise to samples
opts.q = 2; % ytrue_i = |[A*x]_i|^q
opts.nonnegative = true; % threshold samples to be nonnegative

opts = vararg_pair(opts,varargin);

y = abs(A*x).^opts.q;
max_y = max(y(:));
norm_y = norm(y(:));

if nargout >= 2 % keep noise-free measurements
    varargout{1} = y;
end

% add outliers
if opts.outliers_oldstyle && opts.outliers > 0
    outliers_ind = randperm(numel(y),min(numel(y),opts.outliers));
    y(outliers_ind(:)) = 2*max_y;
end
if ~opts.outliers_oldstyle && opts.outliers > 0 % outliers at random locations, uniformly distributed
    outliers_ind = randperm(numel(y),min(numel(y),opts.outliers));
    outliers = (opts.outliers_range(2)-opts.outliers_range(1)).*rand(size(outliers_ind),class(y)) + opts.outliers_range(1);
    y(outliers_ind(:)) = outliers(:).*max_y;
end

% add Gaussian noise
if isfinite(opts.AWGN_SNR)
    AWGN_stddev = 10^(-opts.AWGN_SNR/20)*norm_y/sqrt(numel(y));
    y = y + AWGN_stddev.*randn(size(y),class(y));
end

% add Laplace noise
if isfinite(opts.AWLN_SNR)
    AWLN_stddev = 10^(-opts.AWLN_SNR/20)*norm_y/sqrt(numel(y));
    
    AWLN = rand(size(y),class(y));
    neg = AWLN < 0.5;
    % for < 0.5, U = 1/2*exp(outlier/sqrt(var/2)) =>
    % sqrt(var/2)*log(2U) = outlier
    AWLN(neg) = (AWLN_stddev/sqrt(2)).*log(2.*AWLN(neg));
    % for > 0.5, U = 1-1/2*exp(-outlier/sqrt(var/2)) =>
    % -sqrt(var/2)*log(2(1-U)) = outlier
    AWLN(~neg) = (-AWLN_stddev/sqrt(2)).*log(2.*(1-AWLN(~neg)));
    y = y + AWLN;
end

if opts.nonnegative
    y(y < 0) = 0; % threshold
end

end
