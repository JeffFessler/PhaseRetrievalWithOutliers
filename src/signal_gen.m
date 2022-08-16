function x = signal_gen(N,K,sparse_min,sparse_max)
% x = signal_gen(N,K,sparse_min,sparse_max)
%
% Inputs:
%  N - signal length
%  K - sparsity (# of nonzeros)
%  sparse_min,sparse_max - minimum and maximum magnitude for nonzeros
%
% Outputs:
%  x - complex sparse signal
%
% This function generates a K-sparse length-N complex signal x. The
% non-zeros have magnitudes uniformly distributed between sparse_min and
% sparse_max, and phase uniformly distributed between -pi and pi.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

x = zeros(N,1);
sparsity = randperm(N,K); % length(sparsity) == K
x(sparsity) = (sparse_min + (sparse_max-sparse_min).*rand(size(sparsity))).*exp(complex(0,(2*pi).*rand(size(sparsity))));

end
