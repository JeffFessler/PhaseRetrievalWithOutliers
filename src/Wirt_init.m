function [s,x] = Wirt_init(y,A)
% [s,x] = Wirt_init(y,A)
%
% Inputs:
%  y - measurement vector
%  A - measurement matrix
%
% Outputs:
%  s - magnitudes sqrt(max(y,0)), with phase from (Ax)
%  x - unnormalized vector (optional)
%
% This function initializes x (and majorizer s) for phase retrieval based
% on the measurements y using Wirtinger flow method described in E. J.
% Candes, X. Li, and M. Soltanolkotabi, "Phase retrieval via Wirtinger
% flow: Theory and algorithms," IEEE Trans. Info. Theory, vol. 61, no. 4,
% pp. 1985-2007, Apr. 2015. The approach computes the largest eigenvector
% of A'*diag{y}*A using power iteration, but does not perform any
% normalization, since we're interested in the phase of A*x, which
% shouldn't change without a normalization factor.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

y = max(y,0);
[m,n] = size(A);

niters = 50;

% do power iteration
v0 = complex(randn(n,1),randn(n,1));
while norm(v0) < eps, v0 = complex(randn(n,1),randn(n,1)); end
v = v0./norm(v0);

for ii = 1:niters
    v = col(A'*(y.*reshape(A*v,size(y))));
    while norm(v) < eps, v = complex(randn(n,1),randn(n,1)); end
    v = v./norm(v);
end

x = v; % unnormalized
s = sqrt(y).*exp(complex(0,angle(reshape(A*v,size(y)))));

end
