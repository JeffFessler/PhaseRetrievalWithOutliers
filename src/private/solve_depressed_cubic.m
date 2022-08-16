function ts = solve_depressed_cubic(p,q,varargin)
% ts = solve_depressed_cubic(p,q)
% ts = solve_depressed_cubic(p,q,'real')
% ts = solve_depressed_cubic(p,q,'nonnegative')
%
% Inputs:
%  p,q - coefficients of depressed cubic t^3+pt+q = 0
%
% Outputs:
%  ts - solutions of depressed cubic; if fewer than 3, others are NaN
%
% This function computes the solution to the depressed cubic equation
% t^3+pt+q = 0.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% which roots?
real_only = false; nn_only = false;
if nargin >= 3
    if isequal(varargin{1},'real'), real_only = true; end
    if isequal(varargin{1},'nonnegative'), real_only = true; nn_only = true; end
end

if ~isreal(p) || ~isreal(q)
    warning('ignoring imaginary part(s) of p and/or q');
    p = real(p); q = real(q);
end

[n,ii] = max([numel(p),numel(q)]);
switch ii
    case 1
        sz = size(p);
    case 2
        sz = size(q);
end

% initialize roots to NaN's
ts = NaN([n,3],class(p));

disc = p.^3./27 + q.^2./4;
if isscalar(p), p = p.*ones([n,1]); else p = p(:); end
if isscalar(q), q = q.*ones([n,1]); else q = q(:); end

exp_2pi3 = exp(complex(0,2*pi/3));

% fill roots for p == 0 case (t^3 = -q)
p_zero = (p == 0);
if any(p_zero(:))
    q_nonzero = p_zero & (q ~= 0); % p == 0 & q ~= 0 => 3 roots
    ts(p_zero,1) = -sign(q(p_zero)).*abs(q(p_zero)).^(1/3);
    if ~real_only
        ts(q_nonzero,2) = ts(q_nonzero,1).*exp_2pi3;
        ts(q_nonzero,3) = conj(ts(q_nonzero,2)); 
    end
end

% fill roots for q == 0 case (t = 0 or t = sqrt(-p))
q_zero = ~p_zero & (q == 0);
if any(q_zero(:))
    p_neg = q_zero & (p < 0);
    ts(q_zero,1) = 0;
    if real_only
        ts(p_neg,2) = sqrt(-p(p_neg)); 
        ts(p_neg,3) = -ts(p_neg,2); 
    else
        ts(q_zero,2) = sqrt(-p(q_zero)); 
        ts(q_zero,3) = -ts(q_zero,2); 
    end
end

% fill (all real) roots for p ~= 0, disc == 0 case (t1 = -3q/2p, t2 = 3q/p)
disc_zero = ~p_zero & ~q_zero & (disc == 0);
if any(disc_zero(:))
    ts(disc_zero,1) = 3*q(disc_zero)./p(disc_zero);
    ts(disc_zero,2) = -0.5.*ts(disc_zero,1);
end

% fill (3 real) roots for negative discriminant case
disc_neg = ~p_zero & ~q_zero & (disc < 0);
if any(disc_neg(:))
    u = (-q(disc_neg)./2+sqrt(disc(disc_neg))).^(1/3);
    v = -p(disc_neg)./(3.*u);
    ts(disc_neg,1) = real(u + v);
    ts(disc_neg,2) = real(u.*exp_2pi3 + v.*conj(exp_2pi3));
    ts(disc_neg,3) = real(u.*conj(exp_2pi3) + v.*exp_2pi3);
end

% fill (1 real, 2 cplx) roots for positive discriminant case
disc_pos = ~p_zero & ~q_zero & (disc > 0);
if any(disc_pos(:))
    v = sqrt(disc(disc_pos));
    u = -0.5.*q(disc_pos)+v; 
    u = sign(u).*abs(u).^(1/3);
    v = -0.5.*q(disc_pos)-v; 
    v = sign(v).*abs(v).^(1/3);
    ts(disc_pos,1) = u + v;
    if ~real_only
        ts(disc_pos,2) = u.*exp_2pi3 + v.*conj(exp_2pi3);
        ts(disc_pos,3) = conj(ts(disc_pos,2));
    end
end

% take care of negative real roots
if nn_only
    ts(ts < 0) = NaN;
end

if sz(end) == 1, sz(end) = 3; else sz(end+1) = 3; end
ts = reshape(ts,sz);

end
