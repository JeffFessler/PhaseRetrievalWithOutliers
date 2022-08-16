function xs = solve_cubic(a,b,c,d,varargin)
% xs = solve_cubic(a,b,c,d)
% xs = solve_cubic(a,b,c,d,'real')
% xs = solve_cubic(a,b,c,d,'nonnegative')
%
% Inputs:
%  a,b,c,d - coefficients of general cubic ax^3+bx^2+cx+d = 0
%
% Outputs:
%  xs - solutions of general cubic; if fewer than 3, others are NaN
%
% This function computes the solution to the general cubic equation
% ax^3+bx^2+cx+d = 0.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% which roots?
real_only = false; nn_only = false;
if nargin >= 5
    if isequal(varargin{1},'real'), real_only = true; end
    if isequal(varargin{1},'nonnegative'), real_only = true; nn_only = true; end
end

if ~isreal(a) || ~isreal(b) || ~isreal(c) || ~isreal(d)
    warning('ignoring imaginary part(s) of a, b, c, and/or d');
    a = real(a); b = real(b); c = real(c); d = real(d);
end
[n,ii] = max([numel(a),numel(b),numel(c),numel(d)]);
switch ii
    case 1
        sz = size(a);
    case 2
        sz = size(b);
    case 3
        sz = size(c);
    case 4
        sz = size(d);
end
if isscalar(a), a = a.*ones([n,1]); else a = a(:); end
if isscalar(b), b = b.*ones([n,1]); else b = b(:); end
if isscalar(c), c = c.*ones([n,1]); else c = c(:); end
if isscalar(d), d = d.*ones([n,1]); else d = d(:); end

xs = NaN([n,3],class(a));

% take care of the nontrivial cubics first
a_nonzero = (a ~= 0);

if any(a_nonzero(:))
    p = (3.*a(a_nonzero).*c(a_nonzero) - b(a_nonzero).^2)./(3.*a(a_nonzero).^2);
    q = (2.*b(a_nonzero).^3 - 9.*a(a_nonzero).*b(a_nonzero).*c(a_nonzero) + 27.*a(a_nonzero).^2.*d(a_nonzero))./(27.*a(a_nonzero).^3);

    if real_only
        ts = solve_depressed_cubic(p,q,'real');
    else
        ts = solve_depressed_cubic(p,q);
    end

    % convert back to roots x
    xs(a_nonzero,:) = bsxfun(@minus,ts,b(a_nonzero)./(3.*a(a_nonzero)));
end

% trivial cubics (quadratics)
b_nonzero = ~a_nonzero & (b ~= 0);
if any(b_nonzero(:))
    disc = c.^2 - 4.*b.*d;
    disc_zero = b_nonzero & (disc == 0); % one duplicated real root
    disc_pos = b_nonzero & (disc > 0); % two real roots
    disc_neg = b_nonzero & (disc < 0); % two complex roots
    xs(disc_zero,1) = -c(disc_zero)./(2.*b(disc_zero));
    disc_root = sqrt(disc(disc_pos));
    xs(disc_pos,1) = (-c(disc_pos)+disc_root)./(2.*b(disc_pos));
    xs(disc_pos,2) = (-c(disc_pos)-disc_root)./(2.*b(disc_pos));
    if ~real_only
        disc_root = sqrt(disc(disc_neg));
        xs(disc_neg,1) = (-c(disc_neg)+disc_root)./(2.*b(disc_neg));
        xs(disc_neg,2) = (-c(disc_neg)-disc_root)./(2.*b(disc_neg)); 
    end
end

% trivial cubics (linear)
c_nonzero = ~a_nonzero & ~b_nonzero & (c ~= 0);
if any(c_nonzero(:))
    xs(c_nonzero,1) = -d(c_nonzero)./c(c_nonzero);
end

% take care of negative roots
if nn_only
    xs(xs < 0) = NaN;
end

if sz(end) == 1, sz(end) = 3; else sz(end+1) = 3; end
xs = reshape(xs,sz);

end
