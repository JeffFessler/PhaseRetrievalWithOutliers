function xs = solve_quartic(a,b,c,d,e,varargin)
% xs = solve_quartic(a,b,c,d,e)
% xs = solve_quartic(a,b,c,d,e,'real')
% xs = solve_quartic(a,b,c,d,e,'nonnegative')
%
% Inputs:
%  a,b,c,d,e - coefficients of general quartic ax^4+bx^3+cx^2+dx+e = 0
%
% Outputs:
%  xs - solutions of general quartic; if fewer than 4, others are NaN
%
% This function computes the solution to the general quartic equation
% ax^4+bx^3+cx^2+dx+e = 0.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% which roots?
real_only = false; nn_only = false;
if nargin >= 6
    if isequal(varargin{1},'real'), real_only = true; end
    if isequal(varargin{1},'nonnegative'), real_only = true; nn_only = true; end
end

if ~isreal(a) || ~isreal(b) || ~isreal(c) || ~isreal(d) || ~isreal(e)
    warning('ignoring imaginary part(s) of a, b, c, d, and/or e');
    a = real(a); b = real(b); c = real(c); d = real(d); e = real(e);
end
[n,ii] = max([numel(a),numel(b),numel(c),numel(d),numel(e)]);
switch ii
    case 1
        sz = size(a);
    case 2
        sz = size(b);
    case 3
        sz = size(c);
    case 4
        sz = size(d);
    case 5
        sz = size(e);
end
if isscalar(a), a = a.*ones([n,1]); else a = a(:); end
if isscalar(b), b = b.*ones([n,1]); else b = b(:); end
if isscalar(c), c = c.*ones([n,1]); else c = c(:); end
if isscalar(d), d = d.*ones([n,1]); else d = d(:); end
if isscalar(e), e = e.*ones([n,1]); else e = e(:); end

xs = NaN([n,4],class(a));

% take care of the nontrivial quartics first
a_nonzero = (a ~= 0);

if any(a_nonzero(:))
    p = (8.*a(a_nonzero).*c(a_nonzero)-3.*b(a_nonzero).^2)./(8.*a(a_nonzero).^2);
    q = (b(a_nonzero).^3-4.*a(a_nonzero).*b(a_nonzero).*c(a_nonzero)+8.*a(a_nonzero).^2.*d(a_nonzero))./(8.*a(a_nonzero).^3);
    r = (-3.*b(a_nonzero).^4+256.*a(a_nonzero).^3.*e(a_nonzero)-64.*a(a_nonzero).^2.*b(a_nonzero).*d(a_nonzero)+16.*a(a_nonzero).*b(a_nonzero).^2.*c(a_nonzero))./(256.*a(a_nonzero).^4);

    if real_only
        ts = solve_depressed_quartic(p,q,r,'real');
    else
        ts = solve_depressed_quartic(p,q,r);
    end

    % convert back to roots x
    xs(a_nonzero,:) = bsxfun(@minus,ts,b(a_nonzero)./(4.*a(a_nonzero)));
end

% trivial quartics (cubic)
b_nonzero = ~a_nonzero & (b ~= 0);
if any(b_nonzero(:))
    p = (3.*b(b_nonzero).*d(b_nonzero)-c(b_nonzero).^2)./(3.*b(b_nonzero).^2);
    q = (2.*c(b_nonzero).^3-9.*b(b_nonzero).*c(b_nonzero).*d(b_nonzero)+27.*b(b_nonzero).^2.*e(b_nonzero))./(27.*b(b_nonzero).^3);
    if real_only
        ts = solve_depressed_cubic(p,q,'real');
    else
        ts = solve_depressed_cubic(p,q);
    end
    xs(b_nonzero,1:3) = bsxfun(@minus,ts,c(b_nonzero)./(3.*b(b_nonzero)));
end

% trivial quartics (quadratic)
c_nonzero = ~a_nonzero & ~b_nonzero & (c ~= 0);
if any(c_nonzero(:))
    disc = d.^2 - 4.*c.*e;
    disc_zero = c_nonzero & (disc == 0); % one duplicated real root
    disc_pos = c_nonzero & (disc > 0); % two real roots
    disc_neg = c_nonzero & (disc < 0); % two complex roots
    xs(disc_zero,1) = -d(disc_zero)./(2.*c(disc_zero));
    disc_root = sqrt(disc(disc_pos));
    xs(disc_pos,1) = (-d(disc_pos)+disc_root)./(2.*c(disc_pos));
    xs(disc_pos,2) = (-d(disc_pos)-disc_root)./(2.*c(disc_pos));
    if ~real_only
        disc_root = sqrt(disc(disc_neg));
        xs(disc_neg,1) = (-d(disc_neg)+disc_root)./(2.*c(disc_neg)); 
        xs(disc_neg,2) = (-d(disc_neg)-disc_root)./(2.*c(disc_neg));
    end
end

% trivial quartics (linear)
d_nonzero = ~a_nonzero & ~b_nonzero & ~c_nonzero & (d ~= 0);
if any(d_nonzero(:))
    xs(d_nonzero,1) = -e(d_nonzero)./d(d_nonzero);
end

% take care of negative roots
if nn_only
    xs(xs < 0) = NaN;
end

if sz(end) == 1, sz(end) = 4; else sz(end+1) = 4; end
xs = reshape(xs,sz);

end
