function ts = solve_depressed_quartic(p,q,r,varargin)
% ts = solve_depressed_quartic(p,q,r)
% ts = solve_depressed_quartic(p,q,r,'real')
% ts = solve_depressed_quartic(p,q,r,'nonnegative')
%
% Inputs:
%  p,q,r - coefficients of depressed quartic t^4+pt^2+qt+r = 0
%
% Outputs:
%  ts - solutions of depressed quartic; if fewer than 4, others are NaN
%
% This function computes the solution to the depressed quartic equation
% t^4+pt^2+qt+r = 0.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% which roots?
real_only = false; nn_only = false;
if nargin >= 4
    if isequal(varargin{1},'real'), real_only = true; end
    if isequal(varargin{1},'nonnegative'), real_only = true; nn_only = true; end
end

if ~isreal(p) || ~isreal(q) || ~isreal(r)
    warning('ignoring imaginary part(s) of p, q and/or r');
    p = real(p); q = real(q); r = real(r);
end
[n,ii] = max([numel(p),numel(q),numel(r)]);
switch ii
    case 1
        sz = size(p);
    case 2
        sz = size(q);
    case 3
        sz = size(r);
end
if isscalar(p), p = p.*ones([n,1]); else p = p(:); end
if isscalar(q), q = q.*ones([n,1]); else q = q(:); end
if isscalar(r), r = r.*ones([n,1]); else r = r(:); end

% initialize roots to NaN's
ts = NaN([n,4],class(p));

exp_pi4 = exp(complex(0,pi/4));

% p == q == 0 => t^4 = -r
pq_zero = (p == 0) & (q == 0);
if any(pq_zero(:))
    r_zero = pq_zero & (r == 0); % one root t = 0
    r_pos = pq_zero & (r > 0); % four complex roots
    r_neg = pq_zero & (r < 0); % two real roots, two complex
    ts(r_zero,1) = 0;
    ts(r_neg,1) = (-r(r_neg)).^(1/4);
    ts(r_neg,2) = -ts(r_neg,1);
    if ~real_only
        ts(r_neg,3) = complex(0,ts(r_neg,1)); 
        ts(r_neg,4) = conj(ts(r_neg,3));
        ts(r_pos,1) = r(r_pos).^(1/4).*exp_pi4;
        ts(r_pos,2) = 1j.*ts(r_pos,1);
        ts(r_pos,3) = -ts(r_pos,1);
        ts(r_pos,4) = 1j.*ts(r_pos,3);
    end
end

% p ~= 0, q == 0 => biquadratic t^4+pt^2+r = 0
q_zero = (p ~= 0) & (q == 0);
if any(q_zero(:))
    disc = p.^2-4.*r;
    disc_zero = q_zero & (disc == 0); % real double t^2 = -p/2
    disc_pos = q_zero & (disc > 0); % two real t^2
    disc_neg = q_zero & (disc < 0); % two complex t^2

    % disc == 0 case
    z12 = -p(disc_zero)./2;
    sqrt_z12 = sqrt(z12); if real_only, sqrt_z12(z12 < 0) = NaN; end
    ts(disc_zero,1) = sqrt_z12;
    ts(disc_zero,2) = -sqrt_z12;
    
    % disc > 0 case
    disc_root = sqrt(disc(disc_pos));
    z1 = (-p(disc_pos)+disc_root)./2;
    z2 = (-p(disc_pos)-disc_root)./2;
    sqrt_z1 = sqrt(z1); if real_only, sqrt_z1(z1 < 0) = NaN; end
    sqrt_z2 = sqrt(z2); if real_only, sqrt_z2(z2 < 0) = NaN; end
    ts(disc_pos,1) = sqrt_z1;
    ts(disc_pos,2) = -sqrt_z1;
    ts(disc_pos,3) = sqrt_z2;
    ts(disc_pos,4) = -sqrt_z2;
    
    % disc < 0 case
    if ~real_only
        disc_root = sqrt(disc(disc_neg));
        z1 = (-p(disc_neg)+disc_root)./2;
        z2 = (-p(disc_neg)-disc_root)./2;
        ts(disc_neg,1) = sqrt(z1);
        ts(disc_neg,2) = -ts(disc_neg,1);
        ts(disc_neg,3) = sqrt(z2);
        ts(disc_neg,4) = -ts(disc_neg,3);
    end
end

% p ~= 0, q ~= 0 => Ferrari's method
pq_nonzero = (p ~= 0) & (q ~= 0);
if any(pq_nonzero(:))
    b = 2.5.*p(pq_nonzero); 
    c = 2.*p(pq_nonzero).^2-r(pq_nonzero);
    d = p(pq_nonzero).^3./2-p(pq_nonzero).*r(pq_nonzero)./2-q(pq_nonzero).^2./8;

    % use depressed cubic formula
    y = solve_depressed_cubic(c-b.^2/3,(2.*b.^3-9.*b.*c+27.*d)./27,'real'); % all have at least one real root
    y = y(:,1) - b./3; % convert back

    b_sqd = p(pq_nonzero) + 2.*y; % should not be equal to zero
    b1 = sqrt(b_sqd); if real_only, b1(b_sqd < 0) = NaN; end
    b2 = -b1;
    disc1 = -(3.*p(pq_nonzero)+2.*y+2.*q(pq_nonzero)./b1);
    sqrt_disc1 = sqrt(disc1); if real_only, sqrt_disc1(b_sqd < 0 | disc1 < 0) = NaN; end
    disc2 = -(3.*p(pq_nonzero)+2.*y-2.*q(pq_nonzero)./b1);
    sqrt_disc2 = sqrt(disc2); if real_only, sqrt_disc2(b_sqd < 0 | disc2 < 0) = NaN; end
    ts(pq_nonzero,1) = (b1+sqrt_disc1)./2;
    ts(pq_nonzero,2) = (b1-sqrt_disc1)./2;
    ts(pq_nonzero,3) = (b2+sqrt_disc2)./2;
    ts(pq_nonzero,4) = (b2-sqrt_disc2)./2;
end

% take care of negative real roots
if nn_only
    ts(ts < 0) = NaN;
end

if sz(end) == 1, sz(end) = 4; else sz(end+1) = 4; end
ts = reshape(ts,sz);

end
