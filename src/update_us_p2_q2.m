function [u,F,indicators] = update_us_p2_q2(d,y,s,etas)
% [u,F,indicators] = update_us_p2_q2(d,y,s,etas)
%
% Inputs:
%  d - the sum A*x+d_u from PADMM_1split_p2_q2.m
%  y - measurements
%  s - majorization vector
%  etas - (norm_factor*mu/2)./weights from PADMM_1split_p2_q2.m
%
% Outputs:
%  u - solution of u-update step
%  F - majorizer function value
%  indicators - 1, 2, or 3 for u_+, u_-, or u_{+-} solution
%
% This function computes the minimizer of the majorizer function for
% quadratic data fit term with squared-magnitude measurements (see
% supplement for derivations).
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

sz_u = size(y); % u and y have same size
y = y(:); d = d(:); s = s(:); etas = etas(:);
u = zeros(size(y),class(y));
F = zeros(size(y),class(y));

phase_s = angle(s);
abs_s = abs(s);

if nargout >= 3, indicators = zeros(size(y)); end

%%%%%%%%%%%%%%%%%%
% find the minimum of f1 %
%%%%%%%%%%%%%%%%%%
phase_d = angle(d);
abs_d = abs(d);
p=0.5.*(etas-2*y); q=-0.5*etas.*abs_d;

% get nonnegative real roots of (cubic) derivative of f1
u1 = solve_depressed_cubic(p,q,'nonnegative'); % roots of r^3+pr+q=0
abs_u1 = abs(u1);

% get global minimum
F1_u1 = bsxfun(@minus,abs_u1.^2,y).^2+bsxfun(@times,etas,bsxfun(@minus,abs_u1,abs_d).^2);
F1_u1(isnan(abs_u1)) = NaN;
[F1_u1,iu1] = min(F1_u1,[],2);
abs_u1 = abs_u1(sub2ind(size(abs_u1),(1:size(abs_u1,1)).',iu1));
u1 = abs_u1.*exp(complex(0,phase_d));

% get F2(u1) for the global minimizer u1 of F1
F2_u1 = (y+abs_s.^2-2.*real(conj(s).*u1)).^2+bsxfun(@times,etas,(abs_u1-abs_d).^2);
F2_u1(isnan(abs_u1)) = NaN;

% solution iff u1 satisfies F1(u1) >= F2(u1)
u1_min = (F1_u1 >= F2_u1); 
u(u1_min) = u1(u1_min);
F(u1_min) = F1_u1(u1_min);
if nargout >= 3, indicators(u1_min) = 1; end

%%%%%%%%%%%%%%%%%%
% find the minimum of f2 %
%%%%%%%%%%%%%%%%%%

if any(~u1_min) % keep looking
    u2_min = ~u1_min;
    phase_ds = phase_d(u2_min) - phase_s(u2_min); % shift phase by phase_s
    
    % solve (I think this is correct now)
    if isscalar(etas), etas_u2 = etas; else etas_u2 = etas(u2_min); end
    aux1 = (2.*abs_s(u2_min).*(y(u2_min)+abs_s(u2_min).^2)+etas_u2.*abs_d(u2_min).*cos(phase_ds))./(etas_u2+4.*abs_s(u2_min).^2);
    aux2 = abs_d(u2_min).*sin(phase_ds);
    abs_u2 = sqrt(aux1.^2+aux2.^2);
    phase_u2 = atan2(aux2,aux1);
    u2 = abs_u2.*exp(complex(0,phase_u2+phase_s(u2_min)));
    F2_u2 = (y(u2_min)+abs_s(u2_min).^2-2.*real(conj(s(u2_min)).*u2)).^2+etas_u2.*abs(u2-d(u2_min)).^2;
    F1_u2 = (abs_u2.^2-y(u2_min)).^2+etas_u2.*abs(u2-d(u2_min)).^2;
    
    % solution iff u2 satisfies F2(u2) >= F1(u2)
    u2_test = (F2_u2 >= F1_u2);
    u2_min(u2_min) = u2_test;
    u(u2_min) = u2(u2_test);
    F(u2_min) = F2_u2(u2_test);
    if nargout >= 3, indicators(u2_min) = 2; end
else
    u2_min = false;
end

%%%%%%%%%%%%%%%%%%%%
% find the minimum of f1=f2 %
%%%%%%%%%%%%%%%%%%%%
if any(~u1_min & ~u2_min) % keep looking
    u3_min = ~u1_min & ~u2_min;
    phase_ds = phase_d(u3_min) - phase_s(u3_min); % shift phase by phase_s
    
    if isscalar(etas), etas_u3 = etas; else etas_u3 = etas(u3_min); end
    d_u3 = abs_d(u3_min).*exp(complex(0,phase_ds)); % shifted by phase_s
    y_u3 = y(u3_min);
    abs_s_u3 = abs_s(u3_min);
    c0 = sqrt(2.*y_u3+2.*abs_s_u3.^2);
    K1 = c0.^2+abs_s_u3.^2-y_u3;
    r1 = 2.*c0.*abs_s_u3;
    aux = 2.*c0.*(2.*K1.*abs_s_u3+etas_u3.*(abs_s_u3+d_u3));
    alpha = angle(aux);
    
    r1_nonzero = (r1 ~= 0);
    
    % trivial case (u = s)
    u3 = [abs_s_u3,NaN([length(abs_s_u3),4])];

    % for r1 == 0
    u3(:,2) = c0.*exp(complex(0,alpha))-abs_s_u3;

    if any(r1_nonzero)
        % assemble quartic using second expression from paper
        R = 1./(r1(r1_nonzero).^2); % divide by nonzero r1
        aa = (R.*imag(aux(r1_nonzero)));
        bb = (4+2.*R.*real(aux(r1_nonzero)));
        dd = bb-8;
        ee = -aa;

        xs = solve_quartic(aa,bb,0,dd,ee,'real');
        thetas = 2*atan(xs);
        u3_r1nz = bsxfun(@minus,bsxfun(@times,c0(r1_nonzero),exp(complex(0,thetas))),abs_s_u3(r1_nonzero));
        u3_r1nz(isnan(xs)) = NaN;
        u3(r1_nonzero,2:end) = u3_r1nz;
    end

    % find best u3 out of candidates
    F1_u3 = bsxfun(@minus,abs(u3).^2,y_u3).^2 + bsxfun(@times,etas_u3,abs(bsxfun(@minus,u3,d_u3)).^2);
    [F1_u3,iu3] = min(F1_u3,[],2);
    u3 = u3(sub2ind(size(u3),(1:size(u3,1)).',iu3));
    
    u(u3_min) = u3.*exp(complex(0,phase_s(u3_min)));
    F(u3_min) = F1_u3;
    if nargout >= 3, indicators(u3_min) = 3; end
end

u = reshape(u,sz_u);
if nargout >= 3, indicators = reshape(indicators,sz_u); end

end
