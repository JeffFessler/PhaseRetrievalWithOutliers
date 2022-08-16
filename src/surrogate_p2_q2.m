function f = surrogate_p2_q2(x,s,Ax,y,weights,beta,check_use_y_s)
% f = surrogate_p2_q2(x,s,Ax,y,weights,beta,check_use_y_s)
%
% Inputs:
%  x - signal value
%  s - majorization vector
%  Ax - measurement matrix (A) times x (A*x)
%  y - measurements
%  weights - weights for data fit term
%  beta - regularization parameter
%  check_use_y_s - if true, will use sqrt(y)*exp(i*angle(s)) if abs(s) > sqrt(y) (default is false)
%
% Outputs:
%  f - majorizer value for x
%
% This function computes the value of the majorizer function for quadratic
% data fit term with squared-magnitude measurements.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

hplus = @(t,y) abs(t).^2-y;
phiminus = @(t,s,y) y+abs(s).^2-2.*real(conj(s).*t);
phi = @(t,s,y) max(hplus(t,y),phiminus(t,s,y));

yleq0 = (y <= 0); % use |t|^2-y only
if exist('check_use_y_s','var') && check_use_y_s
    % tighter majorizer for those |s|^2 > y
    use_y_s = (~yleq0 & abs(s).^2 > y);
    s(use_y_s) = sqrt(y(use_y_s)).*exp(complex(0,angle(s(use_y_s))));
end

t = Ax;
if isempty(weights), weights = 1; end
f = sum(col(weights.*(yleq0.*hplus(t,y).^2+(~yleq0).*phi(t,s,y).^2))) + beta*sum(abs(col(x)));

end
