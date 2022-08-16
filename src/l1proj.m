function [x,theta] = l1proj(y,K,w)
% [x,theta] = l1proj(y,K,w)
%
% Inputs:
%  y - signal to project
%  K - radius of l1-ball
%  w - weights (optional), constrains w'*abs(y) <= K
%
% Outputs:
%  x - projected signal
%  theta - learned threshold for soft threshold operation that yields x
%
% This function performs l1 projection based on sorting, which is a
% standard O(N log N) method. Note that other methods are faster on
% average; this method is considered for its simplicity.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

if ~exist('w','var'), w = []; end

yabs = abs(y);
N = length(yabs);
if ~isempty(w)
    if w'*yabs <= K, x = y; return; end
else
    if sum(yabs) <= K, x = y; return; end
end

if ~isempty(w)
    [~,sortinds] = sort(yabs./w,'descend');
    ysort = yabs(sortinds);
    wsort = w(sortinds);
    wsortsqd = cumsum(wsort.^2);
    s = cumsum(wsort.*ysort);
    rho = find(ysort.*wsortsqd > (s - K).*wsort, 1, 'last');
else
    ysort = sort(yabs,'descend');
    wsortsqd = (1:N).';
    s = cumsum(ysort);
    rho = find(ysort.*wsortsqd > s - K, 1, 'last');
end
theta = max(0, (s(rho) - K) / wsortsqd(rho));

if ~isempty(w)
    x = max(yabs-theta.*w,0);
else
    x = max(yabs-theta,0);
end
x = x.*(y./yabs);
x(isnan(x)) = 0;

end
