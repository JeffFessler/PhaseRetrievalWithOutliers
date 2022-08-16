function x = shrink(x,tau)
% perform shrinkage for soft-thresholding
%
% Input:
%  x -- vector to shrink
%  tau -- shrinkage threshold parameter
%
% Output:
%  x -- vector after shrinkage
%
% Copyright 2014-01-31 Daniel Weller, University of Michigan. All rights
% reserved.

xabs = abs(x);
x = x.*(1-tau./xabs);
x(xabs <= tau) = 0;

end
