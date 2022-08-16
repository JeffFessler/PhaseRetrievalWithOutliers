function disc = make_disc(radius,value)
% disc = make_disc(radius,value)
%
% Inputs:
%  radius - radius of disc in pixels (disc is 2*radius pixels wide)
%  value - value of pixels in the disc
%
% Outputs:
%  disc - small image containing disc
%
% This function makes a disc of the specified radius, filled with the
% specified value. Note: around the border, there's a sharp edge, not a
% soft one.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

if ~exist('value','var') || isempty(value), value = 1; end

% make disc
disc = zeros(2*radius);

[x,y] = ndgrid((0:2*radius-1)-(radius-0.5),(0:2*radius-1)-(radius-0.5));
distsqds = x.^2 + y.^2;
disc(distsqds <= radius.^2) = value;

end
