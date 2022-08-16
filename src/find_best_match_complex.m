function [x_best,error] = find_best_match_complex(x,x_true)
% [x_best,error] = find_best_match_complex(x,x_true)
%
% Inputs:
%  x - reconstructed signal
%  x_true - true signal
%
% Outputs:
%  x_best - best version of reconstructed signal
%  error - RMS error of best version of x
%
% This function identifies the best shifted/flipped/phase shifted version
% of x and calculates its RMS error. It now supports 2D signals. The RMS
% error is normalized by the (square root) # of elements in x.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

% get correlations with x_true
sz_x = size(x_true);
x_true_fft = fftn(x_true);
x_fft = fftn(x);
[corr_x_best,corr_off] = max(col(abs(ifftn(x_true_fft.*conj(x_fft)))));
[corr_x_flipped,corr_off_flipped] = max(col(abs(ifftn(x_true_fft.*x_fft))));
corr_offs = cell(1,length(sz_x));
[corr_offs{:}] = ind2sub(sz_x,corr_off);
corr_offs_flipped = cell(1,length(sz_x));
[corr_offs_flipped{:}] = ind2sub(sz_x,corr_off_flipped);
if corr_x_best < corr_x_flipped
    x_best = conj(x);
    for idim = 1:length(sz_x), x_best = flip(x_best,idim); end
    x_best = circshift(x_best,[corr_offs_flipped{:}]); % flipped & shifted
else
    x_best = circshift(x,[corr_offs{:}]-1); % shifted
end

% find angle to correct
cphase = angle(x_best(:)'*x_true(:)); % based on least squares
x_best = x_best.*exp(complex(0,cphase));

error = norm(x_best(:) - x_true(:),2)./sqrt(numel(x_true));

end
