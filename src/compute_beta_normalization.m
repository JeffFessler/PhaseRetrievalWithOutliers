function norm_factor = compute_beta_normalization(y,outliers,outliers_range,AWGN_SNR,AWLN_SNR,p)
% norm_factor = compute_beta_normalization(y,outliers,outliers_range,AWGN_SNR,AWLN_SNR,p)
%
% Inputs:
%  y - measurements
%  outliers - approx. # of outliers
%  outliers_range - approx. range of outliers (unused)
%  AWGN_SNR,AWLN_SNR - approx. SNR (dB) of additive Gaussian or Laplace noise
%  p - data fit term power (p = 1 or 2 typically)
%
% Outputs:
%  norm_factor - normalization factor for data fit term
%
% This function computes a normalization factor for the data fit term based
% on the l2-norm of its gradient. The function assumes the A matrix is
% right unitary.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

numy = numel(y);
if outliers > 0
    y = sort(y(:),'descend');
    yout = y(1:outliers);
    y = y(outliers+1:end);
else
    yout = [];
    y = y(:);
end

if p == 1
    % f(x) = sum_m sqrt((yt_m-y_m)^2+eps)-eps+beta*sum_n sqrt(x_n^2+eps)-eps
    % df(x) = sum_m (2a_m'a_mx)(yt_m-y_m)/sqrt((yt_m-y_m)^2+eps)) + beta*[x_n/sqrt(x_n^2+eps)]
    % df(x) = 2A'*diag((yt-y)./sqrt((yt-y).^2+eps))*Ax + beta*diag(1/sqrt(x.^2+eps))*x
    % norm_factor = ||2A'*diag((yt-y)./sqrt((yt-y).^2+eps))*Ax||_2
    % norm_factor^2 = 4*\sum_m |sign(yt_m-y_m).*a_mx|^2 = 4*\sum(yt_m,has noise)
    % norm_factor^2 = 4*(outliers*avgy + (AWGN or AWLN).*sum(y_m,not outliers))
    avgy = mean(y);
    outliers_noise = outliers*avgy;
    if isfinite(AWGN_SNR) || isfinite(AWLN_SNR)
        regular_noise = sum(y);
    else
        regular_noise = 0;
    end
    norm_factor = 2*sqrt(outliers_noise + regular_noise);
elseif p == 2
    % f(x) = sum_m (yt_m-y_m)^2+beta*sum_n sqrt(x_n^2+eps)-eps
    % df(x) = sum_m 2(2a_m'a_mx)(yt_m-y_m) + beta*[x_n/sqrt(x_n^2+eps)]
    % df(x) = 4A'*diag(yt-y)*Ax + beta*diag(1/sqrt(x.^2+eps))*x
    % norm_factor = ||4*A'*diag(yt-y)*Ax||_2
    % norm_factor^2 = 16*\sum_m (yt_m-y_m)^2*|a_mx|^2
    % norm_factor^2 = 16*(noise_var*sum(y_m,not outliers) + sum((yout-yavg).^2.*yavg))
    
    avgy = mean(y);
    avgy2 = mean(y.^2);
    noise_var = (10^(-AWGN_SNR/10)+10^(-AWLN_SNR/10))*avgy2;
    norm_factor = 4*sqrt((numy-outliers)*noise_var*avgy + sum((yout-avgy).^2)*avgy);
else
    fail('Unsupported choice of p');
end

if norm_factor == 0, norm_factor = 1; end

end
