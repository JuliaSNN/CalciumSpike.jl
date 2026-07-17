function y = gaussianfilter(x, dim, sigma, forward)
%GAUSSIANFILTER Gaussian smoothing filter, symmetric or forward (causal), for a 1-D signal.
%
% y = gaussianfilter(x, dim, sigma, forward)
%
% - x        1-D data (row or column vector)
% - dim      accepted for call-signature compatibility with the pipeline
%            that calls this function; not used — the vector's own
%            orientation (row/column) is auto-detected, same convention as
%            this toolbox's own fn_smooth.m
% - sigma    Gaussian standard deviation, in samples
% - forward  true/1: causal filtering (only past+current samples
%            contribute). false/0: symmetric (past+future), same shape as
%            fn_smooth.m
%
% Edge-renormalized: near the boundaries the kernel is truncated and
% renormalized by the locally-available weight (no zero-padding bias) —
% same principle as fn_smooth.m and CalciumSpike.jl's gaussian_smooth.
%
% Reconstructed: the original `gaussianfilter` this pipeline calls was not
% found (checked this repo's full git history and searched the web — it is
% not a standard MATLAB/toolbox function); this is a from-scratch
% implementation matching its call signature and the forward-mode
% semantics implied by its 4th argument.

if nargin < 4
    forward = false;
end

n = length(x);
half = ceil(4*sigma);
if forward
    offsets = -half:0;
else
    offsets = -half:half;
end
kernel = exp(-offsets.^2 / (2*sigma^2));
kernel = kernel / sum(kernel);

y = zeros(size(x));
for i = 1:n
    s = 0; w = 0;
    for j = 1:length(offsets)
        idx = i + offsets(j);
        if idx < 1 || idx > n
            continue
        end
        s = s + kernel(j) * x(idx);
        w = w + kernel(j);
    end
    y(i) = s / w;
end

end
