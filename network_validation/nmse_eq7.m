function v = nmse_eq7(s, n)
% NMSE as in Eq. 7 of the paper: the mean squared error of s against the
% reference n, divided by the variance of n, over all entries together
% (for network simulations: all regions and time points).
%   v = sum((s - n).^2) / sum((n - mean(n)).^2)
% Used for every error reported in the paper, the Supplement, the response and
% the parameter-recovery report (decided 2026-11-10; before, network and edge
% errors were sum((s - n).^2) / sum(n.^2), now sse_rel.m).
s = s(:);  n = n(:);
v = sum((s - n).^2) / sum((n - mean(n)).^2);
end
