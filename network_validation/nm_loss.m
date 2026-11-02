function loss = nm_loss(N, N_truth, kind, day)
% Loss of a simulated N (regions x time points) against the truth's N, for
% nelder_mead_network_recovery.m (NM_LOSS) and compare_losses_t60.m.
%   nmse    sum (s - n)^2 / sum n^2 over all regions and times (the default)
%   region  the same per region, then averaged over regions: every region
%           counts equally, not only the ones with the most tau
%   log     mean of (log10(s + c) - log10(n + c))^2, c = 1e-3 x the truth's
%           largest N: relative error, with c keeping near-zero values from
%           dominating
%   late    nmse over day >= 20 only (leaves out the early transient)
% day  time points in days (needed for late)
%
%   loss = nm_loss(out.Predicted.N, truth.N, "region")

if nargin < 3 || isempty(kind), kind = "nmse"; end
switch string(kind)
    case "nmse"
        loss = sum((N(:) - N_truth(:)).^2) / sum(N_truth(:).^2);
    case "region"
        r = sum(N_truth.^2, 2) > 0;   % regions the truth never reaches have no scale
        loss = mean(sum((N(r, :) - N_truth(r, :)).^2, 2) ./ sum(N_truth(r, :).^2, 2));
    case "log"
        c = 1e-3 * max(N_truth(:));
        loss = mean((log10(max(N(:), 0) + c) - log10(N_truth(:) + c)).^2);
    case "late"
        k = day >= 20;
        loss = sum((N(:, k) - N_truth(:, k)).^2, 'all') / sum(N_truth(:, k).^2, 'all');
    otherwise
        error('unknown loss %s: nmse, region, log or late', kind);
end
end
