function [seed_dirs, tag] = l2f_runs(run_set)
% Where the network runs of the 8 l2f truths are, per seed (ECl, CA1, RSPd), and
% the suffix of the result files the analyses write for them.
%   "original"    the first runs: ECl in outputs_recovery/, CA1 and RSPd in
%                 outputs_seeds/ (joint flux solve; W1 at lambda2 = 0.01); no suffix
%   "consistent"  the runs in the surrogates' training convention, all seeds in
%                 outputs_consistent/ (per-edge flux solve for the numeric NTM;
%                 W1 at lambda2 = lambda1; solver_check/SOLVER_CHECK.md); suffix "_consistent"
%
%   [seed_dirs, tag] = l2f_runs("consistent");

switch string(run_set)
    case "original"
        seed_dirs = ["outputs_recovery", "outputs_seeds", "outputs_seeds"];
        tag = "";
    case "consistent"
        seed_dirs = ["outputs_consistent", "outputs_consistent", "outputs_consistent"];
        tag = "_consistent";
    otherwise
        error('run_set is "original" or "consistent", not "%s"', run_set);
end
end
