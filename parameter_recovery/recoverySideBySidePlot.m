function recoverySideBySidePlot(outdir)
% For every locked set, the phi_tau (flux) and Delta_tau (w1) recoveries in one
% figure: each model gets two columns next to each other, loss on phi_tau and
% loss on Delta_tau (recoveryOverviewPlot.m with both results).
%
% Reads <outdir>/flux/recovery_curves.mat and <outdir>/w1/recovery_curves.mat
% (synthetic_parameter_recovery_all_models.m). Writes to <outdir>/both/by_set/:
%   single_<set>.png   rows: parameters; columns: model x output
%   joint_<set>.png    rows: parameter pairs; columns: model x output
%
%   recoverySideBySidePlot('recovery_all_models')

    F = load(fullfile(outdir, 'flux', 'recovery_curves.mat'));
    W = load(fullfile(outdir, 'w1', 'recovery_curves.mat'));
    assert(isequal(F.rc.set_ids, W.rc.set_ids), 'flux and w1 were run on different locked sets');
    both_dir = fullfile(outdir, 'both');
    for s = 1:numel(F.rc.set_ids)
        recoveryOverviewPlot({F.rc, W.rc}, s, both_dir);
    end
    fprintf('wrote %d sets to %s\n', numel(F.rc.set_ids), fullfile(both_dir, 'by_set'));
end
