function T = time_single_edge(n, nworkers)
% Time the numeric single-edge data generator (training data of the surrogates),
% for the cost of generating 10^6 samples (collaborator, 11/20). Each sample is
% computed exactly as in the parfor branch of the simulator's
% NetworkTransportModel_SingleEdge_AllParams.m: the same input sampling (rng(42),
% gamma1 and lambda1 log-uniform, delta and epsilon uniform, N1 and N2 Gamma(2, 5e-4)
% as in Table 7), the same tolerances (1e-10) and options, and the same two calls,
% SingleEdgeFluxCalculator_WrapperFun (phi_tau) and SingleEdgeWCalculator
% (Delta_tau). The first n of the generator's samples are timed.
%   nworkers = 0: one after the other on one computational thread.
%   nworkers = p: a parfor over the n samples on a pool of p workers, as the
%     generator (it uses parpool(4)); each worker has one thread.
% Reports the time per sample (summed over workers: the core time) and the
% wall-clock time of the n samples (without starting the pool), and extrapolates
% both to 10^6 samples. Writes results/time_single_edge_<host>_<p>w.csv (one row
% per sample) and prints the summary.
%
%   T = time_single_edge();          % 200 samples, one thread
%   T = time_single_edge(400, 4);    % 400 samples, parfor on 4 workers
% On the cluster: sbatch time_single_edge.slurm (see its header).

if nargin < 1, n = 200; end
if nargin < 2, nworkers = 0; end
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
loadpath = cfg.matdir;

% sampling of NetworkTransportModel_SingleEdge_AllParams.m (first n of its rows)
rng(42);
N = 100000;
logdist = @(x, m) 10.^(log10(x(1)) + (log10(x(2) / x(1)) * rand(m, 1)));
lindist = @(x, m) x(1) + (x(2) - x(1)) * rand(m, 1);
% N1, N2 ~ Gamma(2, 5e-4) as the training data (main text Table 7; the generator
% in MLP_Tau_Transport, "biased widerange"), drawn without the Statistics Toolbox
% as -5e-4*log(U1*U2). (Before 11/21: uniform on [0, 0.005], as the copy of the
% generator in the simulator folder.)
gamdist = @(m) -5e-4 * log(rand(m, 1) .* rand(m, 1));
P = [logdist([5e-5 5e-2], N), logdist([0.01 0.1], N), lindist([10 100], N), lindist([10 100], N), ...
     gamdist(N), gamdist(N)];
P = P(1:n, :);
names = {'gamma1', 'lambda1', 'delta', 'epsilon', 'tau_x0', 'tau_xL'};
opt = {'beta', 1e-6, 'frac', 0.7, 'L_int', 1000, 'L1', 200, 'L2', 200, 'L_ais', 40, 'L_syn', 40, ...
       'resmesh', 'coarse', 'reltol', 1e-10, 'abstol', 1e-10, 'fsolvetol', 1e-10};

t_flux = zeros(n, 1);  t_w1 = zeros(n, 1);
if nworkers > 0
    % a pool with its own job storage folder, so that simultaneous jobs on a
    % cluster node do not share one
    % (the folder must exist before it is assigned; on a SLURM node use the job's
    % own scratch folder if there is one)
    jsl = tempname(getenv_or('TMPDIR', tempdir));
    [ok, msg] = mkdir(jsl);
    assert(ok, 'cannot create the job storage folder %s: %s', jsl, msg);
    c = parcluster('local');
    c.JobStorageLocation = jsl;
    pool = parpool(c, nworkers);
    t_wall = tic;
    parfor i = 1:n
        a = [names; num2cell(P(i, :))];
        t0 = tic;
        out = SingleEdgeFluxCalculator_WrapperFun(loadpath, a{:}, opt{:}, 'sim_no', i, 'use_sr_flux', 0, 'use_sr_w1', 0);
        t_flux(i) = toc(t0);
        t0 = tic;
        SingleEdgeWCalculator(loadpath, out.Predicted.F, a{:}, opt{:}, 'lambda2', P(i, 2), 'gamma2', 0);
        t_w1(i) = toc(t0);
    end
    wall = toc(t_wall);
    delete(pool);
    rmdir(jsl, 's');
else
    maxNumCompThreads(1);
    t_wall = tic;
    for i = 1:n
        a = [names; num2cell(P(i, :))];
        t0 = tic;
        out = SingleEdgeFluxCalculator_WrapperFun(loadpath, a{:}, opt{:}, 'sim_no', i, 'use_sr_flux', 0, 'use_sr_w1', 0);
        t_flux(i) = toc(t0);
        t0 = tic;
        SingleEdgeWCalculator(loadpath, out.Predicted.F, a{:}, opt{:}, 'lambda2', P(i, 2), 'gamma2', 0);
        t_w1(i) = toc(t0);
    end
    wall = toc(t_wall);
end
total = t_flux + t_w1;
[~, host] = system('hostname');  host = strtrim(host);
T = array2table([P, t_flux, t_w1, total], 'VariableNames', [names, {'t_flux_s', 't_w1_s', 't_total_s'}]);
writetable(T, fullfile(here, 'results', sprintf('time_single_edge_%s_%dw.csv', host, nworkers)));
p = max(nworkers, 1);
fprintf('\n%s, MATLAB %s, %d worker(s), %d samples\n', host, version, p, n);
fprintf('per sample (one core): mean %.3g s (phi_tau %.3g, Delta_tau %.3g), median %.3g, max %.3g\n', ...
    mean(total), mean(t_flux), mean(t_w1), median(total), max(total));
fprintf('wall clock for the %d samples: %.1f s (%.3g s per sample)\n', n, wall, wall / n);
fprintf('10^6 samples: %.3g core-hours; wall clock with %d worker(s) %.3g h (%.3g days)\n', ...
    1e6 * mean(total) / 3600, p, 1e6 * wall / n / 3600, 1e6 * wall / n / 86400);
end

function v = getenv_or(name, default)
v = getenv(name);
if isempty(v) || ~isfolder(v), v = default; end
end
