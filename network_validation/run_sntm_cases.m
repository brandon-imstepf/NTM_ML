function run_sntm_cases(manifest_file, first_row, last_row)
% Run rows first_row..last_row of a manifest CSV (see make_manifest.py).
% Each case is saved to <out_dir>/<case_id>.mat; finished cases are skipped,
% so a job can simply be resubmitted after a timeout.

if ischar(first_row) || isstring(first_row), first_row = str2double(first_row); end
if ischar(last_row) || isstring(last_row), last_row = str2double(last_row); end

cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
addpath(fileparts(mfilename('fullpath')));
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end

M = readtable(manifest_file, 'TextType', 'string', 'Delimiter', ',');
last_row = min(last_row, height(M));
surrogates = struct();

for r = first_row:last_row
    c = M(r, :);
    out_file = fullfile(cfg.out_dir, c.case_id + ".mat");
    if isfile(out_file)
        fprintf('[%s] exists, skipping\n', c.case_id);
        continue
    end

    method = char(c.method);
    if ~isfield(surrogates, method)
        surrogates.(method) = make_surrogate(method, cfg);
    end
    sur = surrogates.(method);
    use_sur = ~strcmp(method, 'numeric');
    % Per-case seed load from the manifest when given, else the config default
    init_rescale = cfg.init_rescale;
    if ismember('init_rescale', M.Properties.VariableNames), init_rescale = c.init_rescale; end
    % lambda2 from the manifest when given, else lambda2 = lambda1
    lambda2 = c.lambda1;
    if ismember('lambda2', M.Properties.VariableNames), lambda2 = c.lambda2; end

    fprintf('[%s] method=%s seed=%s gamma=%g lambda=%g lambda2=%g delta=%g epsilon=%g init_rescale=%g subset=%s\n', ...
        c.case_id, method, c.seed, c.gamma1, c.lambda1, lambda2, c.delta, c.epsilon, init_rescale, c.subset);
    t0 = tic;
    p = struct('gamma1', c.gamma1, 'lambda1', c.lambda1, 'lambda2', lambda2, 'delta', c.delta, 'epsilon', c.epsilon);
    out = sntm_simulate(cfg, sur, use_sur, p, c.seed, c.subset, init_rescale, r);
    runtime_s = toc(t0);

    res = struct();
    res.meta = table2struct(c);
    res.N = out.Predicted.N;
    res.M = out.Predicted.M;
    res.F = out.Predicted.F;
    res.W1 = out.Predicted.W1;
    if ~use_sur
        res.EdgeMass = out.Predicted.EdgeMass;
    end
    res.trange = out.Sim.trange;
    res.region_names = out.Sim.region_names;
    res.C = out.Sim.C;
    res.runtime_s = runtime_s;
    res.env = run_environment();
    res.cfg = rmfield(cfg, 'seeds');
    res.cfg.init_rescale = init_rescale;
    res.cfg.lambda2 = lambda2;
    save(out_file, '-struct', 'res', '-v7.3');
    fprintf('[%s] done in %.1f s\n', c.case_id, runtime_s);
end
end

function env = run_environment()
[~, host] = system('hostname');
env.host = strtrim(host);
env.matlab = version;
env.ncores_visible = feature('numcores');
env.max_threads = maxNumCompThreads;
env.slurm_job = getenv('SLURM_JOB_ID');
env.slurm_task = getenv('SLURM_ARRAY_TASK_ID');
if isunix
    [~, cpu] = system('grep -m1 "model name" /proc/cpuinfo');
    env.cpu = strtrim(cpu);
else
    env.cpu = getenv('PROCESSOR_IDENTIFIER');
end
env.timestamp = string(datetime('now'));
end
