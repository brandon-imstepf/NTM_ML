function make_locked_edge_data(sets_file, first_row, last_row, ncores)
% Single-edge datasets for parameter recovery.
%
% WHAT IT DOES
% Each row of sets_file "locks" the four model parameters (gamma1, lambda1,
% delta, epsilon). For that row, the script lays out a square grid of tau
% concentrations at the two ends of one edge (N1 = tau_x0, N2 = tau_xL) and,
% at every grid point, runs the numeric single-edge model to get:
%   F  (FValue)  steady-state flux along the edge
%   W1           edge mass term
% So one locked set = n_grid^2 simulations with the same parameters and
% different (N1, N2). With n_grid = 32 that is 1024 points.
%
% HOW TO CALL IT
%   make_locked_edge_data('locked_sets.csv', 1, 9)      % rows 1-9, one core
%   make_locked_edge_data('locked_sets.csv', 1, 3, 8)   % rows 1-3, 8 cores
%
%   sets_file  CSV with columns set_id, gamma1, lambda1, delta, epsilon, n_grid, N_max
%   first_row, last_row  which data rows of sets_file to run (header not counted)
%   ncores     optional. 0 or omitted = one point at a time. N > 0 = run the
%              grid points in parallel on N MATLAB workers (parfor). The pool
%              of workers stays open afterwards; close it with
%              delete(gcp('nocreate')).
%
% OUTPUT
% One file per locked set, <out_dir>/<set_id>.csv, with the same columns as
% the training data:
%   gamma1, lambda1, delta, epsilon, NRow (= N1), NCol (= N2), FValue, W1
% Simulation settings are the ones in NetworkTransportModel_SingleEdge_AllParams.m.
%
% Optional environment variables (mainly for the cluster):
%   NTM_SIM_ROOT    folder of the simulator (Tau_Transport-jtorok_dev)
%   LOCKED_OUT_DIR  where to write the CSVs (default: outputs/ next to this file)
%   LOCKED_SMOKE=1  quick test: 3 x 3 grid only, written to outputs_smoke/

%% 1. Tidy up the inputs
% On the cluster the arguments arrive as text (e.g. '1'), so convert to numbers.
if ischar(first_row) || isstring(first_row), first_row = str2double(first_row); end
if ischar(last_row) || isstring(last_row), last_row = str2double(last_row); end
if nargin < 4, ncores = 0; end % no 4th argument -> run serially
if ischar(ncores) || isstring(ncores), ncores = str2double(ncores); end

%% 2. Locate the simulator and the output folder
here = fileparts(mfilename('fullpath')); % folder this script lives in
% Default simulator location: laptop path on Windows, cluster path otherwise
if ispc
    default_sim_root = 'C:\Users\USER\Documents\NTM_ML Paper\Tau_Transport-jtorok_dev';
else
    default_sim_root = '/home/bimstepf/endor/vault/Tau_Transport-jtorok_dev';
end
sim_root = getenv_default('NTM_SIM_ROOT', default_sim_root);

% Smoke-test mode writes to a separate folder so it never mixes with real results
smoke = strcmp(getenv('LOCKED_SMOKE'), '1');
if smoke, default_out = fullfile(here, 'outputs_smoke'); else, default_out = fullfile(here, 'outputs'); end
out_dir = getenv_default('LOCKED_OUT_DIR', default_out);
if ~isfolder(out_dir), mkdir(out_dir); end

% Make the simulator functions callable, and point at its data files
addpath(fullfile(sim_root, 'Code', 'Transport_Network'));
matdir = fullfile(sim_root, 'MatFiles');

%% 3. Settings that never change between runs
% Name/value pairs passed straight to the simulator. Same values as
% NetworkTransportModel_SingleEdge_AllParams.m (tight 1e-10 tolerances, as
% used to make the training data).
fixed = {'beta', 1e-6, 'gamma2', 0, 'frac', 0.7, ...
    'L_int', 1000, 'L1', 200, 'L2', 200, 'L_ais', 40, 'L_syn', 40, ...
    'resmesh', 'coarse', 'reltol', 1e-10, 'abstol', 1e-10, 'fsolvetol', 1e-10};

%% 4. Loop over the locked parameter sets
S = readtable(sets_file, 'TextType', 'string', 'Delimiter', ',');
last_row = min(last_row, height(S)); % don't run past the end of the file
for r = first_row:last_row
    s = S(r, :); % this row's set_id, gamma1, lambda1, delta, epsilon, n_grid, N_max

    % 4a. Skip sets that are already finished, so a rerun only does what's missing
    out_file = fullfile(out_dir, s.set_id + ".csv");
    if isfile(out_file)
        fprintf('[%s] exists, skipping\n', s.set_id);
        continue
    end

    % 4b. Build the (N1, N2) grid.
    % Nv = n_grid evenly spaced values from 0 to N_max. ndgrid makes every
    % (N1, N2) combination of them; (:) flattens each into a column, so row k
    % of [N1, N2] is the k-th grid point.
    n_grid = s.n_grid; if smoke, n_grid = 3; end
    Nv = linspace(0, s.N_max, n_grid);
    [N1, N2] = ndgrid(Nv, Nv);
    N1 = N1(:); N2 = N2(:);
    np = numel(N1); % number of grid points = n_grid^2
    fprintf('[%s] gamma1=%g lambda1=%g delta=%g epsilon=%g  %d points, N in [0, %g]\n', ...
        s.set_id, s.gamma1, s.lambda1, s.delta, s.epsilon, np, s.N_max);

    % 4c. Full argument list for the simulator: fixed settings + this set's
    % locked parameters (lambda2 is kept equal to lambda1, as in training)
    a = [fixed, {'gamma1', s.gamma1, 'lambda1', s.lambda1, 'lambda2', s.lambda1, ...
        'delta', s.delta, 'epsilon', s.epsilon}];

    % 4d. Run the simulator at every grid point.
    % F and W1 start as NaN; a point that fails stays NaN.
    [F, W1] = deal(nan(np, 1));
    set_id = s.set_id; % plain variable, so parfor workers can use it
    t0 = tic;
    if ncores > 0
        % --- Parallel: grid points are shared out across ncores workers ---
        % Start a pool of workers if there isn't one of the right size yet
        pool = gcp('nocreate');
        if isempty(pool) || pool.NumWorkers ~= ncores
            delete(pool);
            parpool(ncores);
        end
        % Workers can't print progress in order, so each one sends a message
        % to this DataQueue when it finishes a point, and the main MATLAB
        % session counts them and prints every 100.
        q = parallel.pool.DataQueue;
        progress(set_id, np, t0, true); % reset the counter to 0
        afterEach(q, @(~) progress(set_id, np, t0, false));
        parfor k = 1:np
            [F(k), W1(k)] = edge_point(matdir, a, N1(k), N2(k), k, set_id);
            send(q, k); % "point k done"
        end
    else
        % --- Serial: one point after another ---
        progress(set_id, np, t0, true);
        for k = 1:np
            [F(k), W1(k)] = edge_point(matdir, a, N1(k), N2(k), k, set_id);
            progress(set_id, np, t0, false);
        end
    end

    % 4e. Save in the training-data column layout.
    % The locked parameters are repeated on every row (repmat).
    T = table(repmat(s.gamma1, np, 1), repmat(s.lambda1, np, 1), repmat(s.delta, np, 1), ...
        repmat(s.epsilon, np, 1), N1, N2, F, W1, 'VariableNames', ...
        {'gamma1', 'lambda1', 'delta', 'epsilon', 'NRow', 'NCol', 'FValue', 'W1'});
    % Write to <set_id>.csv.part first, then rename. If MATLAB is killed
    % mid-write, there is no <set_id>.csv, so step 4a won't wrongly skip the set.
    writetable(T, out_file + ".part", 'FileType', 'text');
    movefile(out_file + ".part", out_file);
    fprintf('[%s] done in %.1f s, %d failed points\n', s.set_id, toc(t0), sum(isnan(F) | isnan(W1)));
end
end

%% Helper functions

function [F, W1] = edge_point(matdir, a, n1, n2, k, set_id)
% One simulation: flux F and edge mass W1 at a single (N1, N2) = (n1, n2).
% Two calls, same as the training-data script:
%   1. SingleEdgeFluxCalculator_WrapperFun -> steady-state flux F
%   2. SingleEdgeWCalculator, given that F  -> W1
% use_sr_* = 0 means the numeric model is used, not a surrogate.
% If either call errors, the message is printed and NaN is returned.
F = nan; W1 = nan;
try
    o = SingleEdgeFluxCalculator_WrapperFun(matdir, a{:}, ...
        'tau_x0', n1, 'tau_xL', n2, 'sim_no', k, 'use_sr_flux', 0, 'use_sr_w1', 0);
    F = o.Predicted.F;
    W1 = SingleEdgeWCalculator(matdir, F, a{:}, 'tau_x0', n1, 'tau_xL', n2, 'use_sr_w1', 0);
catch err
    fprintf('[%s] point %d (N1=%g, N2=%g) failed: %s\n', set_id, k, n1, n2, err.message);
end
end

function progress(set_id, np, t0, reset)
% Progress counter. reset = true sets it to 0; otherwise it adds one finished
% point and prints a line every 100 points and at the last one.
% "persistent" means n keeps its value between calls.
persistent n
if reset || isempty(n), n = 0; end
if reset, return, end
n = n + 1;
if mod(n, 100) == 0 || n == np
    fprintf('[%s] %d/%d points, %.0f s\n', set_id, n, np, toc(t0));
end
end

function v = getenv_default(name, default)
% Value of environment variable "name", or "default" if it isn't set
v = getenv(name);
if isempty(v), v = default; end
end
