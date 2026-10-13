function cfg = sntm_config()
% Paths and fixed simulation settings shared by every network-validation run.
% Paths can be overridden with the environment variables NTM_ML_ROOT,
% NTM_SIM_ROOT and SNTM_OUT_DIR (useful on the cluster).

here = fileparts(mfilename('fullpath'));
cfg.ntm_ml_root = getenv_default('NTM_ML_ROOT', fileparts(here));

if ispc
    default_sim_root = 'C:\Users\USER\Documents\NTM_ML Paper\Tau_Transport-jtorok_dev';
else
    default_sim_root = '/home/bimstepf/endor/vault/Tau_Transport-jtorok_dev';
end
cfg.sim_root = getenv_default('NTM_SIM_ROOT', default_sim_root);
cfg.matdir = fullfile(cfg.sim_root, 'MatFiles');
cfg.out_dir = getenv_default('SNTM_OUT_DIR', fullfile(here, 'outputs'));

% Surrogate artifacts
cfg.sr_dir = fullfile(cfg.ntm_ml_root, 'parameter_recovery');
cfg.mlp_weights_file = fullfile(here, 'models', 'mlp_e6_weights.mat'); % from export_mlp_weights.py
cfg.lr_weights_file = fullfile(here, 'models', 'lr_weights.mat'); % from fit_lr_surrogate.py
% Compact copy of gaussian_process-interpolant/saved_models_superset_5000/
% FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc.mat (identical predictive
% mean), with the target scaling of its 5000-sample training subset stored inside
cfg.gp_model_file = fullfile(here, 'models', 'gp_e5k_exact_linear_ardmatern32_fminunc.mat');
cfg.interp_train_size = 10000; % must match interpolant.m
cfg.train_superset_file = fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', ...
    'training_data_e6', 'data_bias_e5_1.csv');
cfg.train_seed = 42;

% Fixed NTM settings (trange: 151 points, t = 1 is 180 days)
cfg.beta = 1e-6;
cfg.gamma2 = 0;
cfg.frac = 0.7;
cfg.L_int = 1000;
cfg.L1 = 200;
cfg.L2 = 200;
cfg.L_ais = 40;
cfg.L_syn = 40;
cfg.trange = [0:0.0025:0.1, 0.105:0.005:0.3, 0.31:0.01:1];
if ~isempty(getenv('SNTM_SMOKE')) % short horizon for pipeline tests
    cfg.trange = 0:0.0025:0.01;
end
cfg.resmesh = 'coarse';
% ODE tolerances; SNTM_RELTOL / SNTM_ABSTOL override them for settings where the
% steady-state solve fails at the defaults (each output stores the cfg it used)
cfg.reltol = str2double(getenv_default('SNTM_RELTOL', '1e-4'));
cfg.abstol = str2double(getenv_default('SNTM_ABSTOL', '1e-4'));
cfg.fsolvetol = 1e-6;
cfg.init_rescale = 1e-2;
cfg.study = 'Hurtado';
cfg.conn_thresh = 0.8;

% Seed codes used in the manifest -> CCF region names expected by NetworkTransportModel
cfg.seeds = containers.Map( ...
    {'ECl', 'CA1', 'DG', 'SUBd', 'RSPd'}, ...
    {'Entorhinal area, lateral part_L', 'Field CA1_L', 'Dentate gyrus_L', ...
     'Subiculum, dorsal part_L', 'Retrosplenial area, dorsal part_L'});
end

function val = getenv_default(name, default)
val = getenv(name);
if isempty(val)
    val = default;
end
end
