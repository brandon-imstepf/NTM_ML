function T = interp_training_data(cfg)
% Training samples of the nearest-neighbor (interp) sNTM: cfg.interp_train_size
% rows [gamma1 lambda1 delta epsilon N1 N2 phi_tau Delta_tau].
%   n <= 10^5  rows randperm(10^5, n) of cfg.train_superset_file after
%              rng(cfg.train_seed) (the GP's sampling; 10^4 before 2026-11-13)
%   n >  10^5  the MLP's training data, cfg.interp_train_files (10 x 10^5 rows,
%              data_bias_e5_1..10.csv, in that order); all of it for n = 10^6,
%              else rows randperm(10^6, n) after rng(cfg.train_seed)
% data_bias_e5_11.csv (the held-out set) is not among them.
%
%   T = interp_training_data(sntm_config());

persistent cache key
n = cfg.interp_train_size;
k = sprintf('%d|%s', n, strjoin(string(cfg.interp_train_files), ','));
if ~isempty(cache) && strcmp(key, k), T = cache; return, end
if n <= 1e5
    S = readmatrix(cfg.train_superset_file);
    rng(cfg.train_seed);
    T = S(randperm(size(S, 1), n), :);
else
    T = cell(numel(cfg.interp_train_files), 1);
    for i = 1:numel(T), T{i} = readmatrix(cfg.interp_train_files{i}); end
    T = vertcat(T{:});
    if n < size(T, 1)
        rng(cfg.train_seed);
        T = T(randperm(size(T, 1), n), :);
    end
end
cache = T;  key = k;
end
