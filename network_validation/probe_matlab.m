% MATLAB-side checks for probe_cluster.slurm: licence checkout, toolboxes,
% and the functions the NTM simulations and surrogates call.
fprintf('version: %s\n', version);
v = ver;
fprintf('installed: %s\n', strjoin({v.Name}, ' | '));
fprintf('licence test  Optimization: %d  Statistics: %d  Parallel: %d  DeepLearning: %d\n', ...
    license('test', 'Optimization_Toolbox'), license('test', 'Statistics_Toolbox'), ...
    license('test', 'Distrib_Computing_Toolbox'), license('test', 'Neural_Network_Toolbox'));
for f = {'fsolve', 'ode45', 'deval', 'readtable', 'readmatrix', 'fitrgp', 'KDTreeSearcher', 'parpool'}
    fprintf('  %-15s %s\n', f{1}, which(f{1}));
end
try
    x = fsolve(@(x) x.^2 - 2, 1, optimset('Display', 'off'));
    fprintf('fsolve works: sqrt(2) = %.6f\n', x);
catch e
    fprintf('fsolve FAILED: %s\n', e.message);
end
try
    [~, y] = ode45(@(t, y) -y, [0 1], 1);
    fprintf('ode45 works: exp(-1) = %.6f\n', y(end));
catch e
    fprintf('ode45 FAILED: %s\n', e.message);
end
try
    m = fitrgp((1:20)', sin((1:20)'));
    fprintf('Statistics Toolbox checkout works: fitrgp ok (%s)\n', class(m));
catch e
    fprintf('Statistics Toolbox FAILED: %s\n', e.message);
end
fprintf('cores: %d  threads: %d\n', feature('numcores'), maxNumCompThreads);
