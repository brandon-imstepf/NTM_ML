function run_local(manifest_file, n_workers)
% Run every row of a manifest on this machine, one case per parallel worker
% (Parallel Computing Toolbox). Finished cases are skipped, so the same call
% resumes an interrupted run. A case that errors is logged and the rest continue.
%
%   run_local('manifests/surrogate.csv')
%   run_local('manifests/numeric.csv', 6)

if nargin < 2, n_workers = feature('numcores'); end
M = readtable(manifest_file, 'TextType', 'string', 'Delimiter', ',');
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= n_workers
    delete(pool);
    pool = parpool(n_workers);
end
% Hand out one case at a time so long cases do not leave workers idle
opts = parforOptions(pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
failed = strings(height(M), 1);
parfor (r = 1:height(M), opts)
    try
        run_sntm_cases(manifest_file, r, r);
    catch err
        failed(r) = M.case_id(r) + ": " + err.message;
        fprintf(2, '[%s] FAILED: %s\n', M.case_id(r), err.message);
    end
end
failed = failed(failed ~= "");
fprintf('%d of %d cases failed\n', numel(failed), height(M));
fprintf('  %s\n', failed);
end
