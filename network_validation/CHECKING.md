# Checking cluster runs

How to check that runs finished, copy them back, and compare two runs.
Commands marked **cluster** run in PuTTY from
`/home/bimstepf/endor/vault/NTM_ML/network_validation`. Commands marked **PC**
run in PowerShell or MATLAB from `C:\Users\USER\Documents\Code\NTM_ML\network_validation`.

## 1. Are the jobs done? (cluster)

```bash
squeue -u bimstepf                       # anything still queued or running
sacct -j 4700 -X --format=JobID%15,State,ExitCode,Elapsed,NodeList
```

`-X` shows one line per array task. Without it, `sacct` also prints every job
step (`.batch`, `.extern`), which is where the hundreds of lines come from.
`State` should be `COMPLETED` with `ExitCode 0:0`. `FAILED`, `TIMEOUT`,
`OUT_OF_MEMORY` or `CANCELLED` mean that case has no output.

Only the tasks that did not complete:

```bash
sacct -j 4700 -X --format=JobID%15,State,ExitCode -n | grep -v COMPLETED
```

## 2. Did every case write its file? (cluster, and again on the PC after copying)

Each manifest row writes `<SNTM_OUT_DIR>/<case_id>.mat`. Count, then list
the missing ones:

```bash
ls outputs_t60_dt0.45/*.mat | wc -l      # should equal the manifest rows (16 here)
tail -n +2 manifests/recovery_truth_mlp_l2f.csv | cut -d, -f1 | while read c; do
  [ -f outputs_t60_dt0.45/$c.mat ] || echo "missing $c"; done
```

`run_sntm_cases` skips a case whose file already exists. A file left over from
an earlier run with other settings is therefore never overwritten. Use a new
`SNTM_OUT_DIR` for every new setting, or move the old folder away first.

## 3. Logs (cluster)

Each array task writes `logs/sntm_<jobid>_<task>.out`.

```bash
grep -il error logs/sntm_4700_*.out      # tasks that printed an error
grep -h "done in" logs/sntm_4700_*.out   # runtime of each case
grep -h "^host" logs/sntm_4700_*.out     # node of each task
```

## 4. Copy back (PC, PowerShell)

```powershell
$R = 'bimstepf@romulus.ucmerced.edu:/home/bimstepf/endor/vault/NTM_ML/network_validation'
pscp -r "$R/outputs_t60_dt0.45" .
```

`pscp -r` into a folder that already exists overwrites files with the same
name and leaves the others. Copies can drop mid-way, so count again (step 2,
or `(Get-ChildItem outputs_t60_dt0.45 *.mat).Count`).

## 5. Look inside one file (PC, MATLAB)

```matlab
r = load('outputs_t60_dt0.45/rec_l2f__numeric__ECl__g0.0005_l0.02_d80_e20.mat');
fieldnames(r)              % N, M, F, W1, EdgeMass, trange, meta, cfg, env, ...
size(r.N)                  % [regions, time points]
r.trange(end) * 180        % last day (trange is in units of 180 days)
r.cfg.dt_max_days          % step limit it ran with ([] = default grid)
r.meta                     % the manifest row: method, gamma1, lambda1, ...
r.env                      % host, MATLAB version, threads, SLURM job, runtime context
```

`r.cfg` and `r.env` record the settings and machine that produced the file.
Check them before trusting a file.

## 6. Compare two runs (PC, MATLAB)

```matlab
compare_runs(file_a, file_b)          % a = reference, b = run to check
compare_runs(file_a, file_b, 60)      % only time points up to day 60
```

Runs on different time grids are compared at the time points they share.
For each of N, M, W1, F (and EdgeMass) it prints:

| column  | meaning |
|---------|---------|
| isequal | 1 = bit-identical |
| max abs | largest absolute difference |
| max rel | max abs / max abs(a) |
| NMSE    | sum (a-b)^2 / sum a^2 |
| first   | first shared day where they differ |

How to read the size of a difference:

- **isequal 1**: the same numbers. This is the only result that shows two
  codes are the same.
- **max rel ~1e-16 to 1e-15**: last-bit rounding. Different order of
  additions, nothing more.
- **NMSE ~1e-3 for numeric runs with different thread counts (or
  machines)**: the numeric sNTM amplifies last-bit differences over 180 days.
  Same code, so this is the floor for comparing numeric runs not made the
  same way.
- **Bigger than that**: a real difference (other settings, time step, code
  change). `first` shows when it starts.

To compare whole folders case by case:

```matlab
f = dir('outputs_t60_dt0.225/*.mat');
for k = 1:numel(f)
    compare_runs(fullfile('outputs_t60_dt0.45', f(k).name), fullfile(f(k).folder, f(k).name));
end
```

## 7. Common comparisons

- **Same code, two simulator folders** (original vs sped up). Same node, same
  thread count, so the answer should be `isequal 1`. Pin the node with
  `sbatch --nodelist=node01 ...`.
- **Repeats** (`outputs_t60_base`, `outputs_t60_repeats/rep2`, `rep3`): the
  same settings run three times. `isequal 1` means the cluster is
  reproducible. If not, compare `env.host` of the two files.
- **180-day run cut to 60 days**: `compare_runs(run180, run60, 60)`. Forward
  Euler steps from one `trange` point to the next, so the first 60 days of a
  180-day run on the same grid are the 60-day run.
- **Time step** (default grid vs 0.45 vs 0.225 days): forward Euler error is
  roughly proportional to the step. So the difference between the 0.45-day
  and 0.225-day runs estimates the error of the 0.225-day run. If that
  difference is small next to the differences between models, 0.225 is fine
  as the truth. The default grid already uses 0.45-day steps up to day 18
  (spacing 0.0025 × 180), then 0.9 days to day 54, then 1.8 days.
