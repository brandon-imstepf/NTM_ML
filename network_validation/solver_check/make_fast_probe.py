"""Make nfc_fast_probe.m from the EDITED NetworkFluxCalculator.m (same flux as the
original, verified bit-identical): records fsolve's exit flag and the residual of
every edge, can solve each edge alone ('per_edge', 1), and skips the mass edge.
  [F, ~, info] = nfc_fast_probe(...)  info.flag, info.res: nroi x nroi (NaN off the solved edges)
"""
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="utf-8").read()

def sub(old, new):
    global t
    assert t.count(old) == 1, old
    t = t.replace(old, new)

sub("function [network_flux,mass_edge] = NetworkFluxCalculator(tau_x0,tau_xL,matdir,varargin)",
    "function [network_flux,mass_edge,info] = nfc_fast_probe(tau_x0,tau_xL,matdir,varargin)\ninfo = struct();")
sub("addParameter(ip, 'flux_parallel', 0)",
    "addParameter(ip, 'flux_parallel', 0)\naddParameter(ip, 'per_edge', 0)")
sub("jobs_A = solve_jobs(jobs_B,jobs_C,jobs_x0,cst,jac_diag,fsolve_options,fsolve_options_jac,ip.Results.flux_parallel);\n"
    "for r = 1:nj\n    network_flux(jobs_app{r},jobs_i(r)) = jobs_A{r};\nend",
    "[jobs_A,jobs_ef,jobs_fv] = solve_jobs(jobs_B,jobs_C,jobs_x0,cst,jac_diag,fsolve_options,fsolve_options_jac,ip.Results.per_edge);\n"
    "info.flag = nan(nroi); info.res = nan(nroi);\n"
    "for r = 1:nj\n    network_flux(jobs_app{r},jobs_i(r)) = jobs_A{r};\n"
    "    info.flag(jobs_app{r},jobs_i(r)) = jobs_ef{r};\n    info.res(jobs_app{r},jobs_i(r)) = jobs_fv{r};\nend\n"
    "mass_edge = [];\nreturn")
sub("""function jobs_A = solve_jobs(jobs_B,jobs_C,jobs_x0,c,jac_diag,opts,opts_jac,par)
% (parfor cannot sit in a function that has nested functions)
nj = numel(jobs_B);
jobs_A = cell(1,nj);
if par
    parfor r = 1:nj
        jobs_A{r} = solve_region(jobs_B{r},jobs_C{r},jobs_x0{r},c,jac_diag,opts,opts_jac);
    end
else
    for r = 1:nj
        jobs_A{r} = solve_region(jobs_B{r},jobs_C{r},jobs_x0{r},c,jac_diag,opts,opts_jac);
    end
end
end""", """function [jobs_A,jobs_ef,jobs_fv] = solve_jobs(jobs_B,jobs_C,jobs_x0,c,jac_diag,opts,opts_jac,per_edge)
% probe: joint solve per region (original) or one scalar solve per edge
nj = numel(jobs_B);
[jobs_A,jobs_ef,jobs_fv] = deal(cell(1,nj));
for r = 1:nj
    if per_edge
        ne = numel(jobs_B{r});
        [A,ef,fv] = deal(zeros(ne,1));
        for e = 1:ne
            [A(e),fv(e),ef(e)] = fsolve(@(a) f_ss_c(a,jobs_B{r}(e),jobs_C{r}(e),c),jobs_x0{r}(e),opts);
        end
        jobs_A{r} = A; jobs_ef{r} = ef; jobs_fv{r} = fv;
    else
        [jobs_A{r},jobs_fv{r},ef] = fsolve(@(A) f_ss_c(A,jobs_B{r},jobs_C{r},c),jobs_x0{r},opts);
        jobs_ef{r} = repmat(ef,numel(jobs_B{r}),1);
    end
end
end""")
open(dst, "w", encoding="utf-8").write(t)
print("wrote", dst)
