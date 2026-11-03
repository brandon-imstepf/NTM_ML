"""Make nfc_probe.m: a copy of the ORIGINAL NetworkFluxCalculator.m with probe options.
  probe_mode  'joint'    (original: one fsolve over all active edges into a region)
              'per_edge' (one scalar fsolve per edge)
  nonneg_n    number of state components ode45 keeps non-negative (original: 1)
  only_target solve only the edges into this region (others stay 0)
  probe_A     candidate fluxes for the active edges into only_target: info.res is the
              shooting residual f_ss there, info.min_n the smallest n along the edge
The mass-edge computation is skipped. Everything else is the original code.
"""
import sys, re
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="utf-8").read()

def sub(old, new, count=1):
    global t
    assert t.count(old) >= 1, old
    t = t.replace(old, new, count)

sub("function [network_flux,mass_edge] = NetworkFluxCalculator(tau_x0,tau_xL,matdir,varargin)",
    "function [network_flux,mass_edge,info] = nfc_probe(tau_x0,tau_xL,matdir,varargin)\ninfo = struct();")
sub("parse(ip, varargin{:});",
    "addParameter(ip, 'probe_mode', 'joint');\naddParameter(ip, 'nonneg_n', 1);\n"
    "addParameter(ip, 'only_target', []);\naddParameter(ip, 'probe_A', []);\nparse(ip, varargin{:});")
sub("options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(n0));",
    "options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:ip.Results.nonneg_n);")
old_loop = """for i = 1:nroi
%     fprintf('ROI %d/%d\\n',i,nroi)"""
sub(old_loop, """targets = 1:nroi;
if ~isempty(ip.Results.only_target), targets = ip.Results.only_target; end
for i = targets""")
sub("""        network_flux(i_app,i) = fsolve(fun_ss,x0,options);
    end
end""", """        if strcmp(ip.Results.probe_mode, 'per_edge')
            rows = find(i_app);
            for e = 1:numel(rows)
                network_flux(rows(e),i) = fsolve(@(A) f_ss(A,tau_x0_i(e),tau_xL_i(e)),0,options);
            end
        else
            [xs, fv, ef, out] = fsolve(fun_ss,x0,options);
            network_flux(i_app,i) = xs;
            info.fval = fv; info.exitflag = ef; info.iterations = out.iterations;
        end
        info.rows = find(i_app); info.B = tau_x0_i; info.C = tau_xL_i;
        if ~isempty(ip.Results.probe_A)
            % k candidate values = the first k edges of the stack, solved as a k-edge system
            A = ip.Results.probe_A(:);  sel = 1:numel(A);
            B = tau_x0_i(sel);  C = tau_xL_i(sel);
            info.res = f_ss(A,B,C);
            info.min_n = min([n_ss_presyn(A,B,xmesh_presyn), n_ss_ais(A,B,xmesh_ais), ...
                n_ss_axon(A,B,xmesh_axon), n_ss_syncleft(A,B,xmesh_syncleft)], [], 2);
        end
    end
end
mass_edge = [];
return""")
open(dst, "w", encoding="utf-8").write(t)
print("wrote", dst)
