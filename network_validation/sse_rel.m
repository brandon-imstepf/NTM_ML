function v = sse_rel(s, n)
% Sum of squared errors relative to the sum of squared reference values,
%   v = sum((s - n).^2) / sum(n.^2),
% the network error used before 2026-11-10 (then called "NMSE" in the code).
% Kept as a secondary column; the reported metric is nmse_eq7.m. For a fixed
% reference n the two differ by a constant factor, sum((n - mean n)^2) / sum(n^2).
s = s(:);  n = n(:);
v = sum((s - n).^2) / sum(n.^2);
end
