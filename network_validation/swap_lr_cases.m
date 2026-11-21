function swap_lr_cases()
% After rerun_lr_cases.m: for every file outputs_lrfix/<folder>/<case>.mat,
% move the old <folder>/<case>.mat to outputs_lr_old/<folder>/ and the new
% file into <folder>/. Reversible: move the files back from outputs_lr_old/.
% Prints, for each case, the largest relative change in N.

here = fileparts(mfilename('fullpath'));
d = dir(fullfile(here, 'outputs_lrfix', '*', '*.mat'));
fprintf('%d rerun files\n', numel(d));
for k = 1:numel(d)
    [~, folder] = fileparts(d(k).folder);
    new = fullfile(d(k).folder, d(k).name);
    old = fullfile(here, folder, d(k).name);
    assert(isfile(old), 'no original for %s/%s', folder, d(k).name);
    a = load(old, 'N');  b = load(new, 'N');
    fprintf('%-20s %-50s max |dN|/max N = %.3g\n', folder, d(k).name, max(abs(a.N(:) - b.N(:))) / max(abs(a.N(:))));
    bak = fullfile(here, 'outputs_lr_old', folder);
    if ~isfolder(bak), mkdir(bak); end
    movefile(old, fullfile(bak, d(k).name));
    movefile(new, old);
end
end
