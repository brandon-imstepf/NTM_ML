function pcshow(pc, varargin)
% Stand-in for the Computer Vision Toolbox pcshow, as far as Brainframe uses
% it: draw the points of pointCloud.m with scatter3 in their colours.
sz = 1;
for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'MarkerSize'), sz = varargin{k + 1}; end
end
if isempty(pc.Location), return, end
scatter3(pc.Location(:, 1), pc.Location(:, 2), pc.Location(:, 3), sz, pc.Color, '.');
end
