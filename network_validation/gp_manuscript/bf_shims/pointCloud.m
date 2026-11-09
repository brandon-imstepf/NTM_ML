function pc = pointCloud(xyz, varargin)
% Stand-in for the Computer Vision Toolbox pointCloud, as far as Brainframe
% uses it: the point locations and their colours (Intensity is ignored).
% Only on the path when that toolbox is missing (glass_brain_images.m).
pc.Location = xyz;
pc.Color = [];
for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'Color'), pc.Color = varargin{k + 1}; end
end
end
