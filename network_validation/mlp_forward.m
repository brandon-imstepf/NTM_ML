function y = mlp_forward(net, X)
% Forward pass of an exported my_nn.NeuralNetwork (Linear/ReLU stack) on raw
% inputs X (n x 6). net comes from export_mlp_weights.py and holds:
%   W{k}, b{k}            layer weights (out x in) and biases (out x 1)
%   x_mean, x_std         input StandardScaler
%   y_mean, y_std         target StandardScaler
% Plain matrix algebra, so no Deep Learning Toolbox is needed on the cluster.

h = ((X - net.x_mean(:).') ./ net.x_std(:).').';
nl = numel(net.W);
for k = 1:nl
    h = net.W{k} * h + net.b{k}(:);
    if k < nl
        h = max(h, 0);
    end
end
y = h(:) * net.y_std + net.y_mean;
end
