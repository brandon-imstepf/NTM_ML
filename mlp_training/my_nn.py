import torch
import torch.nn as nn
from collections import OrderedDict


class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, hidden_layers, output_size):
        super().__init__()

        layers = []
        layers.append(('layer_first', nn.Linear(input_size, hidden_size)))
        layers.append(('relu_1', nn.ReLU()))
        for i in range(hidden_layers-1):
            layers.append((f'layer_hidden_{i+1}', nn.Linear(hidden_size, hidden_size)))
            layers.append((f'relu_{i+2}', nn.ReLU()))
        layers.append(('layer_last', nn.Linear(hidden_size, output_size)))

        model = nn.Sequential(OrderedDict(layers))

        self.linear_relu_stack = model

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits
