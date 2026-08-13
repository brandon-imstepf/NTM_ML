import pandas as pd
import numpy as np
import copy

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

import pickle
import sys

from my_nn import NeuralNetwork
from nn_hyp_utils import model_fit

# ---------------------------------------------------------------------------
# Parse command line arguments
# ---------------------------------------------------------------------------

target = sys.argv[1]
target = [target]

data_save_file = sys.argv[2]
random_state = int(sys.argv[3])
hidden_size_in = int(sys.argv[4])
hidden_layers_in = int(sys.argv[5])
epochs_in = int(sys.argv[6])
batch_size = int(sys.argv[7])
fold_num = int(sys.argv[8])
model_save_file = sys.argv[9]
learn_rate_in = float(sys.argv[10])

# ---------------------------------------------------------------------------
# Early stopping parameters (per paper: min_delta=1e-7, patience=20)
# ---------------------------------------------------------------------------

early_stop_patience = 20
early_stop_min_delta = 1e-7

# ---------------------------------------------------------------------------
# Select device
# ---------------------------------------------------------------------------

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('Using GPU')
else:
    device = torch.device('cpu')
    print('GPU not available, using CPU instead')

# ---------------------------------------------------------------------------
# Dataset indices: 10-fold CV leaves out fold_num for validation
# ---------------------------------------------------------------------------

train_datasets = np.arange(10) + 1
train_datasets = np.delete(train_datasets, np.where(train_datasets == fold_num))

features = ['gamma1', 'lambda1', 'delta', 'epsilon', 'NRow', 'NCol']

# ---------------------------------------------------------------------------
# Fit Z-score scaler on fold 1 training data only (prevents data leakage)
# ---------------------------------------------------------------------------

fit_data_filepath = './training_data/data_bias_e5_1.csv'
scale_data = pd.read_csv(fit_data_filepath)

X_scale = scale_data[features].values
y_scale = scale_data[target].values

PredictorScaler = StandardScaler()
TargetScaler = StandardScaler()

PredScaleFit = PredictorScaler.fit(X_scale)
TargetScaleFit = TargetScaler.fit(y_scale)

# ---------------------------------------------------------------------------
# Define model architecture
# ---------------------------------------------------------------------------

input_size = 6
output_size = 1
hidden_size = hidden_size_in
hidden_layers = hidden_layers_in

epoch_chunk_size = 50
epochs_n = epochs_in

size = (int(epochs_n / epoch_chunk_size), 2, 6)
output_data = np.zeros(size)

model = NeuralNetwork(input_size, hidden_size, hidden_layers, output_size)
model = model.to(device)

criterion = nn.MSELoss()
optimizer = optim.SGD(model.parameters(), lr=learn_rate_in)

# ---------------------------------------------------------------------------
# Prepare test dataset (held-out fold)
# ---------------------------------------------------------------------------

filepath_test = './training_data/data_bias_e5_' + str(fold_num) + '.csv'
test_data = pd.read_csv(filepath_test)

X_test = test_data[features].values
X_test_scaled = PredScaleFit.transform(X_test)
X_test_scaled_ten = torch.tensor(X_test_scaled, dtype=torch.float32).to(device)
y_test = test_data[target].values

# ---------------------------------------------------------------------------
# Prepare validation dataset (dataset 11, separate from all folds)
# ---------------------------------------------------------------------------

filepath_val = './training_data/data_bias_e5_11.csv'
val_data = pd.read_csv(filepath_val)

X_val = val_data[features].values
X_val_scaled = PredScaleFit.transform(X_val)
X_val_scaled_ten = torch.tensor(X_val_scaled, dtype=torch.float32).to(device)
y_val = val_data[target].values

# ---------------------------------------------------------------------------
# Training loop with early stopping
# ---------------------------------------------------------------------------

data_index = 0
min_metric = float('inf')
best_model = copy.deepcopy(model)

# Early stopping state
best_nmse = float('inf')
epochs_without_improvement = 0

for epoch_i in range(epochs_n):

    model.train()

    for dataset_index in train_datasets:

        filepath = './training_data/data_bias_e5_' + str(dataset_index) + '.csv'
        train_data = pd.read_csv(filepath)

        X_train = train_data[features].values
        y_train = train_data[target].values

        X_train = PredScaleFit.transform(X_train)
        y_train = TargetScaleFit.transform(y_train)

        X_train = torch.tensor(X_train, dtype=torch.float32).to(device)
        y_train = torch.tensor(y_train, dtype=torch.float32).to(device)

        dataset = TensorDataset(X_train, y_train)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        for batch_i, (inputs, targets) in enumerate(dataloader):

            outputs = model(inputs)
            loss = criterion(outputs, targets)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    # Print progress
    print(f'Epoch [{epoch_i+1}/{epochs_n}], Loss: {loss.item():.4f}')

    # ----- Evaluate and check early stopping every epoch_chunk_size epochs -----

    if (epoch_i + 1) % epoch_chunk_size == 0:

        print('Saving Data...')

        model.eval()
        with torch.no_grad():

            # Predict test data and get error metrics
            preds_test_scaled = model(X_test_scaled_ten)
            preds_test = TargetScaleFit.inverse_transform(
                preds_test_scaled.cpu().numpy()
            )
            metrics_test = model_fit(
                y_test.flatten(), preds_test.flatten(), verbose=False
            )

            # Predict validation data and get error metrics
            preds_val_scaled = model(X_val_scaled_ten)
            preds_val = TargetScaleFit.inverse_transform(
                preds_val_scaled.cpu().numpy()
            )
            metrics_val = model_fit(
                y_val.flatten(), preds_val.flatten(), verbose=False
            )

        output_data[data_index, 0, :] = metrics_test
        output_data[data_index, 1, :] = metrics_val
        data_index += 1

        with open(data_save_file, 'wb') as f:
            pickle.dump(output_data, f)

        # Track best model by test NMSE (index 4 in model_fit output)
        current_nmse = metrics_test[4]

        if current_nmse < min_metric:
            best_model = copy.deepcopy(model)
            min_metric = current_nmse

            dummy_input = torch.randn((1, 6), dtype=torch.float32).to(device)
            traced_model = torch.jit.trace(best_model.forward, dummy_input)
            traced_model.save(model_save_file + ".pt")

        # ----- Early stopping check -----
        if best_nmse - current_nmse > early_stop_min_delta:
            best_nmse = current_nmse
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += epoch_chunk_size

        if epochs_without_improvement >= early_stop_patience * epoch_chunk_size:
            print(
                f'Early stopping at epoch {epoch_i+1}: '
                f'no NMSE improvement > {early_stop_min_delta} '
                f'for {early_stop_patience} evaluation windows.'
            )
            break
