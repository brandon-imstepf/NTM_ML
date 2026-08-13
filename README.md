# NTM Surrogate Models

Machine learning surrogate models for accelerating the Network Transport Model (NTM) of tau protein spread.

## Project Structure

### `mlp_training/`
Docker-containerized MLP neural network training with hyperparameter grid search.

- `my_nn.py` — MLP architecture (input → ReLU hidden layers → output)
- `nn_hyp_utils.py` — Error metric utility function (R², Pearson R, MSE, NMSE)
- `nn-hyp-tuning.py` — Main training script with 10-fold CV, early stopping (patience=20, min_delta=1e-7), and model saving
- `shapley_analysis.ipynb` — SHAP feature importance analysis on trained MLP models
- `Dockerfile`, `run_script.sh`, `run_container.sh` — Docker build and run files
- `requirements.txt` — Python dependencies
- `training_data/` — Single-edge NTM training CSVs (`data_bias_e5_{1-11}.csv`; files 1-10 are CV folds, file 11 is held-out validation)
- `output_data/` — Saved models and training metrics (populated at runtime)

### `linear_regression/`
Classical linear regression model training, forward feature selection, and feature importance analysis.

- `classical_regressions.ipynb` — Trains four regression variants (Linear, Quadratic, Joint, All) with 10-fold CV and bootstrapped test performance
- `feature_selection_and_importance.ipynb` — Forward stepwise feature selection on the 'All' model and OLS t-tests for feature importance
- `training_data/` — Place `bias_train_e3.csv` and `bias_val_e3.csv` here
- `outputs/` — Saved results (populated at runtime)

### `visualizations/`
Notebooks for reproducing paper figures. Each assumes input data is in `./data/`.

- `figure4_regression_graph.ipynb` — Figure 4: Linear regression forward feature selection NMSE curves
- `figure6_nn_heatmaps.ipynb` — Figure 6: MLP hyperparameter grid search heatmaps + training loss
- `figure8b_figure11_network_sim.ipynb` — Figures 8b & 11: Network-wide NTM tau trajectories and R² scatter plots
- `figure9_bio_figure.ipynb` — Figure 9: Anterograde vs. retrograde transport regime comparison
- `figure10_shapley.ipynb` — Figure 10: MLP Shapley analysis beeswarm plots and mean |SHAP| values
- `data/` — Place simulation `.mat` files, trained model `.pt` files, and training CSVs here

## Training Data Format

All CSVs share the same column format:
```
gamma1, lambda1, delta, epsilon, NRow, NCol, FValue, W1
```
- `gamma1`–`NCol`: 6 input features (biophysical NTM parameters + nodal tau concentrations)
- `FValue`: tau flux φ_τ (prediction target)
- `W1`: edge mass correction Δ_τ (prediction target)

## Usage

**MLP Training (Docker):**
```bash
cd mlp_training
docker build -t ntm-mlp .
./run_container.sh
```

**Linear Regression:** Run `classical_regressions.ipynb` then `feature_selection_and_importance.ipynb` in Jupyter.

**Visualizations:** Run the respective notebook after placing required data files in `visualizations/data/`.
