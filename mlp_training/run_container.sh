#! /usr/bin/env bash

# Usage: run_container.sh
# Arguments to nn-hyp-tuning.py:
#   1: target         - 'FValue' (phi_tau) or 'W1' (Delta_tau)
#   2: data_save_file - path for pickle output of metrics
#   3: random_state   - random seed
#   4: hidden_size    - nodes per hidden layer
#   5: hidden_layers  - number of hidden layers
#   6: epochs         - max training epochs
#   7: batch_size     - mini-batch size
#   8: fold_num       - held-out fold (1-10)
#   9: model_save     - path prefix for saved model
#  10: learn_rate     - SGD learning rate

docker run --mount type=bind,src=./output_data,dst=/app/output_data \
    nbarron00/nn-grid-search-e6 \
    'FValue' ./output_data/test_output 42 80 4 250 200 1 ./output_data/test_model 0.01
