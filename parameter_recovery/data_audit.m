% dataset audit. Suzanne asks "Which cases (and how many) parameter
% combinations can we recover." So this script looks for how many unique
% combinations there are, so we know the upper limit. 


% load datasets 
filepaths = ["data_bias_e5_1.csv", 
    "data_bias_e5_2.csv",...
    "data_bias_e5_3.csv",...
    "data_bias_e5_4.csv",...
    "data_bias_e5_5.csv",...
    "data_bias_e5_6.csv",...
    "data_bias_e5_7.csv",...
    "data_bias_e5_8.csv",...
    "data_bias_e5_9.csv",...
    "data_bias_e5_1-.csv",];


all_files = cell(10,1);
for i = 1:10
    fname = sprintf('data_bias_e5_%d.csv', i);
    filepath = fullfile('..','data/data_neural_network_model/training_data_e6/', fname);
    all_files{i} = readmatrix(filepath);
end
all_data = vertcat(all_files{:});

% subset data based off seed
%training_seed = 42;
%rng(training_seed)
%dataset_size = 5000;
%training_idx = randperm(100000,dataset_size);
%training_dataset = training_superset(training_idx,:);

[unique_rows, ~, group_idx] = unique(training_dataset(:,1:4), 'rows');
counts = accumarray(group_idx, 1);