% dataset audit. Suzanne asks "Which cases (and how many) parameter
% combinations can we recover." So this script looks for how many unique
% combinations there are, so we know the upper limit. 

% load datasets 
all_files = cell(9,1);   % only 1–9 exist
for i = 1:9
    fname = sprintf('data_bias_e5_%d.csv', i);
    filepath = fullfile('..','data/data_neural_network_model/training_data_e6/', fname);

    if exist(filepath, 'file')
        all_files{i} = readmatrix(filepath);
    else
        warning('File not found: %s', filepath);
        all_files{i} = [];
    end
end

% concatenate all data
all_data = vertcat(all_files{:});

% subset data based off seed (if needed)
% training_seed = 42;
% rng(training_seed)
% dataset_size = 5000;
% training_idx = randperm(size(all_data,1), dataset_size);
% training_dataset = all_data(training_idx,:);

% ensure training_dataset exists
training_dataset = all_data;   % remove this if you define it earlier

% find unique parameter combinations
[unique_rows, ~, group_idx] = unique(training_dataset(:,1:4), 'rows');
counts = accumarray(group_idx, 1);
