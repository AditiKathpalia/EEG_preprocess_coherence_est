close all; clear all;

% load the EEG dataset
load('D:\Postdoc\Only_on_hardisk\Depression_EEG_data\Barbora_data\concatenated\A01_loaded_concatenated_big');

%load EEG cap and label info
load('cap_and_label_info.mat');

data = preprocessing_all_subj(data,cap,labels);         % the variable 'data' contains preprocessed data in EEGLAB format



