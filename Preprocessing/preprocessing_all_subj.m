function [data] = preprocessing_all_subj(data,cap,labels)

%Pass data of all subjects here in order to preprocess as follows

subjecti=1;

%% Pass one subgect at a time
while subjecti < length(data)+1

%for subjecti = 1:length(data)
    
    fprintf('\n Subject no: %3.0f \n', subjecti);
    if ~isempty(data(subjecti).EEG) && rem(subjecti,2)==1 && subjecti < 289     % As odd numbered subjects (visit 1) below 288 have longer recordings
        data(subjecti).EEG_prepro = preprocess2(data(subjecti).EEG,cap,labels);
    elseif ~isempty(data(subjecti).EEG) && rem(subjecti,2)==0 && subjecti > 289 % As even numbered subjects (visit 2) above 289 have longer recordings
        data(subjecti).EEG_prepro = preprocess2(data(subjecti).EEG,cap,labels);
    elseif ~isempty(data(subjecti).EEG)
        data(subjecti).EEG_prepro = preprocess(data(subjecti).EEG,cap,labels);
    end
    
    subjecti=subjecti+1;

end
data = rmfield(data,'EEG');

end