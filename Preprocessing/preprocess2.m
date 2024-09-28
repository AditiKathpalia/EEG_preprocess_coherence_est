%To preprocess EEG of one subject

function [EEG_clean] = preprocess2(EEG,cap,labels)  

if ~isfield(EEG,'event')
    EEG.event = [];
end
EEG.etc = [];
EEG.icachansind = [];

EEG.event = struct('latency',[],'duration',[],'channel',[],'type',[],'code',[]);
EEG = eeg_checkset(EEG);


% Remove any channels in excess of the required 19 electrodes, common
% across all subjects
%chan_to_retain=[1:21];
chans_to_remove=remove_channels(EEG,labels);
EEG = pop_select(EEG, 'nochannel', chans_to_remove);    %EEG = pop_select(EEG, 'channel', chan_to_retain);    %?? names of these channels dont have location
EEG = eeg_checkset(EEG);

% Set location of the remaining channels
EEG.chanlocs = cap;
EEG = eeg_checkset(EEG);

% Downsample to 250Hz, if sampling rate is 1000Hz
if EEG.srate==1000
    data1=downsample(EEG.data',4);
    EEG=rmfield(EEG,'data');
    EEG.data=data1';
    EEG.srate=250;
    EEG = eeg_checkset( EEG );
end

% Remove first 30s of data and anything after 10 mins as that might be
% affected by luminiscence given to the subjects
% pts_to_remove=30*EEG.srate;
% EEG = pop_select( EEG,'time',[1 pts_to_remove; EEG.pnts-pts_to_remove+1 EEG.pnts] );
EEG = pop_select( EEG,'time',[30 570] );
EEG_clean = eeg_checkset( EEG );

%set average reference 
EEG = pop_reref( EEG, []);
EEG = eeg_checkset( EEG );

%Filter data
EEG = pop_eegfiltnew(EEG, 1, 40);       %Check what kind of filtering is happening - maybe change filter type for causality analysis
EEG = eeg_checkset( EEG );

cleaning=1;
if cleaning==1 
%Remove bad channels, peaks and windows - compare methods
%ASR - toolbox - interpolates bad signal
%Data.rawevents{i}=EEG.event;
EEG= clean_windows(EEG, 0.2,[],2);       % Check how many go missing, make window length as 2s
EEG_clean = eeg_checkset( EEG );

end



end

