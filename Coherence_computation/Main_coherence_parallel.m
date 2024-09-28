% Main code to obtain weighted imaginary coherence

close all; clear all;

%load preprocessed EEG data
%load('C:\Aditi\Aditi_prepro_EEG\Prepro_EEG_Nov22\Prepro_EEG_cl_win_2s_0_2ch.mat');
load('D:\Postdoc\Only_on_hardisk\Depression_EEG_data\EEG_data_to_share\Prepro_EEG_Nov22_v2.mat')

%No_of_patients=length(data);

%Parallelize over 12 patients and hence run the code for these 12 patients
%at a time
No_of_patients=12;

trans=1000;

LEN_to_compute=500;         %2s of data
rep_times=8;                % Sum coherences over these many 2s segments
No_of_scales=78;            % No. of frequency scales
n_surr=30;
const=1.65;                 % Constant for checking if z-score for coherences is significant
N_channels=19;
ch_count=171;               % No. of channel pairs for which coherences are estimated
first_patient_no=1;

actual_patient_no=first_patient_no;

% Put the requisite data into a cell
for patient_no=1:No_of_patients
    if ~isempty(data(actual_patient_no).EEG_prepro.data)
        Patient_data_cell{patient_no}=data(actual_patient_no).EEG_prepro.data;
    else
        Patient_data_cell{patient_no}=[];
    end
    actual_patient_no=actual_patient_no+1;
end

clear data;

%Parallel processing for each subject
tic

spmd(No_of_patients)
    
    %  for patient_no=1:No_of_patients
    Patient_data=Patient_data_cell{labindex};
    Patient_data_cell={};               %As each worker now has its own dataset, clear the cell containing all
    [Coh_seg, z_sc_seg, sig_coh_seg]=coh_process_each_subject(Patient_data,trans,LEN_to_compute,rep_times,No_of_scales,n_surr,const,N_channels,ch_count);
    % Coh_pat{labindex}=Coh_seg_seg; mean_coh_s_pat=mean_coh_s_seg{labindex}; std_coh_s_pat{labindex}=std_coh_s_seg; z_sc_pat{labindex}=z_sc_seg; sig_coh_pat{labindex}=sig_coh_seg;
    %   save(['Patient_' num2str(labindex) '.mat'],'Coh_seg_seg','mean_coh_s_seg','std_coh_s_seg','z_sc_seg','sig_coh_seg')
    %  end
    
end

toc


% save the results
actual_patient_no=first_patient_no;

for patient_no=1:No_of_patients
    Coh=Coh_seg{patient_no};
    z_sc=z_sc_seg{patient_no};
    sig_coh=sig_coh_seg{patient_no};
   % save(['C:\Aditi\Aditi_prepro_EEG\Coherence_analysis_res_Nov22\res_patient' num2str(actual_patient_no) '.mat'],'Coh','z_sc','sig_coh','const','n_surr');
    save(['D:\Postdoc\EEG_coherence_expts\Codes_to_share\For_coherence_computation\Coherence_analysis_res_test\res_patient' num2str(actual_patient_no) '.mat'],'Coh','z_sc','sig_coh','const','n_surr');
    actual_patient_no=actual_patient_no+1;
    clear Coh z_sc sig_coh
end

