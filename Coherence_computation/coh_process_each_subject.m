function [Coh_seg, z_sc_seg, sig_coh_seg]=coh_process_each_subject(Patient_data,trans,LEN_to_compute,rep_times,No_of_scales,n_surr,const,N_channels,ch_count)

% Function to process weighted imaginary coherences for a single subject
% Returns coherence for chosen different set of scales between each pair of
% channels considered separately for different time segments.
% For more details, please refer:
% Paluš, M., Kathpalia, A., & Brunovský, M. (2023, May). EEG connectivity
% in treatment of major depressive disorder: Tackling the conductivity
% effects. In 2023 14th International Conference on Measurement (pp.
% 76-79). IEEE.
% For description of inputs to this function, please see 'Main_coherence_parallel'
% Outputs:
% Coh_seg is the matrix of computed weighted imaginary coherences and has the
% dimension of the No. of channel pairs*No. of scales*No. of time segments
% in data
% z_sc is the matrix consisting of z scores for each scale computed based
% on the surrogates. It is of the same dimension as 'Coh_seg'
% sig_coh_seg is the matrix consisting of significant weighted imaginary
% coherences. It is the same dimension as 'Coh_seg'. In this the values in
% Coh_seg, which are non significant based on 'z_sc' are set to zero and
% the rest remain the same as in 'Coh_seg'

%clear data
seg_no=0;
i=trans+1;

if ~isempty(Patient_data)
    
    while true                                          % Till the data lasts
        remaining_len=size(Patient_data,2)-i-trans;     % Also, remove trans elements from the end
        if remaining_len > LEN_to_compute*rep_times
            X=Patient_data(:,i:i+LEN_to_compute*rep_times-1);       %Data for one segment = 2s*8 times
            seg_no=seg_no+1;
            %  del_coh_reptimes(:,:,seg_no)=zeros(19,19);
            for count=1:rep_times
                for channel_no=1:N_channels
                    X_s=X(channel_no, (count-1)*LEN_to_compute+1:count*LEN_to_compute); % Data for 1 count/ 2s
                    [s_a]=FourierTransform(X_s',n_surr);         % Generate surrogates
                    [WAVE(:,:,channel_no,count)]=wavelet_comp(X_s, No_of_scales);       %compute wavelet coefficients from signal
                    for sur_no=1:n_surr
                        [WAVE_s(:,:,channel_no,count,sur_no)]=wavelet_comp(s_a(:,sur_no)', No_of_scales);       %compute wavelet coefficients from surrogates
                    end
                    clear X_s s_a
                    % X_s=[]; s_a=[];
                end
            end
            Coh_ch=zeros(ch_count,No_of_scales);
            %                 mean_coh_s_ch=zeros(ch_count,No_of_scales);
            %                 std_coh_s_ch=zeros(ch_count,No_of_scales);
            z_sc_ch=zeros(ch_count,No_of_scales);
            sig_coh_ch=zeros(ch_count,No_of_scales);
            ch_loop_no=0;
            for channel_no1= 1:N_channels
                for channel_no2= channel_no1+1:N_channels
                    ch_loop_no=ch_loop_no+1;
                    Coh_scale=zeros(No_of_scales,1);
                    mean_coh_s=zeros(No_of_scales,1);
                    std_coh_s=zeros(No_of_scales,1);
                    z_sc=zeros(No_of_scales,1);
                    sig_coh=zeros(No_of_scales,1);
                    
                    for scale_level= 1:No_of_scales
                        %  numb=numb+1;
                        w1_tot=0;
                        w2_tot=0;
                        w3_tot=0;
                        w1s_tot=zeros(n_surr,1);
                        w2s_tot=zeros(n_surr,1);
                        w3s_tot=zeros(n_surr,1);
                        Coh_s_seg=zeros(n_surr,1);
                        for count=1:rep_times
                            w1_s=zeros(n_surr,1);
                            w2_s=zeros(n_surr,1);
                            w3_s=zeros(n_surr,1);
                            WAVE1=process_wave(WAVE(scale_level,:,channel_no1,count));
                            WAVE2=process_wave(WAVE(scale_level,:,channel_no2,count));
                            [w1,w2,w3]=coherence_comp(WAVE1,WAVE2);
                            w1_tot=w1_tot+w1;
                            w2_tot=w2_tot+w2;
                            w3_tot=w3_tot+w3;
                            clear w1 w2 w3 WAVE1 WAVE2
                            % w1=[]; w2=[]; w3=[]; WAVE1=[]; WAVE2=[];
                            % Total_coh(channel_no1,channel_no2)=Total_coh(channel_no1,channel_no2)+Coh;
                            
                            %     sur process and coh comp!
                            for sur_no=1:n_surr
                                WAVE1_s=process_wave(WAVE_s(scale_level,:,channel_no1,count,sur_no));
                                WAVE2_s=process_wave(WAVE_s(scale_level,:,channel_no2,count,sur_no));
                                [w1_s(sur_no),w2_s(sur_no),w3_s(sur_no)]=coherence_comp(WAVE1_s,WAVE2_s);
                                clear WAVE1_s WAVE2_s
                                % WAVE1_s=[]; WAVE2_s=[];
                            end
                            w1s_tot=w1s_tot+w1_s;
                            w2s_tot=w2s_tot+w2_s;
                            w3s_tot=w3s_tot+w3_s;
                            clear w1_s w2_s w3_s
                            %    w1_s=[]; w2_s=[]; w3_s=[];
                        end
                        Coh_scale(scale_level)=weight_imag_coh(w1_tot,w2_tot,w3_tot);
                        for sur_no=1:n_surr
                            Coh_s_seg(sur_no)=weight_imag_coh(w1s_tot(sur_no),w2s_tot(sur_no),w3s_tot(sur_no));
                        end
                        
                        mean_coh_s(scale_level)=mean(Coh_s_seg);
                        std_coh_s(scale_level)=std(Coh_s_seg);
                        % qty_comp=mean_coh_s(scale_level)+std_coh_s(scale_level)*const;
                        z_sc(scale_level)=(Coh_scale(scale_level)-mean_coh_s(scale_level))/(std_coh_s(scale_level));
                        %   if Coh_seg(scale_level)> qty_comp
                        if z_sc(scale_level)> const
                            sig_coh(scale_level)=Coh_scale(scale_level);
                        end
                        % clear WAVE1 WAVE2
                        clear Coh_s_seg
                        % qty_comp=[]; Coh_s_seg=[];
                    end
                    Coh_ch(ch_loop_no,:)=Coh_scale;
                    %   mean_coh_s_ch(ch_loop_no,:)=mean_coh_s;
                    %  std_coh_s_ch(ch_loop_no,:)=std_coh_s;
                    z_sc_ch(ch_loop_no,:)=z_sc;
                    sig_coh_ch(ch_loop_no,:)=sig_coh;
                    
                    % del_coh(channel_no1,channel_no2)=del_coh(channel_no1,channel_no2)/numb;
                    clear Coh_scale mean_coh_s std_coh_s z_sc sig_coh
                    % Coh_scale=[]; mean_coh_s=[]; std_coh_s=[]; z_sc=[]; sig_coh=[];
                end
            end
            Coh_seg(:,:,seg_no)=Coh_ch;
            %   mean_coh_s_seg(:,:,seg_no)=mean_coh_s_ch;
            %  std_coh_s_seg(:,:,seg_no)=std_coh_s_ch;
            z_sc_seg(:,:,seg_no)=z_sc_ch;
            sig_coh_seg(:,:,seg_no)=sig_coh_ch;
            clear Coh_seg_ch sig_coh_ch z_sc_ch WAVE WAVE_s
            %  Coh_seg_ch=[]; mean_coh_s_ch=[]; std_coh_s_ch=[]; z_sc_ch=[]; WAVE=[]; WAVE_s=[];
            %   del_coh_reptimes(:,:,seg_no)=del_coh_reptimes(:,:,seg_no)+del_coh;
            
        else
            break
        end
        i=i+LEN_to_compute*rep_times;
    end
    
else
    Coh_seg=[];
    z_sc_seg=[];
    sig_coh_seg=[];
end
%  clear Patient_data
