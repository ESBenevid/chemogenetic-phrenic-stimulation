clear all
close all
clc

%% load data
root = '~\Fig_1_2_3_WT_ChAT_Cre_mice\data'; % change to directory where data files are located
cd(root);
strain = "Chat-Cre"; %choose to analyze wild-type or ChAT-Cre mouse data; type = "WT" for wild-type or "Chat-Cre" for ChAT-Cre data
save_data = 'no';
animals = dir;

if strcmp(strain, "Chat-Cre")
    an_o_interest = {'11_20BF','12_20BF','13_20BF','14_20BF','15_20BF','41_20BF','42_20BF','43_20BF','45_20BF'};
    animals = animals(ismember({animals.name},an_o_interest));
    disp("Analyzing Chat-Cre data")
elseif strcmp(strain, "WT")
    an_o_interest = {'2_20BF','20_20BF','21_20BF','22_20BF','23_20BF','25_20BF','26_20BF','35_20BF',...
        '36_20BF','37_20BF','39_20BF'};
    animals = animals(ismember({animals.name},an_o_interest));
    disp("Analyzing wildtype data")
else
end


for a = 1:numel(animals)
    cd([root '\' animals(a).name '\mat_export'])
    files = dir('*.mat');
    
    for f = 1:numel(files)
        load(files(f).name)
        
        if strcmp(files(f).name,'21_20BF_20211011_full_trace.mat') || strcmp(files(f).name,'35_20BF_20211201_full_trace.mat')...
                || strcmp(files(f).name,'22_20BF_20211011_full_trace.mat') || strcmp(files(f).name,'36_20BF_20211202_full_trace.mat')...
                || strcmp(files(f).name,'37_20BF_20211203_full_trace.mat') || strcmp(files(f).name,'13_20BF_20210607_full_trace.mat')...
                || strcmp(files(f).name,'12_20BF_20210609_full_trace.mat') || strcmp(files(f).name,'12_20BF_20210609_full_trace.mat')
            continue
        else
        end
        
        clear pre_CTA_tongue pre_CTA_ldia pre_CTA_rdia AUC_tongue_alt AUC_ldia_alt AUC_rdia_alt
        %% define sampling rates
        if strcmp(animals(a).name, '2_20BF')
            SR = 5000;
            ms = SR/1000;
            
            %% filter EMG traces
            [B,A] = butter(2, [100 1000]/(SR/2));
            
            L_Dia.values = resample(L_Dia.values,SR,1/L_Dia.interval);
            R_Dia.values = resample(R_Dia.values,SR,1/R_Dia.interval);
            
            filt_ldia = filtfilt(B,A,L_Dia.values);
            filt_rdia = filtfilt(B,A,R_Dia.values);
            filt_tongue = filtfilt(B,A,L_Dia.values);
            
        else
            tongue_fs = 1/Tongue.interval;
            SR = 1/L_Dia.interval;
            SR = 1/R_Dia.interval;
            ms = SR/1000;
            
            %% filter EMG traces
            [B,A] = butter(2, [100 1000]/(tongue_fs/2));
            filt_tongue = filtfilt(B,A,Tongue.values);
            filt_ldia = filtfilt(B,A,L_Dia.values);
            filt_rdia = filtfilt(B,A,R_Dia.values);
        end
        %% rectify and smooth EMG traces
        rec_tongue = abs(filt_tongue);
        rec_ldia = abs(filt_ldia);
        rec_rdia = abs(filt_rdia);
        
        med_tongue = movmedian(rec_tongue,50*ms);
        med_ldia = movmedian(rec_ldia,50*ms);
        med_rdia = movmedian(rec_rdia,50*ms);
        
        smooth_tongue = smooth(med_tongue,50*ms);
        smooth_ldia = smooth(med_ldia,50*ms);
        smooth_rdia = smooth(med_rdia,50*ms);
        
        %% identify trace portions of interest
        if f == 1
            times_o_interest = [420 600];
            data_row = 1;
        elseif f==2
            times_o_interest = [420 600];
            data_row = 2;
        elseif f == 3
            %times_o_interest = [1 600; 601 1200; 1201 1800];
            times_o_interest = [120 300; 720 900; 1620 1800; 3420 3600;...
                5220 5400];
            data_row = [3 4 5 6 7];
        else
        end
        
        %% loop through times of interest
        clear rising_epochs_tongue falling_epochs_tongue
        
        for t = 1:size(times_o_interest,1)
            if (times_o_interest(t,2)*tongue_fs) > length(smooth_tongue)
                continue
            else
                smooth_tongue_trim = smooth_tongue((times_o_interest(t,1)*tongue_fs):times_o_interest(t,2)*tongue_fs);
                smooth_ldia_trim = smooth_ldia((times_o_interest(t,1)*SR):times_o_interest(t,2)*SR);
                smooth_rdia_trim = smooth_rdia((times_o_interest(t,1)*SR):times_o_interest(t,2)*SR);
            end
            
            %% define rising and falling epochs for tongue EMG
            if strcmp(animals(a).name, '2_20BF')
                threshold = 0.02;
            elseif strcmp(animals(a).name, '13_20BF')
                threshold = 0.0015;
            else
                threshold = 0.005;
            end
            above_thresh = find(smooth_tongue_trim > threshold);
            spaces = diff(above_thresh);
            gaps = find(spaces > 5);
            rising_epochs_tongue = above_thresh([1; gaps+1]);
            falling_epochs_tongue = above_thresh([gaps; end]);

                        
            %% calculation pre-CTAs
            clear pre_CTA_tongue pre_CTA_ldia pre_CTA_rdia LDia_tonic RDia_tonic Tongue_tonic pre_CTA_ldia_sub pre_CTA_rdia_sub...
                pre_CTA_tongue_sub RDia_AUC LDia_AUC Tongue_AUC RDia_pk_amp LDia_pk_amp Tongue_pk_amp
            for r = 1:numel(rising_epochs_tongue)
                start = rising_epochs_tongue(r) - 250*ms;
                stop = rising_epochs_tongue(r) + 450*ms;
                
                if (start < 0) || (stop > numel(smooth_tongue_trim))
                    continue
                else
                end
                
                
                pre_CTA_tongue(r,:) = smooth_tongue_trim(start:stop);
                pre_CTA_ldia(r,:) = smooth_ldia_trim(start:stop);
                pre_CTA_rdia(r,:) = smooth_rdia_trim(start:stop);
            end
            
            
            RDia_tonic = mean(pre_CTA_rdia(:,1:20),2);
            LDia_tonic = mean(pre_CTA_ldia(:,1:20),2);
            Tongue_tonic = mean(pre_CTA_tongue(:,1:20),2);
            
            pre_CTA_ldia_sub = pre_CTA_ldia - LDia_tonic;
            pre_CTA_rdia_sub = pre_CTA_rdia - RDia_tonic;
            pre_CTA_tongue_sub = pre_CTA_tongue - Tongue_tonic;
            
            RDia_AUC = trapz(pre_CTA_rdia_sub,2);
            LDia_AUC = trapz(pre_CTA_ldia_sub,2);
            Tongue_AUC = trapz(pre_CTA_tongue_sub,2);
            
            RDia_pk_amp = max((pre_CTA_rdia - min(pre_CTA_rdia,[],2)),[],2);
            LDia_pk_amp = max((pre_CTA_ldia - min(pre_CTA_ldia,[],2)),[],2);
            Tongue_pk_amp = max((pre_CTA_tongue - min(pre_CTA_tongue,[],2)),[],2);

            
            %% Tongue CTA/AUC/PK2PK
            CTA_tongue(data_row(t),:,a) = mean(pre_CTA_tongue);
            avg_Tongue_tonic(data_row(t),a) = mean(Tongue_tonic);
            AUC_tongue(data_row(t),a) = mean(Tongue_AUC);
            PK2PK_tongue(data_row(t),a) = mean(Tongue_pk_amp);
            
            CTA_tongue(CTA_tongue == 0) = nan;
            avg_Tongue_tonic(avg_Tongue_tonic == 0) = nan;
            AUC_tongue(AUC_tongue == 0) = nan;
            PK2PK_tongue(PK2PK_tongue == 0) = nan;
            
            %% LDia CTA/AUC/PK2PK
            CTA_ldia(data_row(t),:,a) = mean(pre_CTA_ldia);
            avg_LDia_tonic(data_row(t),a) = mean(LDia_tonic);
            AUC_ldia(data_row(t),a) = mean(LDia_AUC);
            PK2PK_ldia(data_row(t),a) = mean(LDia_pk_amp);
            
            CTA_ldia(CTA_ldia == 0) = nan;
            avg_LDia_tonic(avg_LDia_tonic == 0) = nan;
            AUC_ldia(AUC_ldia == 0) = nan;
            PK2PK_ldia(PK2PK_ldia == 0) = nan;
            
            %% RDia CTA/AUC/PK2PK/resp rate
            CTA_rdia(data_row(t),:,a) = mean(pre_CTA_rdia);
            avg_RDia_tonic(data_row(t),a) = mean(RDia_tonic);
            AUC_rdia(data_row(t),a) = mean(RDia_AUC);
            PK2PK_rdia(data_row(t),a) = mean(RDia_pk_amp);
            resp_rate(data_row(t),a) = numel(rising_epochs_tongue)/((length(smooth_tongue_trim)/tongue_fs)/60);

            
            CTA_rdia(CTA_rdia == 0) = nan;
            avg_RDia_tonic(avg_RDia_tonic == 0) = nan;
            AUC_rdia(AUC_rdia == 0) = nan;
            PK2PK_rdia(PK2PK_rdia == 0) = nan;
            resp_rate(resp_rate == 0) = nan;

            
            %% simple threshold crossing using MPP line
            
            vars_o_interest = ['smooth_ldia_trim'; 'smooth_rdia_trim'];
            for m = 1:size(vars_o_interest,1)
                
                muscle_o_interest = eval(vars_o_interest(m,:));
                
                [PKS_RecMedSm, LOCS_RecMedSm] = findpeaks(muscle_o_interest,SR,'MinPeakDistance',1.2);
                PeakBurstRecMedSm = mean(PKS_RecMedSm);
                
                thresh = ((PeakBurstRecMedSm - min(muscle_o_interest))*0.5)+min(muscle_o_interest);
                high_pos = find(muscle_o_interest > thresh);
                high_diff = diff(high_pos);
                gaps = find(high_diff > 50*ms)+1;
                Rising_epochs = high_pos([1; gaps]);
                Rising_epochs = Rising_epochs(5:end-5);
                Falling_epochs =  high_pos([gaps-1;length(high_pos)]);
                Falling_epochs = Falling_epochs(5:end-5);
                Falling_epochs = Falling_epochs(Falling_epochs > Rising_epochs(1));
                Rising_epochs = Rising_epochs(Rising_epochs < Falling_epochs(end));
                
                %% calc first differential and first inflection points
                first_diff = diff(muscle_o_interest);
                
                thresh = 0;
                high_pos = find(first_diff > thresh);
                high_diff = diff(high_pos);
                gaps = find(high_diff > 25*ms)+1;
                FirstDeriv_Rising_epochs = high_pos([1; gaps]);
                FirstDeriv_Rising_epochs = FirstDeriv_Rising_epochs(5:end-5);
                FirstDeriv_Falling_epochs =  high_pos([gaps-1;length(high_pos)]);
                FirstDeriv_Falling_epochs = FirstDeriv_Falling_epochs(5:end-5);
                FirstDeriv_Falling_epochs = FirstDeriv_Falling_epochs(FirstDeriv_Falling_epochs > FirstDeriv_Rising_epochs(1));
                FirstDeriv_Rising_epochs = FirstDeriv_Rising_epochs(FirstDeriv_Rising_epochs < FirstDeriv_Falling_epochs(end));
                
                %% find closest inflection point to the simple threshold crossings
                New_RisingEpochs = nan(1,numel(Rising_epochs));
                New_FallingEpochs = nan(1,numel(Falling_epochs));
                
                for e = 1:numel(Rising_epochs)
                    temp = find(FirstDeriv_Rising_epochs < Rising_epochs(e));
                    New_RisingEpochs(e) = FirstDeriv_Rising_epochs(temp(end));
                    
                    temp = find(FirstDeriv_Rising_epochs > Falling_epochs(e));
                    New_FallingEpochs(e) = FirstDeriv_Rising_epochs(temp(1));
                end
                
                
                %% calculate Ti, Te,
                for e = 1:numel(New_RisingEpochs)-1
                    Ti_vec(e) = (New_FallingEpochs(e)-New_RisingEpochs(e))/SR;
                    Te_vec(e) = (New_RisingEpochs(e+1)-New_FallingEpochs(e))/SR;
                end
                
                if strcmp(vars_o_interest(m,:), 'smooth_ldia_trim')
                    Avg_Ti_L_Dia(data_row(t),a) = mean(Ti_vec);
                    Avg_Te_L_Dia(data_row(t),a) = mean(Te_vec);
                    
                    Avg_Ti_L_Dia(Avg_Ti_L_Dia == 0) = nan;
                    Avg_Te_L_Dia(Avg_Te_L_Dia == 0) = nan;
                elseif strcmp(vars_o_interest(m,:), 'smooth_rdia_trim')
                    Avg_Ti_R_Dia(data_row(t),a) = mean(Ti_vec);
                    Avg_Te_R_Dia(data_row(t),a) = mean(Te_vec);
                    
                    Avg_Ti_R_Dia(Avg_Ti_R_Dia == 0) = nan;
                    Avg_Te_R_Dia(Avg_Te_R_Dia == 0) = nan;
                else
                end
            end
        end
    end
end

%% normalize data to baseline to individual baseline
for a = 1:numel(animals)
    norm_fact_tongue = AUC_tongue(1,a);
    norm_AUC_tongue(:,a) = (AUC_tongue(:,a)./norm_fact_tongue)*100;
    
    norm_fact_ldia = AUC_ldia(1,a);
    norm_AUC_ldia(:,a) = (AUC_ldia(:,a)./norm_fact_ldia)*100;
    
    norm_fact_rdia = AUC_rdia(1,a);
    norm_AUC_rdia(:,a) = (AUC_rdia(:,a)./norm_fact_rdia)*100;
    
    norm_fact_tongue_pk2pk = PK2PK_tongue(1,a);
    norm_PK2PK_tongue(:,a) = (PK2PK_tongue(:,a)./norm_fact_tongue_pk2pk)*100;
    
    norm_fact_pk2pk_ldia = PK2PK_ldia(1,a);
    norm_PK2PK_ldia(:,a) = (PK2PK_ldia(:,a)./norm_fact_pk2pk_ldia)*100;
    
    norm_fact_pk2pk_rdia = PK2PK_rdia(1,a);
    norm_PK2PK_rdia(:,a) = (PK2PK_rdia(:,a)./norm_fact_pk2pk_rdia)*100;
    
    norm_fact_tonic_rdia = avg_RDia_tonic(1,a);
    norm_tonic_rdia(:,a) = (avg_RDia_tonic(:,a)./norm_fact_tonic_rdia)*100;
    
    norm_fact_tonic_ldia = avg_LDia_tonic(1,a);
    norm_tonic_ldia(:,a) = (avg_LDia_tonic(:,a)./norm_fact_tonic_ldia)*100;
         
end


%% find max left and right diaphragm avg response for each animal
[norm_left_AUC_max, left_I] = max(norm_AUC_ldia(3:4,:));

[norm_right_AUC_max, right_I] = max(norm_AUC_rdia(3:4,:));


left_right_max_AUC = [norm_left_AUC_max; norm_right_AUC_max];

[overall_max_AUC, overall_I] = max(left_right_max_AUC);


for o = 1:numel(overall_I)
    if overall_I(o) == 1
        AUC_max_saline(o) = norm_AUC_ldia(1,o);
    elseif overall_I(o) == 2
        AUC_max_saline(o) = norm_AUC_rdia(1,o);
    else
    end
end


norm_left_PK2PK_max = max(norm_PK2PK_ldia(2:4,:));
norm_right_PK2PK_max = max(norm_PK2PK_rdia(2:4,:));

left_right_max_PK2PK = [norm_left_PK2PK_max;norm_right_PK2PK_max];

overall_max_PK2PK = max(left_right_max_PK2PK);

%% save important variables
if strcmp(save_data, 'yes')
    if strcmp(strain, "Chat-Cre")
        animals_chat = animals;
        
        AUC_ldia_chat = AUC_ldia;
        AUC_rdia_chat = AUC_rdia;
        
        PK2PK_ldia_chat = PK2PK_ldia;
        PK2PK_rdia_chat = PK2PK_rdia;
        
        tonic_rdia_chat = avg_RDia_tonic;
        tonic_ldia_chat = avg_LDia_tonic;
        
        norm_AUC_ldia_chat = norm_AUC_ldia;
        norm_AUC_rdia_chat = norm_AUC_rdia;
        
        norm_PK2PK_ldia_chat = norm_PK2PK_ldia;
        norm_PK2PK_rdia_chat = norm_PK2PK_rdia;
        
        norm_tonic_rdia_chat = norm_tonic_rdia;
        norm_tonic_ldia_chat = norm_tonic_ldia;
        
        AUC_max_saline_chat = AUC_max_saline;
        overall_max_AUC_chat = overall_max_AUC;
        
        resp_rate_chat = resp_rate;
        
        Avg_Ti_L_Dia_chat = Avg_Ti_L_Dia;
        Avg_Ti_R_Dia_chat = Avg_Ti_R_Dia;
        Avg_Te_L_Dia_chat = Avg_Te_L_Dia;
        Avg_Te_R_Dia_chat = Avg_Te_R_Dia;
        
        cd('')  % enter directory where data will be saved
        save('chat_cre_mouse_data.mat', 'animals_chat', 'SR', 'ms', 'AUC_ldia_chat', 'AUC_rdia_chat', 'PK2PK_ldia_chat',...
            'PK2PK_rdia_chat', 'tonic_ldia_chat', 'tonic_rdia_chat','norm_AUC_ldia_chat' ,'norm_AUC_rdia_chat',...
            'norm_PK2PK_ldia_chat', 'norm_PK2PK_rdia_chat','norm_tonic_rdia_chat', 'norm_tonic_ldia_chat',...
            'AUC_max_saline_chat', 'overall_max_AUC_chat', 'resp_rate_chat','Avg_Ti_L_Dia_chat','Avg_Ti_R_Dia_chat','Avg_Te_L_Dia_chat', 'Avg_Te_R_Dia_chat')
    elseif strcmp(strain, "WT")
        animals_wt = animals;
  
        AUC_ldia_wt = AUC_ldia;
        AUC_rdia_wt = AUC_rdia;
        
        PK2PK_ldia_wt = PK2PK_ldia;
        PK2PK_rdia_wt = PK2PK_rdia;
        
        tonic_rdia_wt = avg_RDia_tonic;
        tonic_ldia_wt = avg_LDia_tonic;
        
        norm_AUC_ldia_wt = norm_AUC_ldia;
        norm_AUC_rdia_wt = norm_AUC_rdia;

        norm_PK2PK_ldia_wt = norm_PK2PK_ldia;
        norm_PK2PK_rdia_wt = norm_PK2PK_rdia;
        
        norm_tonic_rdia_wt = norm_tonic_rdia;
        norm_tonic_ldia_wt = norm_tonic_ldia;
        
        AUC_max_saline_wt = AUC_max_saline;
        overall_max_AUC_wt = overall_max_AUC;
        
        resp_rate_wt = resp_rate;
        
        Avg_Ti_L_Dia_wt = Avg_Ti_L_Dia;
        Avg_Ti_R_Dia_wt = Avg_Ti_R_Dia;
        Avg_Te_L_Dia_wt = Avg_Te_L_Dia;
        Avg_Te_R_Dia_wt = Avg_Te_R_Dia;
        
        cd('') % enter directory where data will be saved
        save('wt_mouse_data.mat', 'animals_wt', 'SR', 'ms', 'norm_AUC_ldia_wt', 'norm_AUC_rdia_wt', 'AUC_ldia_wt',...
            'AUC_rdia_wt', 'PK2PK_ldia_wt', 'PK2PK_rdia_wt', 'tonic_rdia_wt', 'tonic_ldia_wt', 'norm_PK2PK_ldia_wt',...
            'norm_PK2PK_rdia_wt','norm_tonic_rdia_wt', 'norm_tonic_ldia_wt', 'AUC_max_saline_wt', 'overall_max_AUC_wt',...
            'resp_rate_wt','Avg_Ti_L_Dia_wt','Avg_Ti_R_Dia_wt','Avg_Te_L_Dia_wt', 'Avg_Te_R_Dia_wt')
        
    else
    end
    
elseif strcmp(save_data, 'no')
else
end

