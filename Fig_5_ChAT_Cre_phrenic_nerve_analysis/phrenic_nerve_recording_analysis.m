clear all
close all
clc

%%
root = '~\Fig_5_ChAT_Cre_phrenic_nerve_analysis\data'; %change to directory where data is located
cd(root)

animals = dir;
animals = animals(3:end); %ignores hidden files
down_sample = 'no';
disp_test_figs = 'no';

for an = 1:numel(animals)
    
    an_name = animals(an).name;
    cd([root '\' an_name '\mat_export'])
    
    files = dir;
    files = files(3:end);
    for f = 1:numel(files)-1
        clear lpks llocs rpks rlocs SBP BeatTiming DBP
        
        fname = files(f).name;
        load(fname)
        
        %% define standard variables
        phrenic_Fs = 1./LPhr_R.interval;
        ABP_Fs = 1./Art_Pres.interval;
        ms = phrenic_Fs/1000;
        
        if strcmp(down_sample, 'yes')
            %% downsample raw data
            downsample_fact = 10;
            
            lphr_raw_vals = downsample(LPhr_R.values,downsample_fact);
            rphr_raw_vals = downsample(RPhr_R.values,downsample_fact);
            
            clear LPhr_R.values RPhr_R.values
            
            new_Fs = phrenic_Fs./downsample_fact;
            
            phrenic_Fs = new_Fs;
            
            ms = phrenic_Fs/1000;
            
        elseif strcmp(down_sample, 'no')
            lphr_raw_vals = LPhr_R.values;
            rphr_raw_vals = RPhr_R.values;
            
        else
        end
        
        %% filter phrenic traces
        [B,A] = butter(2,[10 .95*(phrenic_Fs/2)]/(phrenic_Fs/2), "bandpass");
        
        filt_lphr = filtfilt(B,A,lphr_raw_vals);
        filt_rphr = filtfilt(B,A,rphr_raw_vals);
        
        %% rectify and integrate
        rec_lphr = abs(filt_lphr);
        rec_rphr = abs(filt_rphr);
        
        med_lphr = movmedian(rec_lphr, 50*ms);
        med_rphr = movmedian(rec_rphr, 50*ms);
        
        medsmooth_lphr = smooth(med_lphr,50*ms);
        medsmooth_rphr = smooth(med_rphr,50*ms);
        
        %% identify trace portions of interest
        % blood gases sampled at 20, 40, 60, 80, 100 mins 
        
        if f == 1
            times_o_interest = [419 599];
            data_row = 1;
        elseif f == 2
            times_o_interest = [419 599];
            data_row = 2;
        elseif f == 3
            times_o_interest = [120 300; 720 900; 1620 1800; 3420 3600;...
                5220 5400];
            data_row = [3 4 5 6 7];
        elseif f == 4
            if length(medsmooth_lphr) >= length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_rphr)/phrenic_Fs)-180 (length(medsmooth_rphr)/phrenic_Fs)];
            elseif length(medsmooth_lphr) < length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_lphr)/phrenic_Fs)-180 (length(medsmooth_lphr)/phrenic_Fs)];
            else
                keyboard
            end
            
            data_row = 8;
        elseif f == 5
            if length(medsmooth_lphr) >= length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_rphr)/phrenic_Fs)-180 (length(medsmooth_rphr)/phrenic_Fs)];
            elseif length(medsmooth_lphr) < length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_lphr)/phrenic_Fs)-180 (length(medsmooth_lphr)/phrenic_Fs)];
            else
                keyboard
            end
            data_row = 9;
        elseif f == 6
            if length(medsmooth_lphr) >= length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_rphr)/phrenic_Fs)-180 (length(medsmooth_rphr)/phrenic_Fs)];
            elseif length(medsmooth_lphr) < length(medsmooth_rphr)
                times_o_interest = [(length(medsmooth_lphr)/phrenic_Fs)-180 (length(medsmooth_lphr)/phrenic_Fs)];
            else
                keyboard
            end
            data_row = 10;
        else
            
        end
        
        %% loop through times of interest
        for t = 1:size(times_o_interest,1)
            
            smooth_lphr_trim = medsmooth_lphr(floor(times_o_interest(t,1)*phrenic_Fs): floor(times_o_interest(t,2)*phrenic_Fs));
            smooth_rphr_trim = medsmooth_rphr(floor(times_o_interest(t,1)*phrenic_Fs): floor(times_o_interest(t,2)*phrenic_Fs));
            
            Art_Pres_trim = Art_Pres.values(floor(times_o_interest(t,1)*ABP_Fs): floor(times_o_interest(t,2)*ABP_Fs));
            
            %% calculate suggested MPPs
            mult_fact = 0.25;
            
            [N,X] = hist(smooth_lphr_trim,1000);
            N = smooth(N,100);
            [~,LOCS] = findpeaks(N,'MinPeakProminence',1);
            RecMedSm_sug_MPP_left = X(LOCS(end))*mult_fact;
            
            [N,X] = hist(smooth_rphr_trim,1000);
            N = smooth(N,100);
            [~,LOCS] = findpeaks(N,'MinPeakProminence',1);
            RecMedSm_sug_MPP_right = X(LOCS(end))*mult_fact;
            
            % adjust mean peak prominance if calculated values does not
            % accurately detect peaks
            if strcmp(fname, "01_28_21AP_20230508_BL.mat") && data_row(t) == 1
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.8; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "01_29_21AP_20230801_BL.mat") && data_row(t) == 1
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.7; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_33_21AP_20230816_j60.mat") && data_row(t) == 7
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.8; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.8; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_33_21AP_20230816_j60.mat") && data_row(t) == 6
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.6;
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_33_21AP_20230816_j60.mat") && data_row(t) == 5
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.5; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1.1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_33_21AP_20230816_j60.mat") && data_row(t) == 4
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.8; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_33_21AP_20230816_j60.mat") && data_row(t) == 3
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.5;
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "01_33_21AP_20230816_BL.mat") && data_row(t) == 1
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1;
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.7; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_30_21AP_20230728_j60.mat") && data_row(t) == 7
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.5; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.3; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_30_21AP_20230728_j60.mat") && data_row(t) == 6
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.7; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_30_21AP_20230728_j60.mat") && data_row(t) == 5
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.1; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.9; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_30_21AP_20230728_j60.mat") && data_row(t) == 3
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.75; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1.25; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "02_30_21AP_20230728_saline.mat") && data_row(t) == 2
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.75; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_34_20230726_j60.mat") && data_row(t) == 7
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.25; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1.25; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_34_20230726_j60.mat") && data_row(t) == 4
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.25; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1.25; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_34_20230726_j60.mat") && data_row(t) == 3
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1; %good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.8; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "03_34_20230726_j60.mat") && data_row(t) == 6
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1.1; % good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1.35; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "04_30_21AP_20230728_BLhypoxia.mat") && data_row(t) == 8
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.2; % good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.4; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "05_28_21AP_20230508_j60hypoxia.mat") && data_row(t) == 9
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.8; % good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*1; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "04_29_21AP_20230801_BLhypoxia.mat") && data_row(t) == 8
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*1; % good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.6; %good
                disp([fname '_' num2str(data_row(t))])
            elseif strcmp(fname, "05_33_21AP_20230816_j60hypoxia.mat") && data_row(t) == 9
                RecMedSm_sug_MPP_left = RecMedSm_sug_MPP_left*0.5; % good
                RecMedSm_sug_MPP_right = RecMedSm_sug_MPP_right*0.5; %good
                disp([fname '_' num2str(data_row(t))])
            else
            end
            %% avg peak amplitude
            [lpks, llocs] = findpeaks(smooth_lphr_trim, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_left,'MinPeakDistance',1);
            [rpks, rlocs] = findpeaks(smooth_rphr_trim, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_right,'MinPeakDistance',1);
            
            if strcmp(disp_test_figs, 'yes')
                figure('Name', [fname ' data row: ' num2str(data_row(t))])
                subplot(2,1,1)
                findpeaks(smooth_lphr_trim, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_left,'MinPeakDistance',1)
                title('Left phrenic')
                subplot(2,1,2)
                findpeaks(smooth_rphr_trim, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_right,'MinPeakDistance',1)
                title('Right phrenic')
                linkaxes
                set(gcf,'pos',[ 1          41        1920         963])
                
            elseif strcmp(disp_test_figs, 'no')
            else
            end
            
            avg_left_pk_amp(data_row(t),an) = mean(lpks);
            avg_right_pk_amp(data_row(t),an) = mean(rpks);
            
            
            %% avg tonic
            [l_localmin, llocsmin] = findpeaks(smooth_lphr_trim*-1, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_left,'MinPeakDistance',1);
            [r_localmin, rlocsmin] = findpeaks(smooth_rphr_trim*-1, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_right,'MinPeakDistance',1);
            
            if strcmp(disp_test_figs, 'yes')
                figure('Name', [fname ' tonic'])
                subplot(2,1,1)
                findpeaks(smooth_lphr_trim*-1, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_left,'MinPeakDistance',1)
                title('Left phrenic Tonic')
                subplot(2,1,2)
                findpeaks(smooth_rphr_trim*-1, phrenic_Fs,'MinPeakProminence',RecMedSm_sug_MPP_right,'MinPeakDistance',1)
                title('Right phrenic Tonic')
                linkaxes
                set(gcf,'pos',[ 1          41        1920         963])
                
            elseif strcmp(disp_test_figs, 'no')
            else
            end
            
            avg_left_tonic(data_row(t),an) = mean((l_localmin*-1));
            avg_right_tonic(data_row(t),an) = mean((r_localmin*-1));
            
            %% instanteous respiratory rate
            RespRate_l(data_row(t),an) = numel(llocs)/((llocs(end)-llocs(1))/60);
            RespRate_r(data_row(t),an) = numel(rlocs)/((rlocs(end)-rlocs(1))/60);
            
            %% mean arterial blood pressure
            % calculate suggested MPPs
            MPP = 35;
            PkPk_90 = ((max(Art_Pres_trim)*0.95)-(min(Art_Pres_trim)*1.05));
            if PkPk_90 < MPP
                MPP = PkPk_90*0.8;
            else
            end
            
            [SBP, BeatTiming] = findpeaks(Art_Pres_trim,ABP_Fs,'MinPeakProminence',MPP,'MinPeakDistance',0.07);
            avg_SBP(data_row(t),an) = mean(SBP);
            % figure(2220)
            % findpeaks(Art_Pres.values,ABP_Fs,'MinPeakProminence',MPP)
            HR(data_row(t),an) = numel(BeatTiming)/((BeatTiming(end)-BeatTiming(1))/60);
            [DBP] = findpeaks(-Art_Pres_trim,ABP_Fs,'MinPeakProminence',MPP);
            avg_DBP(data_row(t),an) = mean(-DBP);
            MAP(data_row(t),an)= avg_DBP(data_row(t),an) + 1/3*(avg_SBP(data_row(t),an)-avg_DBP(data_row(t),an));
            
        end
    end
end


%% plotting raw values

%pk amplitudes
mean_left_pk_amp = mean(avg_left_pk_amp(1:7,:),2);
sem_left_pk_amp = std(avg_left_pk_amp(1:7,:),[],2)/sqrt(size(avg_left_pk_amp(1:7,:),2));

mean_right_pk_amp = mean(avg_right_pk_amp(1:7,:),2);
sem_right_pk_amp = std(avg_right_pk_amp(1:7,:),[],2)/sqrt(size(avg_right_pk_amp(1:7,:),2));
 
%tonic activity
mean_left_tonic = mean(avg_left_tonic(1:7,:),2);
sem_left_tonic = std(avg_left_tonic(1:7,:),[],2)/sqrt(size(avg_left_tonic(1:7,:),2));

mean_right_tonic = mean(avg_right_tonic(1:7,:),2);
sem_right_tonic = std(avg_right_tonic(1:7,:),[],2)/sqrt(size(avg_right_tonic(1:7,:),2));

%blood pressure parameters
mean_HR = mean(HR(1:7,:),2);
mean_DBP = mean(avg_DBP(1:7,:),2);
mean_SBP = mean(avg_SBP(1:7,:),2);
mean_MAP = mean(MAP(1:7,:),2);

sem_HR = std(HR(1:7,:),[],2)/sqrt(size(HR(1:7,:),2));
sem_DBP = std(avg_DBP(1:7,:),[],2)/sqrt(size(avg_DBP(1:7,:),2));
sem_SBP = std(avg_SBP(1:7,:),[],2)/sqrt(size(avg_SBP(1:7,:),2));
sem_MAP = std(MAP(1:7,:),[],2)/sqrt(size(MAP(1:7,:),2));

% respiratory rate
mean_left_rr = mean(RespRate_l(1:7,:),2);
sem_left_rr = std(RespRate_l(1:7,:),[],2)/sqrt(size(RespRate_l(1:7,:),2));

mean_right_rr = mean(RespRate_r(1:7,:),2);
sem_right_rr = std(RespRate_r(1:7,:),[],2)/sqrt(size(RespRate_r(1:7,:),2));

all_resp_rate(:,:,1) = RespRate_l(1:7,:);
all_resp_rate(:,:,2) = RespRate_r(1:7,:);

avg_resp_rate = mean(mean(all_resp_rate,3),2);
sd_resp_rate = std(mean(all_resp_rate,3),[],2);

sem_resp_rate = sd_resp_rate/sqrt(9);

%pk amplitudes
norm_fact_l_pk = avg_left_pk_amp(1,:);
norm_left_pk_vals = (avg_left_pk_amp(1:7,:)./norm_fact_l_pk)*100;

norm_mean_left_pk_amp = mean(norm_left_pk_vals,2);
norm_sem_left_pk_amp = std(norm_left_pk_vals,[],2)/sqrt(size(norm_left_pk_vals,2));

norm_fact_r_pk = avg_right_pk_amp(1,:);
norm_right_pk_vals = (avg_right_pk_amp(1:7,:)./norm_fact_r_pk)*100;

norm_mean_right_pk_amp = mean(norm_right_pk_vals,2);
norm_sem_right_pk_amp = std(norm_right_pk_vals,[],2)/sqrt(size(norm_right_pk_vals,2));

%tonic activity
norm_fact_l_tonic = avg_left_tonic(1,:);
norm_left_tonic_vals = (avg_left_tonic(1:7,:)./norm_fact_l_tonic)*100;

norm_mean_left_tonic = mean(norm_left_tonic_vals,2);
norm_sem_left_tonic = std(norm_left_tonic_vals,[],2)/sqrt(size(norm_left_tonic_vals,2));

norm_fact_r_tonic = avg_right_tonic(1,:);
norm_right_tonic_vals = (avg_right_tonic(1:7,:)./norm_fact_r_tonic)*100;

norm_mean_right_tonic = mean(norm_right_tonic_vals,2);
norm_sem_right_tonic = std(norm_right_tonic_vals,[],2)/sqrt(size(norm_right_tonic_vals,2));

%% subpanel plot for pub
figure

subplot(3,3,1)
hold on
scatter(1:length(mean_left_pk_amp), mean_left_pk_amp, [], [0.8500 0.3250 0.0980],  'filled')
errorbar(1:length(mean_left_pk_amp), mean_left_pk_amp, sem_left_pk_amp, sem_left_pk_amp,  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
scatter(1:length(mean_right_pk_amp), mean_right_pk_amp, [], [0 0.4470 0.7410], 'filled')
errorbar(1:length(mean_right_pk_amp), mean_right_pk_amp, sem_right_pk_amp, sem_right_pk_amp, 'color', [0 0.4470 0.7410], 'CapSize', 12)
ylabel('Phrenic nerve peak-to-peak amplitude (a.u.)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])

subplot(3,3,2)
hold on
scatter(1:length(norm_mean_left_pk_amp(2:end)), norm_mean_left_pk_amp(2:end), [], [0.8500 0.3250 0.0980],  'filled')
errorbar(1:length(norm_mean_left_pk_amp(2:end)), norm_mean_left_pk_amp(2:end), norm_sem_left_pk_amp(2:end), norm_sem_left_pk_amp(2:end),  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
scatter(1:length(norm_mean_right_pk_amp(2:end)), norm_mean_right_pk_amp(2:end), [], [0 0.4470 0.7410], 'filled')
errorbar(1:length(norm_mean_right_pk_amp(2:end)), norm_mean_right_pk_amp(2:end), norm_sem_right_pk_amp(2:end), norm_sem_right_pk_amp(2:end), 'color', [0 0.4470 0.7410], 'CapSize', 12)
ylabel('Normalized Phrenic nerve peak-to-peak amplitude (% of baseline)')
set(gca,'XTick', [1:6], 'XTickLabel', {'SL', '5','15','30','60','90'})
xlim([0 7])

subplot(3,3,3)
hold on
scatter(1:length(mean_left_tonic), mean_left_tonic, [], [0.8500 0.3250 0.0980],  'filled')
errorbar(1:length(mean_left_tonic), mean_left_tonic, sem_left_tonic, sem_left_tonic,  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
scatter(1:length(mean_right_tonic), mean_right_tonic, [], [0 0.4470 0.7410], 'filled')
errorbar(1:length(mean_right_tonic), mean_right_tonic, sem_right_tonic, sem_right_tonic, 'color', [0 0.4470 0.7410], 'CapSize', 12)
ylabel('Phrenic nerve tonic activity (a.u.)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])

subplot(3,3,4)
hold on
scatter(1:length(norm_mean_left_tonic(2:end)), norm_mean_left_tonic(2:end), [], [0.8500 0.3250 0.0980],  'filled')
errorbar(1:length(norm_mean_left_tonic(2:end)), norm_mean_left_tonic(2:end), norm_sem_left_tonic(2:end), norm_sem_left_tonic(2:end),  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
scatter(1:length(norm_mean_right_tonic(2:end)), norm_mean_right_tonic(2:end), [], [0 0.4470 0.7410], 'filled')
errorbar(1:length(norm_mean_right_tonic(2:end)), norm_mean_right_tonic(2:end), norm_sem_right_tonic(2:end), norm_sem_right_tonic(2:end), 'color', [0 0.4470 0.7410], 'CapSize', 12)
ylabel('Phrenic nerve tonic activity (% of baseline)')
set(gca,'XTick', [1:6], 'XTickLabel', {'SL', '5','15','30','60','90'})
xlim([0 7])

subplot(3,3,5)
hold on

scatter(1:length(mean_SBP), mean_SBP, [], 'k',  'filled')
errorbar(1:length(mean_SBP), mean_SBP, sem_SBP, sem_SBP,  'color', 'k', 'CapSize', 12)
ylabel('Systolic blood pressure (mmHg)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([100 200])

subplot(3,3,6)
hold on

scatter(1:length(mean_DBP), mean_DBP, [], 'k',  'filled')
errorbar(1:length(mean_DBP), mean_DBP, sem_DBP, sem_DBP,  'color', 'k', 'CapSize', 12)
ylabel('Diastolic blood pressure (mmHg)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([0 100])


subplot(3,3,7)
hold on

scatter(1:length(mean_HR), mean_HR, [], 'k',  'filled')
errorbar(1:length(mean_HR), mean_HR, sem_HR, sem_HR,  'color', 'k', 'CapSize', 12)
ylabel('Heart rate (bpm)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([300 450])

subplot(3,3,8)
hold on

scatter(1:length(mean_MAP), mean_MAP, [], 'k',  'filled')
errorbar(1:length(mean_MAP), mean_MAP, sem_MAP, sem_MAP,  'color', 'k', 'CapSize', 12)
ylabel('Mean arterial blood pressure (mmHg)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([50 150])

subplot(3,3,8)
hold on

scatter(1:length(mean_MAP), mean_MAP, [], 'k',  'filled')
errorbar(1:length(mean_MAP), mean_MAP, sem_MAP, sem_MAP,  'color', 'k', 'CapSize', 12)
ylabel('Mean arterial blood pressure (mmHg)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([50 150])

subplot(3,3,9)
hold on

scatter(1:length(avg_resp_rate), avg_resp_rate, [], 'k', 'filled')
errorbar(1:length(avg_resp_rate), avg_resp_rate, sem_resp_rate, sem_resp_rate, 'color', 'k', 'CapSize', 12)
ylabel('Respiratory Rate (bpm)')
set(gca,'XTick', [1:7], 'XTickLabel', {'BL', 'SL', '5','15','30','60','90'})
xlim([0 8])
ylim([0 50])

































































