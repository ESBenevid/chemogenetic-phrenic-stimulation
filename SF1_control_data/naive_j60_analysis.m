clear all
close all
clc

cd('~\SF1_control_data\data'); % change directory to where data files are located

files = dir;
files = files(3:end);

for i = 1:numel(files)
    fname = files(i).name;
    load(fname)
    %% filter EMG
    Fs = 1/L_Dia.interval;
    ms = Fs/1000;
    
    [B,A] = butter(2, [100 1000]/(Fs/2), "bandpass");
    
    filt_L_Dia = filtfilt(B,A,L_Dia.values);
    filt_R_Dia = filtfilt(B,A,R_Dia.values);
    
    %% rec int
    rec_L_Dia = abs(filt_L_Dia);
    rec_R_Dia = abs(filt_R_Dia);
    
    med_L_Dia = movmedian(rec_L_Dia, 150*ms);
    med_R_Dia = movmedian(rec_R_Dia, 150*ms);
    
    medsmooth_L_Dia = smooth(med_L_Dia,150*ms);
    medsmooth_R_Dia = smooth(med_R_Dia,150*ms);
    
    
    %% peak detection
    % calculate suggested MPPs
    if i == 1
        mult_fact = 0.3;
    elseif i == 2
        mult_fact = 0.65;
    else
        mult_fact = 0.40;
    end
    
    [N,X] = hist(medsmooth_L_Dia,10000);
    N = smooth(N,1000);
    [~,LOCS] = findpeaks(N,'MinPeakProminence',1);
    sug_MPP_left = X(LOCS(end))*mult_fact;

    
    [N,X] = hist(medsmooth_R_Dia,10000);
    N = smooth(N,1000);
    [~,LOCS] = findpeaks(N,'MinPeakProminence',1);
    sug_MPP_right = X(LOCS(end))*mult_fact;

    
    [R_Dia_pks,R_Dia_locs] = findpeaks(medsmooth_R_Dia, 'MinPeakProminence', sug_MPP_right);
    [L_Dia_pks,L_Dia_locs] = findpeaks(medsmooth_L_Dia, 'MinPeakProminence', sug_MPP_left);
    
    mean_L_Dia(i) = nanmean(L_Dia_pks);
    mean_R_Dia(i) = nanmean(R_Dia_pks);
    
    clear R_Dia_pks R_Dia_locs L_Dia_pks L_Dia_locs filt_L_Dia filt_R_Dia rec_L_Dia rec_R_Dia med_L_Dia med_R_Dia medsmooth_L_Dia medsmooth_R_Dia
end


ldia_reshape = reshape(mean_L_Dia,5,3);
rdia_reshape = reshape(mean_R_Dia,5,3);

norm_fact_l = ldia_reshape(:,1);
norm_fact_r = rdia_reshape(:,1);

norm_ldia = (ldia_reshape./norm_fact_l)*100;
norm_rdia = (rdia_reshape./norm_fact_r)*100;


figure
subplot(2,2,1)
hold on
bar(1,nanmean(ldia_reshape(:,1)),'Facecolor','w', 'EdgeColor','k')
errorbar(1,nanmean(ldia_reshape(:,1)),[],nanstd(ldia_reshape(:,1)),'k','Capsize',15)
scatter(repmat(1,5,1), ldia_reshape(:,1),[],'k','filled','jitter','on')
bar(2,nanmean(ldia_reshape(:,2)),'Facecolor','w', 'EdgeColor','k')
errorbar(2,nanmean(ldia_reshape(:,2)),[],nanstd(ldia_reshape(:,2)),'k','Capsize',15)
scatter(repmat(2,5,1), ldia_reshape(:,2),[],'k','filled','jitter','on')
bar(3,nanmean(ldia_reshape(:,3)),'Facecolor','w', 'EdgeColor','k')
errorbar(3,nanmean(ldia_reshape(:,3)),[],nanstd(ldia_reshape(:,3)),'k','Capsize',15)
scatter(repmat(3,5,1), ldia_reshape(:,3),[],'k','filled','jitter','on')
xlim([0.5 3.5])
ylabel('Left Dia Peak amplitude (mV)')
set(gca,'XTick', [1:3], 'XTickLabels', {'Baseline','Saline', 'J60'})

subplot(2,2,2)
hold on
bar(1,nanmean(rdia_reshape(:,1)),'Facecolor','w', 'EdgeColor','k')
errorbar(1,nanmean(rdia_reshape(:,1)),[],nanstd(rdia_reshape(:,1)),'k','Capsize',15)
scatter(repmat(1,5,1), rdia_reshape(:,1),[],'k','filled','jitter','on')
bar(2,nanmean(rdia_reshape(:,2)),'Facecolor','w', 'EdgeColor','k')
errorbar(2,nanmean(rdia_reshape(:,2)),[],nanstd(rdia_reshape(:,2)),'k','Capsize',15)
scatter(repmat(2,5,1), rdia_reshape(:,2),[],'k','filled','jitter','on')
bar(3,nanmean(rdia_reshape(:,3)),'Facecolor','w', 'EdgeColor','k')
errorbar(3,nanmean(rdia_reshape(:,3)),[],nanstd(rdia_reshape(:,3)),'k','Capsize',15)
scatter(repmat(3,5,1), rdia_reshape(:,3),[],'k','filled','jitter','on')
xlim([0.5 3.5])
ylabel('Right Dia Peak amplitude (mV)')
set(gca,'XTick', [1:3], 'XTickLabels', {'Baseline','Saline', 'J60'})

subplot(2,2,3)
hold on
bar(1,nanmean(norm_ldia(:,2)),'Facecolor','w', 'EdgeColor','k')
errorbar(1,nanmean(norm_ldia(:,2)),[],nanstd(norm_ldia(:,2)),'k','Capsize',15)
scatter(repmat(1,5,1), norm_ldia(:,2),[],'k','filled','jitter','on')
bar(2,nanmean(norm_ldia(:,3)),'Facecolor','w', 'EdgeColor','k')
errorbar(2,nanmean(norm_ldia(:,3)),[],nanstd(norm_ldia(:,3)),'k','Capsize',15)
scatter(repmat(2,5,1), norm_ldia(:,3),[],'k','filled','jitter','on')
xlim([0.5 2.5])
ylabel('Left Dia Peak amplitude (% BL)')
set(gca,'XTick', [1:2], 'XTickLabels', {'Saline', 'J60'})

subplot(2,2,4)
hold on
bar(1,nanmean(norm_rdia(:,2)),'Facecolor','w', 'EdgeColor','k')
errorbar(1,nanmean(norm_rdia(:,2)),[],nanstd(norm_rdia(:,2)),'k','Capsize',15)
scatter(repmat(1,5,1), norm_rdia(:,2),[],'k','filled','jitter','on')
bar(2,nanmean(norm_rdia(:,3)),'Facecolor','w', 'EdgeColor','k')
errorbar(2,nanmean(norm_rdia(:,3)),[],nanstd(norm_rdia(:,3)),'k','Capsize',15)
scatter(repmat(2,5,1), norm_rdia(:,3),[],'k','filled','jitter','on')
xlim([0.5 2.5])
ylabel('Right Dia Peak amplitude (% BL)')
set(gca,'XTick', [1:2], 'XTickLabels', {'Saline', 'J60'})



