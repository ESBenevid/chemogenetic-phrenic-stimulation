clear all
close all
clc

%%
cd('~\Fig_1_2_3_WT_ChAT_Cre_mice\data_to_generate_figs')% change to directory where data files are located


load('chat_cre_mouse_data')
load('wt_mouse_data')
%% find max left and right diaphragm avg response for each animal
animals = {animals_chat.name animals_wt.name};

% raw vals
AUC_ldia{1} = AUC_ldia_chat;
AUC_ldia{2} = AUC_ldia_wt;

AUC_rdia{1} = AUC_rdia_chat;
AUC_rdia{2} = AUC_rdia_wt;

PK2PK_ldia{1} = PK2PK_ldia_chat;
PK2PK_ldia{2} = PK2PK_ldia_wt;

PK2PK_rdia{1} = PK2PK_rdia_chat;
PK2PK_rdia{2} = PK2PK_rdia_wt;

tonic_ldia{1} = tonic_ldia_chat;
tonic_ldia{2} = tonic_ldia_wt;

tonic_rdia{1} = tonic_rdia_chat;
tonic_rdia{2} = tonic_rdia_wt;

resp_rate{1} = resp_rate_chat;
resp_rate{2} = resp_rate_wt;

Avg_Te_ldia{1} = Avg_Te_L_Dia_chat;
Avg_Te_ldia{2} = Avg_Te_L_Dia_wt;

Avg_Te_rdia{1} = Avg_Te_R_Dia_chat;
Avg_Te_rdia{2} = Avg_Te_R_Dia_wt;

Avg_Ti_ldia{1} = Avg_Ti_L_Dia_chat;
Avg_Ti_ldia{2} = Avg_Ti_L_Dia_wt;

Avg_Ti_rdia{1} = Avg_Ti_R_Dia_chat;
Avg_Ti_rdia{2} = Avg_Ti_R_Dia_wt;


% normalized vars
norm_AUC_ldia{1} = norm_AUC_ldia_chat;
norm_AUC_ldia{2} = norm_AUC_ldia_wt;

norm_AUC_rdia{1} = norm_AUC_rdia_chat;
norm_AUC_rdia{2} = norm_AUC_rdia_wt;

norm_PK2PK_ldia{1} = norm_PK2PK_ldia_chat;
norm_PK2PK_ldia{2} = norm_PK2PK_ldia_wt;

norm_PK2PK_rdia{1} = norm_PK2PK_rdia_chat;
norm_PK2PK_rdia{2} = norm_PK2PK_rdia_wt;

norm_tonic_ldia{1} = norm_tonic_ldia_chat;
norm_tonic_ldia{2} = norm_tonic_ldia_wt;

norm_tonic_rdia{1} = norm_tonic_rdia_chat;
norm_tonic_rdia{2} = norm_tonic_rdia_wt;

%% Plot normalized mouse data
fig_names = {'Chat-Cre mice norm','WT mice norm'};
for i = 1:numel(fig_names)
    % AUC
    avg_norm_AUC_ldia = nanmean(norm_AUC_ldia{i}(2:end,:),2);
    avg_norm_AUC_rdia = nanmean(norm_AUC_rdia{i}(2:end,:),2);
    SEM_norm_AUC_ldia = (nanstd(norm_AUC_ldia{i}(2:end,:),[],2)/sqrt(size(norm_AUC_ldia{i}(2:end,:),2)));
    SEM_norm_AUC_rdia = (nanstd(norm_AUC_rdia{i}(2:end,:),[],2)/sqrt(size(norm_AUC_rdia{i}(2:end,:),2)));
    
    figure('Name', fig_names{i})
    subplot(2,2,1)
    hold on
    scatter(1:length(avg_norm_AUC_ldia), avg_norm_AUC_ldia, 50, [0.8500 0.3250 0.0980],  'filled')
    errorbar(1:length(avg_norm_AUC_ldia), avg_norm_AUC_ldia, SEM_norm_AUC_ldia, SEM_norm_AUC_ldia,  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
    scatter(1:length(avg_norm_AUC_rdia), avg_norm_AUC_rdia, 50, [0 0.4470 0.7410], 'filled')
    errorbar(1:length(avg_norm_AUC_rdia), avg_norm_AUC_rdia, SEM_norm_AUC_rdia, SEM_norm_AUC_rdia, 'color', [0 0.4470 0.7410], 'CapSize', 12)
    
    plot([0 7], [100 100], 'k--')
    ylabel('Diaphragm EMG AUC (% of baseline)')
    set(gca,'XTick', [1:6], 'XTickLabel', {'SL' '5','15','30','60','90'})
    xlim([0 7])
    ylim([60 220])
    
    
    % PK2PK amp
    avg_norm_PK2PK_ldia = nanmean(norm_PK2PK_ldia{i}(2:end,:),2);
    avg_norm_PK2PK_rdia = nanmean(norm_PK2PK_rdia{i}(2:end,:),2);
    SEM_norm_PK2PK_ldia = (nanstd(norm_PK2PK_ldia{i}(2:end,:),[],2)/sqrt(size(norm_PK2PK_ldia{i}(2:end,:),2)));
    SEM_norm_PK2PK_rdia = (nanstd(norm_PK2PK_rdia{i}(2:end,:),[],2)/sqrt(size(norm_PK2PK_rdia{i}(2:end,:),2)));
    
    subplot(2,2,2)
    hold on
    scatter(1:length(avg_norm_PK2PK_ldia), avg_norm_PK2PK_ldia, 50, [0.8500 0.3250 0.0980],  'filled')
    errorbar(1:length(avg_norm_PK2PK_ldia), avg_norm_PK2PK_ldia, SEM_norm_PK2PK_ldia, SEM_norm_PK2PK_ldia,  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
    scatter(1:length(avg_norm_PK2PK_rdia), avg_norm_PK2PK_rdia, 50, [0 0.4470 0.7410], 'filled')
    errorbar(1:length(avg_norm_PK2PK_rdia), avg_norm_PK2PK_rdia, SEM_norm_PK2PK_rdia, SEM_norm_PK2PK_rdia, 'color', [0 0.4470 0.7410], 'CapSize', 12)
    
    plot([0 7], [100 100], 'k--')
    ylabel('Diaphragm Peak to Peak amplitude (% of baseline)')
    set(gca,'XTick', [1:6], 'XTickLabel', {'SL' '5','15','30','60','90'})
    xlim([0 7])
    ylim([70 185])
    
    % Tonic
    avg_norm_tonic_ldia = nanmean(norm_tonic_ldia{i}(2:end,:),2);
    avg_norm_tonic_rdia = nanmean(norm_tonic_rdia{i}(2:end,:),2);
    SEM_norm_tonic_ldia = (nanstd(norm_tonic_ldia{i}(2:end,:),[],2)/sqrt(size(norm_tonic_ldia{i}(2:end,:),2)));
    SEM_norm_tonic_rdia = (nanstd(norm_tonic_rdia{i}(2:end,:),[],2)/sqrt(size(norm_tonic_rdia{i}(2:end,:),2)));
    
    subplot(2,2,3)
    hold on
    scatter(1:length(avg_norm_tonic_ldia), avg_norm_tonic_ldia, 50, [0.8500 0.3250 0.0980],  'filled')
    errorbar(1:length(avg_norm_tonic_ldia), avg_norm_tonic_ldia, SEM_norm_tonic_ldia, SEM_norm_tonic_ldia,  'color', [0.8500 0.3250 0.0980], 'CapSize', 12)
    scatter(1:length(avg_norm_tonic_rdia), avg_norm_tonic_rdia, 50, [0 0.4470 0.7410], 'filled')
    errorbar(1:length(avg_norm_tonic_rdia), avg_norm_tonic_rdia, SEM_norm_tonic_rdia, SEM_norm_tonic_rdia, 'color', [0 0.4470 0.7410], 'CapSize', 12)

    plot([0 7], [100 100], 'k--')
    ylabel('Diaphragm Tonic activity (% of baseline)')
    set(gca,'XTick', [1:6], 'XTickLabel', {'SL' '5','15','30','60','90'})
    xlim([0 7])
    ylim([50 575])
    
    % Resp Rate
    resp_rate_new = nanmean(resp_rate{i},2);
    SEM_resp_rate_new = (nanstd(resp_rate{i},[],2)/sqrt(size(resp_rate{i},2)));
    
    subplot(2,2,4)
    hold on
    scatter(1:length(resp_rate_new), resp_rate_new, [], [0 0 0],  'filled')
    errorbar(1:length(resp_rate_new), resp_rate_new, SEM_resp_rate_new, SEM_resp_rate_new,  'color', [0 0 0], 'CapSize', 12)
    
    ylabel('Respiratory Rate (bpm)')
    set(gca,'XTick', [1:7], 'XTickLabel', {'BL','SL' '5','15','30','60','90'})
    xlim([0 8])
    ylim([25 55])
end


%% Bar charts to compare WT vs Chat-Cre at 30 min post J60 timepoint

% Left AUC
norm_AUC_ldia_wt_thirty = nanmean(norm_AUC_ldia_wt(4,:),2);
norm_AUC_ldia_chat_thirty = nanmean(norm_AUC_ldia_chat(4,:),2);
SEM_norm_AUC_ldia_wt_thirty = (nanstd(norm_AUC_ldia_wt(4,:),[],2)/sqrt(size(norm_AUC_ldia_wt(4,:),2)));
SEM_norm_AUC_ldia_chat_thirty = (nanstd(norm_AUC_ldia_chat(4,:),[],2)/sqrt(size(norm_AUC_ldia_chat(4,:),2)));

figure('Name', "30 mins post-J60 timepoint")
subplot(3,3,1)
hold on
bar(1,norm_AUC_ldia_wt_thirty,'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_AUC_ldia_wt,2)), norm_AUC_ldia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_AUC_ldia_wt_thirty,[],SEM_norm_AUC_ldia_wt_thirty,'k', 'CapSize', 12)
bar(2,norm_AUC_ldia_chat_thirty,'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_AUC_ldia_chat,2)), norm_AUC_ldia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_AUC_ldia_chat_thirty,[],SEM_norm_AUC_ldia_chat_thirty,'k', 'CapSize', 12)
ylabel('Left diaphragm EMG AUC (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')

% Left PK2PK
norm_PK2PK_ldia_wt_thirty = nanmean(norm_PK2PK_ldia_wt(4,:),2);
norm_PK2PK_ldia_chat_thirty = nanmean(norm_PK2PK_ldia_chat(4,:),2);
SEM_norm_PK2PK_ldia_wt_thirty = (nanstd(norm_PK2PK_ldia_wt(4,:),[],2)/sqrt(size(norm_PK2PK_ldia_wt(4,:),2)));
SEM_norm_PK2PK_ldia_chat_thirty = (nanstd(norm_PK2PK_ldia_chat(4,:),[],2)/sqrt(size(norm_PK2PK_ldia_chat(4,:),2)));


subplot(3,3,2)
hold on
bar(1,norm_PK2PK_ldia_wt_thirty,'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_PK2PK_ldia_wt,2)), norm_PK2PK_ldia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_PK2PK_ldia_wt_thirty,[],SEM_norm_PK2PK_ldia_wt_thirty,'k', 'CapSize', 12)
bar(2,norm_PK2PK_ldia_chat_thirty,'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_PK2PK_ldia_chat,2)), norm_PK2PK_ldia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_PK2PK_ldia_chat_thirty,[],SEM_norm_PK2PK_ldia_chat_thirty,'k', 'CapSize', 12)
ylabel('Left diaphragm PK2PK amplitude (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')

% Left Tonic
norm_tonic_ldia_wt_thirty = nanmean(norm_tonic_ldia_wt(4,:),2);
norm_tonic_ldia_chat_thirty = nanmean(norm_tonic_ldia_chat(4,:),2);
SEM_norm_tonic_ldia_wt_thirty = (nanstd(norm_tonic_ldia_wt(4,:),[],2)/sqrt(size(norm_tonic_ldia_wt(4,:),2)));
SEM_norm_tonic_ldia_chat_thirty = (nanstd(norm_tonic_ldia_chat(4,:),[],2)/sqrt(size(norm_tonic_ldia_chat(4,:),2)));


subplot(3,3,3)
hold on
bar(1,norm_tonic_ldia_wt_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_tonic_ldia_wt,2)), norm_tonic_ldia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_tonic_ldia_wt_thirty,[],SEM_norm_tonic_ldia_wt_thirty,'k', 'CapSize', 12)
bar(2,norm_tonic_ldia_chat_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_tonic_ldia_chat,2)), norm_tonic_ldia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_tonic_ldia_chat_thirty,[],SEM_norm_tonic_ldia_chat_thirty,'k', 'CapSize', 12)
ylabel('Left diaphragm tonic activity (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')
ylim([0 500])

%% Right Side - second row

% Right AUC
norm_AUC_rdia_wt_thirty = nanmean(norm_AUC_rdia_wt(4,:),2);
norm_AUC_rdia_chat_thirty = nanmean(norm_AUC_rdia_chat(4,:),2);
SEM_norm_AUC_rdia_wt_thirty = (nanstd(norm_AUC_rdia_wt(4,:),[],2)/sqrt(size(norm_AUC_rdia_wt(4,:),2)));
SEM_norm_AUC_rdia_chat_thirty = (nanstd(norm_AUC_rdia_chat(4,:),[],2)/sqrt(size(norm_AUC_rdia_chat(4,:),2)));

subplot(3,3,4)
hold on
bar(1,norm_AUC_rdia_wt_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_AUC_rdia_wt,2)), norm_AUC_rdia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_AUC_rdia_wt_thirty,[],SEM_norm_AUC_rdia_wt_thirty,'k','CapSize', 12)
bar(2,norm_AUC_rdia_chat_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_AUC_rdia_chat,2)), norm_AUC_rdia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_AUC_rdia_chat_thirty,[],SEM_norm_AUC_rdia_chat_thirty,'k', 'CapSize', 12)
ylabel('Right diaphragm EMG AUC (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')

% Right PK2PK
norm_PK2PK_rdia_wt_thirty = nanmean(norm_PK2PK_rdia_wt(4,:),2);
norm_PK2PK_rdia_chat_thirty = nanmean(norm_PK2PK_rdia_chat(4,:),2);
SEM_norm_PK2PK_rdia_wt_thirty = (nanstd(norm_PK2PK_rdia_wt(4,:),[],2)/sqrt(size(norm_PK2PK_rdia_wt(4,:),2)));
SEM_norm_PK2PK_rdia_chat_thirty = (nanstd(norm_PK2PK_rdia_chat(4,:),[],2)/sqrt(size(norm_PK2PK_rdia_chat(4,:),2)));


subplot(3,3,5)
hold on
bar(1,norm_PK2PK_rdia_wt_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_PK2PK_rdia_wt,2)), norm_PK2PK_rdia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_PK2PK_rdia_wt_thirty,[],SEM_norm_PK2PK_rdia_wt_thirty,'k', 'CapSize', 12)
bar(2,norm_PK2PK_rdia_chat_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_PK2PK_rdia_chat,2)), norm_PK2PK_rdia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_PK2PK_rdia_chat_thirty,[],SEM_norm_PK2PK_rdia_chat_thirty,'k', 'CapSize', 12)
ylabel('Right diaphragm PK2PK amplitude (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')

% Right Tonic
norm_tonic_rdia_wt_thirty = nanmean(norm_tonic_rdia_wt(4,:),2);
norm_tonic_rdia_chat_thirty = nanmean(norm_tonic_rdia_chat(4,:),2);
SEM_norm_tonic_rdia_wt_thirty = (nanstd(norm_tonic_rdia_wt(4,:),[],2)/sqrt(size(norm_tonic_rdia_wt(4,:),2)));
SEM_norm_tonic_rdia_chat_thirty = (nanstd(norm_tonic_rdia_chat(4,:),[],2)/sqrt(size(norm_tonic_rdia_chat(4,:),2)));


subplot(3,3,6)
hold on
bar(1,norm_tonic_rdia_wt_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(norm_tonic_rdia_wt,2)), norm_tonic_rdia_wt(4,:), 'k', 'jitter', 'on')
errorbar(1,norm_tonic_rdia_wt_thirty,[],SEM_norm_tonic_rdia_wt_thirty,'k', 'CapSize', 12)
bar(2,norm_tonic_rdia_chat_thirty, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(norm_tonic_rdia_chat,2)), norm_tonic_rdia_chat(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,norm_tonic_rdia_chat_thirty,[],SEM_norm_tonic_rdia_chat_thirty,'k', 'CapSize', 12)
ylabel('Right diaphragm tonic activity (% of BL)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})
plot([0 3], [100 100], 'k--')


%% Respiratory rate
resp_rate_wt = nanmean(resp_rate{2}(4,:),2);
SEM_resp_rate_wt = (nanstd(resp_rate{2}(4,:),[],2)/sqrt(size(resp_rate{2}(4,:),2)));

resp_rate_chat = nanmean(resp_rate{1}(4,:),2);
SEM_resp_rate_chat = (nanstd(resp_rate{1}(4,:),[],2)/sqrt(size(resp_rate{1}(4,:),2)));

subplot(3,3,8)
hold on

bar(1,resp_rate_wt, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(1,1,size(resp_rate{2},2)), resp_rate{2}(4,:), 'k', 'jitter', 'on')
errorbar(1,resp_rate_wt,[],SEM_resp_rate_wt,'k', 'CapSize', 12)
bar(2,resp_rate_chat, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0])
scatter(repmat(2,1,size(resp_rate{1},2)), resp_rate{1}(4,:), 'k', 'filled', 'jitter', 'on')
errorbar(2,resp_rate_chat,[],SEM_resp_rate_chat,'k', 'CapSize', 12)
ylabel('Respiratory Rate (bpm)')
set(gca, 'XTick', [1:2], 'XTickLabel', {'Wildtype', 'ChAT-Cre'})

