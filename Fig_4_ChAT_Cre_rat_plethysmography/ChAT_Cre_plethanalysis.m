clc
clear all
close all
tic

%%
DATAOUTPUT = [];
trace_plots = 'no';
summary_plots = 'yes';
root = '~\Fig_4_ChAT_Cre_rat_plethysmography\data'; % change directory is where data files are located
cd(root)
count = 0;

%%
animals = dir;
animals = animals(3:end);

for an = 1:numel(animals)
    
    mult_fact = 1.1;
    
    cd([root '\' animals(an).name '\mat_export'])
    files = dir;
    files = files(3:end);
    
    count = 0;
    for f = 1:numel(files)
        clear Ti_vec Te_vec
        count = count + 1;
        %% change to directory and load data
        fname = files(f).name;
        underbreaks = strfind(fname, '_');
        
        load(fname)
        
        
        if strcmp(fname(underbreaks(end)+1:end-4),'salineBL')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = 'BL';
            fig_base = 10;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'salineinfusion')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = 'infusion';
            fig_base = 20;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'saline05')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = '0-5 mins';
            fig_base = 30;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'saline1015')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = '10-15 mins';
            fig_base = 40;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'saline2530')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = '25-30 mins';
            fig_base = 50;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'salinemaxchemo')
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = 'maxchemo';
            fig_base = 60;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j60BL')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = 'BL';
            fig_base = 70;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j60infusion')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = 'infusion';
            fig_base = 80;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j6005')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = '0-5 mins';
            fig_base = 90;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j601015')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = '10-15 mins';
            fig_base = 100;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j602530')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = '25-30 mins';
            fig_base = 110;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j60maxchemo')
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = 'maxchemo';
            fig_base = 120;
        else
            disp('problem')
            keyboard
        end
        
        %% calculate basic parameters
        sr = 1/Vt.interval;
        ms = sr/1000;
        [b,a] = butter(2,[0.5 10]/(sr/2));
        filteredvt = filtfilt(b,a,NormVt.values);
        
        samples_to_examine = 1:numel(filteredvt);
        
        
        % calculate suggested MPPs on oversmoothed trace
        [N,X] = hist(filteredvt(samples_to_examine)*-1,1000);
        N = smooth(N,100);
        [~,LOCS] = findpeaks(N,'MinPeakProminence',2);
        sug_MPP = abs(X(LOCS(end))*mult_fact);
        
        % hand determine MPP for subset of traces where autodetecting did a
        % bad job
        if strcmp(fname, '02_24_21AP_20230317_salineinfusion.mat') %good
            sug_MPP = 0.26;
        elseif strcmp(fname, '03_25_21AP_20230404_saline05.mat') %good
            sug_MPP = 0.28;
        elseif strcmp(fname, '10_25_21AP_20230404_j602530.mat') %good
            sug_MPP = 0.30;
        elseif strcmp(fname, '08_29_21AP_20230711_j6005.mat') %good
            sug_MPP = 0.19;
        elseif strcmp(fname, '05_30_21AP_20230710_saline2530.mat') %good
            sug_MPP = 0.25;
        elseif strcmp(fname, '06_30_21AP_20230710_j60BL.mat') %good
            sug_MPP = 0.26;
        elseif strcmp(fname, '07_30_21AP_20230710_j60infusion.mat') %good
            sug_MPP = 0.24;
        elseif strcmp(fname, '09_30_21AP_20230710_j601015.mat') %good
            sug_MPP = 0.25;
        elseif strcmp(fname, '10_30_21AP_20230710_j602530.mat') %good
            sug_MPP = 0.2;
        elseif strcmp(fname, '04_31_21AP_20230707_saline1015.mat') %good
            sug_MPP = 0.12;
        elseif strcmp(fname, '05_32_21AP_20230712_saline2530.mat') %good
            sug_MPP = 0.2;
        elseif strcmp(fname, '08_33_21AP_20230802_j6005.mat') %good
            sug_MPP = 0.25;
        elseif strcmp(fname, '09_34_21AP_20230706_j601015.mat')
            sug_MPP = 0.2;
        else
        end
        
        if strcmp(trace_plots, 'yes')
            %find peaks and locs on oversmoothed trace
            figure(fig_base+an)
            ax = get(gcf,'children');
            subplot(2,1,numel(ax)+1)
            findpeaks(filteredvt(samples_to_examine)*-1,sr,'minPeakProminence',sug_MPP)
            title(fname,'interpreter','none')
            ylabel('NormVT');
            linkaxes
            set(gcf,'pos',[ 1          41        1920         963])
        elseif strcmp(trace_plots, 'no')
        else
        end
        
        [I_PKS, I_LOCS] = findpeaks(filteredvt(samples_to_examine)*-1, sr,'minPeakProminence',sug_MPP);
        [E_PKS, E_LOCS] = findpeaks(filteredvt(samples_to_examine), sr, 'minPeakProminence',sug_MPP);
        
        
        duration = I_LOCS(end)-I_LOCS(1); %duration in seconds
        rr(f,an) = numel(I_PKS)/(duration/60); %breaths/min
        
        
        % instantaneous resp rate
        %calculate instantaneous resp rate
        time_RR_diff = diff(I_LOCS);
        resp_rate{f,an} = 60./time_RR_diff;
        
        if numel(I_PKS) > numel(E_PKS)
            I_PKS = I_PKS(1:numel(E_PKS));
        elseif numel(E_PKS) > numel(I_PKS)
            E_PKS = E_PKS(1:numel(I_PKS));
            
        else
        end
        Vt_vec = I_PKS + E_PKS;
        
        normtidalvolume(f,an) = mean(Vt_vec);
        normve(f,an) = normtidalvolume(f,an)*rr(f,an);
        
        %% calculate VT
        
        filteredrealvt = filtfilt(b,a,Vt.values);
        
        vtsamples_to_examine = 1:numel(filteredrealvt);
        
        [N,X] = hist(filteredrealvt(vtsamples_to_examine)*-1,1000);
        N = smooth(N,100);
        [~,LOCS] = findpeaks(N,'MinPeakProminence',2);
        sug_MPPvt = abs(X(LOCS(end))*mult_fact);
        
       
        % hand determine MPP for subset of traces where autodetecting did a
        % bad job
        if strcmp(fname, '02_24_21AP_20230317_salineinfusion.mat') %good
            sug_MPPvt = 0.15;
        elseif strcmp(fname, '03_25_21AP_20230404_saline05.mat') %good
            sug_MPPvt = 0.1;
        elseif strcmp(fname, '10_25_21AP_20230404_j602530.mat') %good
            sug_MPPvt = 0.1;
        elseif strcmp(fname, '08_29_21AP_20230711_j6005.mat') %good
            sug_MPPvt = 0.14;
        elseif strcmp(fname, '05_30_21AP_20230710_saline2530.mat') %good
            sug_MPPvt = 0.1;
        elseif strcmp(fname, '05_30_21AP_20230710_j60BL.mat') %good
            sug_MPPvt = 0.13;
        elseif strcmp(fname, '06_30_21AP_20230710_j60BL.mat') %good
            sug_MPPvt = 0.08;
        elseif strcmp(fname, '07_30_21AP_20230710_j60infusion.mat') %good
            sug_MPPvt = 0.1;
        elseif strcmp(fname, '09_30_21AP_20230710_j601015.mat') %good
            sug_MPPvt = 0.09;
        elseif strcmp(fname, '10_30_21AP_20230710_j602530.mat') %good
            sug_MPPvt = 0.09;
        elseif strcmp(fname, '04_31_21AP_20230707_saline1015.mat') %good
            sug_MPPvt = 0.08;
        elseif strcmp(fname, '05_32_21AP_20230712_saline2530.mat') %good
            sug_MPPvt = 0.08;
        elseif strcmp(fname, '08_33_21AP_20230802_j6005.mat') %good
            sug_MPPvt = 0.1;
        elseif strcmp(fname, '09_34_21AP_20230706_j601015.mat')
            sug_MPPvt = 0.08;
            
        else
        end
        
        if strcmp(trace_plots, 'yes')
            figure(fig_base+an)
            ax = get(gcf,'children');
            subplot(2,1,numel(ax)+1)
            findpeaks(filteredrealvt(vtsamples_to_examine)*-1,sr,'minPeakProminence',sug_MPPvt)
            title(fname,'interpreter','none')
            ylabel('VT');
            linkaxes
            set(gcf,'pos',[ 1          41        1920         963])
        elseif strcmp(trace_plots, 'no')
        else
        end
        
        [Ir_PKS, Ir_LOCS] = findpeaks(filteredrealvt(vtsamples_to_examine)*-1,sr,'minPeakProminence',sug_MPPvt);
        [Er_PKS, Er_LOCS] = findpeaks(filteredrealvt(vtsamples_to_examine),sr,'minPeakProminence',sug_MPPvt);
        
        
        duration = Ir_LOCS(end)-Ir_LOCS(1); %duration in seconds

        if numel(Ir_PKS) > numel(Er_PKS)
            Ir_PKS = Ir_PKS(1:numel(Er_PKS));
        elseif numel(Er_PKS) > numel(Ir_PKS)
            Er_PKS = Er_PKS(1:numel(Ir_PKS));
            
        else
        end
        
        Vtr_vec = Ir_PKS + Er_PKS;
        VTR_vec_cell{f,an} = Vtr_vec;
        
        tidalvolume = mean(Vtr_vec);
        ve = tidalvolume*rr(f,an);
        
        %% determine Ti and Te
        if Er_LOCS(1) < Ir_LOCS(1)
            Er_LOCS = Er_LOCS(2:end);
        else
        end
        
        if numel(Ir_LOCS) > numel(Er_LOCS)
            Ir_LOCS = Ir_LOCS(1:numel(Er_LOCS));
        elseif numel(Er_LOCS) > numel(Ir_LOCS)
            Er_LOCS = Er_LOCS(1:numel(Ir_LOCS));
        else
        end
        
        for y = 1:numel(Er_LOCS)-1
            Ti_vec(y) = Er_LOCS(y) - Ir_LOCS(y);
            Te_vec(y) = Ir_LOCS(y+1) - Er_LOCS(y);
        end
        
        avg_Ti = nanmean(Ti_vec);
        avg_Te = nanmean(Te_vec);
        
        %% determine metabolism
        if strcmp(segment,'Maxchemo')
            metab = NaN;
            vevco2 = NaN;
        elseif strcmp(segment,'Hypercapnia')
            metab = NaN;
            vevco2 = NaN;
        else
            metab = mean(VCO2.values);
            vevco2 = normve(f,an)/metab;
        end
        
        %% data output
        data = num2cell([rr(f,an), normtidalvolume(f,an), normve(f,an), vevco2, metab, tidalvolume, ve, avg_Ti, avg_Te]);
        data = [anid, fname, segment, timepoint, data];
        
        DATAOUTPUT = [DATAOUTPUT;data];
    end
end

%% make variables to make plotting easier

%raw vars
normtv_saline = normtidalvolume(1:5,:);
normtv_j60 = normtidalvolume(6:10,:);

avg_raw_normtv_saline = mean(normtv_saline,2);
avg_raw_normtv_j60 = mean(normtv_j60,2);

raw_normtv_saline_SEM = (std(normtv_saline, [], 2)./sqrt(numel(animals)));
raw_normtv_j60_SEM = (std(normtv_j60, [], 2)./sqrt(numel(animals)));

rr_saline = rr(1:5,:);
rr_j60 = rr(6:10,:);

avg_raw_rr_saline = mean(rr_saline,2);
avg_raw_rr_j60 = mean(rr_j60,2);

raw_rr_saline_SEM = (std(rr_saline, [], 2)./sqrt(numel(animals)));
raw_rr_j60_SEM = (std(rr_j60, [], 2)./sqrt(numel(animals)));

ve_saline = normve(1:5,:);
ve_j60 = normve(6:10,:);

avg_raw_ve_saline = mean(ve_saline,2);
avg_raw_ve_j60 = mean(ve_j60,2);

raw_ve_saline_SEM = (std(ve_saline, [], 2)./sqrt(numel(animals)));
raw_ve_j60_SEM = (std(ve_j60, [], 2)./sqrt(numel(animals)));

%normalized vars
norm_fact_tv_s = normtv_saline(1,:);
norm_normtv_saline = (normtv_saline./norm_fact_tv_s)*100;

norm_fact_tv_j60 = normtv_j60(1,:);
norm_normtv_j60 = (normtv_j60./norm_fact_tv_j60)*100;

avg_normtv_saline = mean(norm_normtv_saline,2);
avg_normtv_j60 = mean(norm_normtv_j60,2);

normtv_saline_SEM = (std(norm_normtv_saline, [], 2)./sqrt(numel(animals)));
normtv_j60_SEM = (std(norm_normtv_j60, [], 2)./sqrt(numel(animals)));

rr_norm_fact_saline = rr_saline(1,:);
norm_rr_saline = (rr_saline./rr_norm_fact_saline)*100;

rr_norm_fact_j60 = rr_j60(1,:);
norm_rr_j60 = (rr_j60./rr_norm_fact_j60)*100;

avg_rr_saline = mean(norm_rr_saline,2);
avg_rr_j60 = mean(norm_rr_j60,2);

rr_saline_SEM = (std(norm_rr_saline, [], 2)./sqrt(numel(animals)));
rr_j60_SEM = (std(norm_rr_j60, [], 2)./sqrt(numel(animals)));

ve_norm_fact_saline = ve_saline(1,:);
norm_ve_saline = (ve_saline./ve_norm_fact_saline)*100;

ve_norm_fact_j60 = ve_j60(1,:);
norm_ve_j60 = (ve_j60./ve_norm_fact_j60)*100;

avg_ve_saline = mean(norm_ve_saline,2);
avg_ve_j60 = mean(norm_ve_j60,2);

ve_saline_SEM = (std(norm_ve_saline, [], 2)./sqrt(numel(animals)));
ve_j60_SEM = (std(norm_ve_j60, [], 2)./sqrt(numel(animals)));



if strcmp(summary_plots, 'yes')
    %% figure for publication plot
    figure
    subplot(2,3,1)
    hold on
    plot(avg_raw_normtv_saline,'b')
    errorbar(1:5, avg_raw_normtv_saline, raw_normtv_saline_SEM, raw_normtv_saline_SEM)
    scatter(1:5, avg_raw_normtv_saline, 'filled', 'b')
    plot(avg_raw_normtv_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_raw_normtv_j60, raw_normtv_j60_SEM, raw_normtv_j60_SEM)
    scatter(1:5, avg_raw_normtv_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), normtv_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), normtv_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Tidal volume (ml/kg)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    %% plot raw respiratory rate
    subplot(2,3,2)
    hold on
    plot(avg_raw_rr_saline,'b')
    errorbar(1:5, avg_raw_rr_saline, raw_rr_saline_SEM, raw_rr_saline_SEM)
    scatter(1:5, avg_raw_rr_saline, 'filled', 'b')
    plot(avg_raw_rr_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_raw_rr_j60, raw_rr_j60_SEM, raw_rr_j60_SEM)
    scatter(1:5, avg_raw_rr_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), rr_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), rr_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Respiratory rate (bpm)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    %% plot raw minute ventilation
    subplot(2,3,3)
    hold on
    plot(avg_raw_ve_saline,'b')
    errorbar(1:5, avg_raw_ve_saline, raw_ve_saline_SEM, raw_ve_saline_SEM)
    scatter(1:5, avg_raw_ve_saline, 'filled', 'b')
    plot(avg_raw_ve_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_raw_ve_j60, raw_ve_j60_SEM, raw_ve_j60_SEM)
    scatter(1:5, avg_raw_ve_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), ve_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), ve_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Minute ventilation (ml/min)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    
    
    
    %% plot normalized tidal volume
    subplot(2,3,4)
    hold on
    plot(avg_normtv_saline,'b')
    errorbar(1:5, avg_normtv_saline, normtv_saline_SEM, normtv_saline_SEM)
    scatter(1:5, avg_normtv_saline, 'filled', 'b')
    plot(avg_normtv_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_normtv_j60, normtv_j60_SEM, normtv_j60_SEM)
    scatter(1:5, avg_normtv_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), norm_normtv_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), norm_normtv_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Normalized tidal volume (% baseline)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    %% plot normalized respiratory rate
    subplot(2,3,5)
    hold on
    plot(avg_rr_saline,'b')
    errorbar(1:5, avg_rr_saline, rr_saline_SEM, rr_saline_SEM)
    scatter(1:5, avg_rr_saline, 'filled', 'b')
    plot(avg_rr_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_rr_j60, rr_j60_SEM, rr_j60_SEM)
    scatter(1:5, avg_rr_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), norm_rr_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), norm_rr_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Normalized Respiratory rate (% baseline)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    %% plot normalized minute ventilation
    subplot(2,3,6)
    hold on
    plot(avg_ve_saline,'b')
    errorbar(1:5, avg_ve_saline, ve_saline_SEM, ve_saline_SEM)
    scatter(1:5, avg_ve_saline, 'filled', 'b')
    plot(avg_ve_j60,'Color', [0.8500 0.3250 0.0980])
    errorbar(1:5, avg_ve_j60, ve_j60_SEM, ve_j60_SEM)
    scatter(1:5, avg_ve_j60, [], [0.8500 0.3250 0.0980], 'filled')
    for n = 1:5
        scatter(repmat(n, 1, numel(animals)), norm_ve_saline(n,:),[], 'b', 'filled', 'jitter', 'on')
        scatter(repmat(n, 1, numel(animals)), norm_ve_j60(n,:),[], [0.8500 0.3250 0.0980], 'filled', 'jitter', 'on')
    end
    xlim([0.75 5.25])
    ylabel('Normalized Minute ventilation (% baseline)')
    set(gca, 'XTick', [1:5], 'XTicklabel', {"Baseline", "Infusion", "0-5 mins", "10-15 mins", "25-30 mins"})
    
    
elseif strcmp(summary_plots, 'no')
    
else
end










