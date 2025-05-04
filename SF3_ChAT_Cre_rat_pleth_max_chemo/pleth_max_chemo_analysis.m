clear all
close all
clc

DATAOUTPUT = [];
trace_plots = 'no';
root = '~\SF3_ChAT_Cre_rat_pleth_max_chemo\data'; % change to directory where data files are located
cd(root)
count = 0;

%%
animals = dir;
animals = animals(3:end);

for an = 1:numel(animals)
    
    mult_fact = 1.1;
    
    
    
    cd([root '\' animals(an).name])
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
            % dur = 5*60;
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = 'BL';
            fig_base = 10;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'salinemaxchemo')
            %dur = 3*60;
            anid = animals(an).name(1:2);
            segment = 'Saline';
            timepoint = 'maxchemo';
            fig_base = 60;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j60BL')
            % dur = 5*60;n
            anid = animals(an).name(1:2);
            segment = 'J60';
            timepoint = 'BL';
            fig_base = 70;
        elseif strcmp(fname(underbreaks(end)+1:end-4),'j60maxchemo')
            %dur = 3*60;
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
        
        if strcmp(fname, '06_30_21AP_20230710_j60BL.mat') %good
            sug_MPP = 0.26;
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
        
        if strcmp(fname, '05_30_21AP_20230710_j60BL.mat') %good
            sug_MPPvt = 0.13;
        elseif strcmp(fname, '06_30_21AP_20230710_j60BL.mat') %good
            sug_MPPvt = 0.08;
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
        %
        
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
        

        %% data output
        data = num2cell([rr(f,an), normtidalvolume(f,an), normve(f,an)]);
        data = [anid, fname, segment, timepoint, data];
        
        DATAOUTPUT = [DATAOUTPUT;data];
    end
end

%% norm vars
norm_tv = (normtidalvolume([2 4],:)./normtidalvolume([1 3],:))*100;
norm_rr = (rr([2 4],:)./rr([1 3],:))*100;
norm_ve = (normve([2 4],:)./normve([1 3],:))*100;

%% Plotting
    %tidal volume
    figure
    subplot(2,3,1)
    hold on
    bar(1, mean(normtidalvolume(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(normtidalvolume(2,:),2),1), normtidalvolume(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(normtidalvolume(2,:)), [], (std(normtidalvolume(2,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(normtidalvolume(4,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(normtidalvolume(4,:),2),1), normtidalvolume(4,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(normtidalvolume(4,:)), [], (std(normtidalvolume(4,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Tidal Volume (ml/kg)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    %resp rate
    subplot(2,3,2)
    hold on
    bar(1, mean(rr(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(rr(2,:),2),1), rr(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(rr(2,:)), [], (std(rr(2,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(rr(4,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(rr(4,:),2),1), rr(4,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(rr(4,:)), [], (std(rr(4,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Respiratory Rate (bpm)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    %min ventilation
    subplot(2,3,3)
    hold on
    bar(1, mean(normve(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(normve(2,:),2),1), normve(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(normve(2,:)), [], (std(normve(2,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(normve(4,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(normve(4,:),2),1), normve(4,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(normve(4,:)), [], (std(normve(4,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Minute Ventilation (ml/min)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    
    % norm vals
        subplot(2,3,4)
    hold on
    bar(1, mean(norm_tv(1,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(norm_tv(1,:),2),1), norm_tv(1,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(norm_tv(1,:)), [], (std(norm_tv(1,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(norm_tv(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(norm_tv(2,:),2),1), norm_tv(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(norm_tv(2,:)), [], (std(norm_tv(2,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Tidal Volume (% BL)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    %resp rate
    subplot(2,3,5)
    hold on
    bar(1, mean(norm_rr(1,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(norm_rr(1,:),2),1), norm_rr(1,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(norm_rr(1,:)), [], (std(norm_rr(1,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(norm_rr(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(norm_rr(2,:),2),1), norm_rr(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(norm_rr(2,:)), [], (std(norm_rr(2,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Respiratory Rate (% BL)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    %min ventilation
    subplot(2,3,6)
    hold on
    bar(1, mean(norm_ve(1,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(1, size(norm_ve(1,:),2),1), norm_ve(1,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(1, mean(norm_ve(1,:)), [], (std(norm_ve(1,:))./numel(animals)), 'k', 'Capsize', 15)
    bar(2,mean(norm_ve(2,:)),'EdgeColor', [0 0 0], 'FaceColor', [1 1 1])
    scatter(repmat(2, size(norm_ve(2,:),2),1), norm_ve(2,:), [], 'k', 'filled', 'jitter', 'on')
    errorbar(2, mean(norm_ve(2,:)), [], (std(norm_ve(2,:))./numel(animals)), 'k', 'Capsize', 15)
    ylabel('Hypercapnic Hypoxia - Minute Ventilation (% BL)')
    set(gca, 'XTick', [1:2], 'XTicklabel', {"Saline", "J60"})
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    