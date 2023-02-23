%% Individual Time-lag

%% Data Preparation
clear
clc
% load data
load 'DATA(3conditions)_rev210330.mat'

% Data Settings
EEG    = DATA(3).EEG;
SPEECH = DATA(3).SPEECH;
INDEX  = DATA(3).INDEX;

%% Individual Time-lag Analysis (only positive lags)

% ======================= Parameter Settings =========================== %
% Decoder Parameters 
Dir    = -1; % 1 for forward, -1 for backward
fs     = 64;
dur    = 60;
lambda = 10;

% Time-lag (tau) params
tmin    = 0;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 500;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );

% Experimental Params
Nchans  = length(EEG(1).chanlocs); % Number of Channels
Ntrials = size(EEG(1).data,1); % half of trials for each direction

% Init mtx for results
Accs   = zeros(length(EEG), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AAD = zeros(length(EEG), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AUD = zeros(length(EEG), length(t), 2); % NSub X NLag X Att. or Unatt.


% Loop over Subjects
for SubIdx = 1:length(EEG)

        % initialize SampledTrials cell array
        eeg_cell        = cell( Ntrials, 1 ); % cell for eeg data
        attended_cell   = cell( Ntrials, 1 ); % cell for attended speech
        unattended_cell = cell( Ntrials, 1 ); % cell for unattended speech
        
        % loop over SampledTrials
        for i = 1:Ntrials
            % TRIALS
            eeg_cell{i}        = squeeze(EEG(SubIdx).data(i, :, :))';
            attended_cell{i}   = SPEECH(INDEX(SubIdx).a(i),:)';
            unattended_cell{i} = SPEECH(INDEX(SubIdx).u(i),:)';
        end % The end of Ntrials loop
        
        
        % individual time-lags
        for ti = 1:length(t)
            % Attention Decoding
            [stats,stats1,stats2,~] = mTRFattncrossval(attended_cell, unattended_cell, eeg_cell, ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);
            
            Accs(SubIdx,ti,1)   = stats.acc;
            Rs_AAD(SubIdx,ti,1) = mean(stats1.r);
            Rs_AAD(SubIdx,ti,2) = mean(stats2.r);
            Sub_AAD_w{SubIdx, ti} = mean(stats.w, 2);
            
            
            % Unattention Decoding
            [stats,stats1,stats2,~] = mTRFattncrossval(unattended_cell, attended_cell, eeg_cell, ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);
            
            Accs(SubIdx,ti,2)   = stats.acc;
            Rs_AUD(SubIdx,ti,1) = mean(stats1.r);
            Rs_AUD(SubIdx,ti,2) = mean(stats2.r);
            Sub_AUD_w{SubIdx, ti} = mean(stats.w, 2);
            
        end % end time-lag
        
        disp(['Subject' num2str(SubIdx) ' finished!'])
        
end % end of Subject loop

% decoder weight array
Ave_AAD_w = zeros( Nchans+1, ti);
Ave_AUD_w = zeros( Nchans+1, ti);
for i = 1:ti
    
    Ave_AAD_w(:,i) = mean(cat(2, Sub_AAD_w{:,i}), 2);
    Ave_AUD_w(:,i) = mean(cat(2, Sub_AUD_w{:,i}), 2);
    
end

%% Plotting

% Data Curation
Ave_Accs   = squeeze(mean(Accs,1));
Ave_Rs_AAD = squeeze(mean(Rs_AAD,1));
Ave_Rs_AUD = squeeze(mean(Rs_AUD,1));

% figure settings
figure(3)
clf
set(gcf, 'color', 'w')

% plot
subplot(1,2,1) % AAD
plot(t,Ave_Accs(:,1)*100,'-o','LineWidth',5,'MarkerSize',20)
xticks(round(t));
set(gca, 'FontSize', 15)
%ylim([30,100])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Accuracy (%) ', 'FontSize', 25, 'FontWeight', 'bold')
title('AADC: AAD Performance', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,2,2) % AUD
plot(t,Ave_Accs(:,2)*100,'-*','LineWidth',5,'MarkerSize',20,'Color',[0.8500 0.3250 0.0980])
xticks(round(t));
%ylim([30,100])
set(gca, 'FontSize', 15)
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Accuracy (%) ', 'FontSize', 25, 'FontWeight', 'bold')
title('AADC: AUD Performance', 'FontSize', 30, 'FontWeight', 'bold')


% figure settings
figure(4)
clf
set(gcf, 'color', 'w')

% plot
subplot(1,2,1) % AAD
plot(t,Ave_Rs_AAD(:,1),'-o','LineWidth',5,'MarkerSize',15)
hold on
plot(t,Ave_Rs_AAD(:,2),'-*','LineWidth',5,'MarkerSize',15)
xticks(round(t)); yticks(-0.1:0.01:0.1);
set(gca, 'FontSize', 15, 'ylim', [-0.1, 0.1])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_att', 'r\_unatt'})
title('AAK: Auditory Attention Decoder', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,2,2) % AUD
plot(t,Ave_Rs_AUD(:,1),'-*','LineWidth',5,'MarkerSize',15,'Color',[0.8500 0.3250 0.0980])
hold on
plot(t,Ave_Rs_AUD(:,2),'-o','LineWidth',5,'MarkerSize',15,'Color',[0 0.4470 0.7410])
xticks(round(t)); yticks(-0.05:0.005:0.05);
set(gca, 'FontSize', 15, 'ylim', [-0.05, 0.05])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_unatt', 'r\_att'})
title('AAK: Auditory Unattention Decoder', 'FontSize', 30, 'FontWeight', 'bold')

%% Topography of the decoder weghts averaged over all subjects
% OPENBCI cap
% load channel location
chanlocs = readlocs('Standard-10-20-Cap19.locs');
chanlocs([1 2 18 19]) = [];
% relocation channel location
Chanloc_ob = [chanlocs(3), chanlocs(8), chanlocs(7), chanlocs(9), chanlocs(11), chanlocs(15), chanlocs(13), chanlocs(1), chanlocs(5), ...
                    chanlocs(2), chanlocs(4), chanlocs(6), chanlocs(10), chanlocs(12), chanlocs(14)];
              
for i = 1:ti
    figure(i)
    set(gcf, 'color', 'w')
    subplot(121)
    topoplot(Ave_AAD_w(3:end,i),  Chanloc_ob);
    subplot(122)
    topoplot(Ave_AUD_w(3:end,i),  Chanloc_ob);
    
    sgtitle(t(i));
    saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\IndiTopoG_'+t(i)+'.jpg');

end

