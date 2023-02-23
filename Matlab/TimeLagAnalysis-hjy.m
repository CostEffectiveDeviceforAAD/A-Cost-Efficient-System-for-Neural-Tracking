%% Individual Time-lag

%% Data Preparation

% load data
load 'DATA(3conditions)_rev210330.mat'

% Data Settings
EEG    = DATA(3).EEG;
SPEECH = DATA(3).SPEECH;
INDEX  = DATA(3).INDEX;
clear DATA

%% Data - jiyeon

load 'AllEEG_seg.mat'   % 1 by sub cell=DATA
load 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\Allspeech.mat'   % 60by3840 -Allspeech


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
tmax    = 250;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );

% Experimental Params
Nchans  = 15; % Number of Channels
Ntrials = 30; % half of trials for each direction

% Init mtx for results
Accs   = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AAD = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AUD = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.


% Loop over Subjects
for SubIdx = 1:length(DATA)

        % initialize SampledTrials cell array
        eeg_cell        = cell( Ntrials, 1 ); % cell for eeg data
        attended_cell   = cell( Ntrials, 1 ); % cell for attended speech
        unattended_cell = cell( Ntrials, 1 ); % cell for unattended speech
        
        % loop over SampledTrials
        for i = 1:Ntrials
            % TRIALS
            eeg_cell{i}        = DATA{SubIdx}(:,:,i);
            attended_cell{i}   = Allspeech(30+i,:)';
            unattended_cell{i} = Allspeech(i,:)';
        end % The end of Ntrials loop
        
        
        % individual time-lags
        for ti = 1:length(t)
            % Attention Decoding
            [stats,stats1,stats2,~] = mTRFattncrossval(attended_cell, unattended_cell, eeg_cell, ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);
            % train 만
            [stats_trn,stats1_trn,stats2_trn,~] = mTRFattncrossval(attended_cell(1:14), unattended_cell(1:14), eeg_cell(1:14), ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);
            
            Accs(SubIdx,ti,1)   = stats.acc;
            Accs_trn(SubIdx,ti,1)   = stats_trn.acc;
            Rs_AAD(SubIdx,ti,1) = mean(stats1.r);
            Rs_AAD(SubIdx,ti,2) = mean(stats2.r);
            Sub_AAD_w{SubIdx, ti} = mean(stats.w, 2);
            
            
            % Unattention Decoding
            [stats,stats1,stats2,~] = mTRFattncrossval(unattended_cell, attended_cell, eeg_cell, ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);
            [stats_trn,stats1_trn,stats2_trn,~] = mTRFattncrossval(unattended_cell(1:14), attended_cell(1:14), eeg_cell(1:14), ...
                fs, Dir, t(ti), t(ti), lambda, 'verbose', 0);   
            
            Accs(SubIdx,ti,2)   = stats.acc;
            Accs_trn(SubIdx,ti,2)   = stats_trn.acc;
            Rs_AUD(SubIdx,ti,1) = mean(stats1.r);
            Rs_AUD(SubIdx,ti,2) = mean(stats2.r);
            Sub_AUD_w{SubIdx, ti} = mean(stats.w, 2); 
            
        end % end time-lag
        
        disp(['Subject' num2str(SubIdx) ' finished!'])
        
        
end % end of Subject loop

Ave_AAD_w = zeros( Nchans+1, ti);
Ave_AUD_w = zeros( Nchans+1, ti);

for i = 1:ti
    
    Ave_AAD_w(:,i) = mean(cat(2, Sub_AAD_w{:,i}), 2);
    Ave_AUD_w(:,i) = mean(cat(2, Sub_AUD_w{:,i}), 2);
    
end


%% Plotting

% Data Curation
Ave_Accs   = squeeze(mean(Accs,1));
Ave_Accs_trn   = squeeze(mean(Accs_trn,1));
Ave_Rs_AAD = squeeze(mean(Rs_AAD,1));
Ave_Rs_AUD = squeeze(mean(Rs_AUD,1));

% figure settings
figure(1)
clf
set(gcf, 'color', 'w')

% plot
subplot(1,2,1) % AAD
plot(t,Ave_Accs(:,1)*100,'-o','LineWidth',5,'MarkerSize',20)
xticks(round(t));
%set(gca, 'FontSize', 15, [40, 90])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Accuracy (%) ', 'FontSize', 25, 'FontWeight', 'bold')
title('AAD Performance', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,2,2) % AUD
plot(t,Ave_Accs(:,2)*100,'-*','LineWidth',5,'MarkerSize',20,'Color',[0.8500 0.3250 0.0980])
xticks(round(t));
set(gca, 'FontSize', 15, 'ylim', [40, 90])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Accuracy (%) ', 'FontSize', 25, 'FontWeight', 'bold')
title('AUD Performance', 'FontSize', 30, 'FontWeight', 'bold')


% figure settings
figure(2)
clf
set(gcf, 'color', 'w')

% plot
subplot(1,2,1) % AAD
plot(t,Ave_Rs_AAD(:,1),'-o','LineWidth',5,'MarkerSize',15)
hold on
plot(t,Ave_Rs_AAD(:,2),'-*','LineWidth',5,'MarkerSize',15)
xticks(round(t)); yticks(0:0.01:0.1);
set(gca, 'FontSize', 15, 'ylim', [0, 0.1])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_att', 'r\_unatt'})
title('Auditory Attention Decoder', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,2,2) % AUD
plot(t,Ave_Rs_AUD(:,1),'-*','LineWidth',5,'MarkerSize',15,'Color',[0.8500 0.3250 0.0980])
hold on
plot(t,Ave_Rs_AUD(:,2),'-o','LineWidth',5,'MarkerSize',15,'Color',[0 0.4470 0.7410])
xticks(round(t)); yticks(0:0.005:0.05);
set(gca, 'FontSize', 15, 'ylim', [0, 0.05])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_unatt', 'r\_att'})
title('Auditory Unattention Decoder', 'FontSize', 30, 'FontWeight', 'bold')

%% Topography of the decoder weghts averaged over all subjects

chanlocs = readlocs('Standard-10-20-Cap19.locs');
chanlocs([1 2 18 19]) = [];
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
end

%% Individual Time-lag Analysis (including negative lags)

% ======================= Parameter Settings =========================== %
% Decoder Parameters 
Dir    = -1; % 1 for forward, -1 for backward
fs     = 64;
dur    = 60;
lambda = 10;

% Time-lag (tau) params
tmin    = -250;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 250;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );

% Experimental Params
% Experimental Params
Nchans  = 15; % Number of Channels
Ntrials = 30; % half of trials for each direction

% Init mtx for results
Accs   = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AAD = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.
Rs_AUD = zeros(length(DATA), length(t), 2); % NSub X NLag X Att. or Unatt.


% Loop over Subjects
for SubIdx = 1:length(DATA)

        % initialize SampledTrials cell array
        eeg_cell        = cell( Ntrials, 1 ); % cell for eeg data
        attended_cell   = cell( Ntrials, 1 ); % cell for attended speech
        unattended_cell = cell( Ntrials, 1 ); % cell for unattended speech
        
        % loop over SampledTrials
        for i = 1:Ntrials
            % TRIALS
            eeg_cell{i}        = DATA{SubIdx}(:,:,i);
            attended_cell{i}   = Allspeech(30+i,:)';
            unattended_cell{i} = Allspeech(i,:)';
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
figure(1)
clf
set(gcf, 'color', 'w')

% plot
plot(t,Ave_Accs(:,1)*100,'-o','LineWidth',5,'MarkerSize',20)
hold on
plot(t,Ave_Accs(:,2)*100,'-*','LineWidth',5,'MarkerSize',20,'Color',[0.8500 0.3250 0.0980])
xticks(round(t));
legend({'AAD', 'AUD'})
set(gca, 'FontSize', 15, 'ylim', [40, 100])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Accuracy (%) ', 'FontSize', 25, 'FontWeight', 'bold')
title('AAD & AUD Performance', 'FontSize', 30, 'FontWeight', 'bold')


% figure settings
figure(2)
clf
set(gcf, 'color', 'w')

% plot
subplot(2,1,1) % AAD
plot(t,Ave_Rs_AAD(:,1),'-o','LineWidth',5,'MarkerSize',15)
hold on
plot(t,Ave_Rs_AAD(:,2),'-*','LineWidth',5,'MarkerSize',15)
xticks(round(t)); yticks(-0.01:0.01:0.07);
set(gca, 'FontSize', 15, 'ylim', [-0.01, 0.07])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_att', 'r\_unatt'})
title('Auditory Attention Decoder', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,1,2) % AUD
plot(t,Ave_Rs_AUD(:,1),'-*','LineWidth',5,'MarkerSize',15,'Color',[0.8500 0.3250 0.0980])
hold on
plot(t,Ave_Rs_AUD(:,2),'-o','LineWidth',5,'MarkerSize',15,'Color',[0 0.4470 0.7410])
xticks(round(t)); yticks(-0.01:0.01:0.05);
set(gca, 'FontSize', 15, 'ylim', [-0.01, 0.05])
xlabel(' Time-lag (ms) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Correlation ', 'FontSize', 25, 'FontWeight', 'bold')
legend({'r\_unatt', 'r\_att'})
title('Auditory Unattention Decoder', 'FontSize', 30, 'FontWeight', 'bold')



