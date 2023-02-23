%% SIMULATION

% import .cdt data
LoadCurryDataFile


%% Epoching EEG data based on triggers

% event length
num_triggers = size(events,2);

% init vector
six   = zeros(1, num_triggers);
seven = zeros(1, num_triggers);

% find trigger 6 & 7
for i = 1:num_triggers
    if events(2,i)  == 6
       six(i) = 1;
    elseif events(2,i) == 7
       seven(i) = 1;
    end    
end

% indeces
six_idx = find(six);
seven_idx = find(seven);

% indeces for test session 2
six_idx_4   = six_idx(end-3:end);
seven_idx_4 = seven_idx(end-3:end);

t4_idx_trig6 = six_idx_4(six_idx_4 < seven_idx_4); 
t4_idx_trig7 = seven_idx_4(six_idx_4 > seven_idx_4);
t4_idx = [t4_idx_trig6, t4_idx_trig7];

% trial indeces
trial_idx = [six_idx(1:13) seven_idx(1:13), t4_idx];
trial_idx = sort(trial_idx);


% second peep latencies (jitter corrected)
secPeep_pnt = zeros(1, length(trial_idx));
for i = 1:length(trial_idx)
    temp = events(1,trial_idx(i)+2);
    secPeep_pnt(i) = temp;
end


% load results
load('models.mat')
load('test1_Result.mat')
load('test2_Result.mat')
load('test3_Result.mat')

% f_onTrials
f_onTrials = [models.f_onTrials(1:46:end), test1_Result.f_onTrials(1:46:end), ...
              test2_Result.f_onTrials(1:46:end), test3_Result.f_onTrials(1:46:end)];


% EEG epoching w/ clipping
epoched_eeg = zeros(30, 68, 60*1000);
for i = 1:30
    temp = data(1:68, secPeep_pnt(i)+2000+f_onTrials(i):secPeep_pnt(i)+2000+f_onTrials(i)+60000-1);
    epoched_eeg(i,:,:) = temp;
end


% EEG epoching w/o clipping
epoched_eeg_off = zeros(30, 68, 60*1000);
for i = 1:30
    temp = data(1:68, secPeep_pnt(i)+2000:secPeep_pnt(i)+2000+60000-1);
    epoched_eeg_off(i,:,:) = temp;
end


% save epoched eeg
save('epoched_eeg', 'epoched_eeg')
save('epoched_eeg_off', 'epoched_eeg_off')
save('f_onTrials', 'f_onTrials')



%% import epoched data

cd '/Users/seung-cheolbaek/Desktop/ExperimentsData/OnlineAAD/EEG/pL02'

load('epoched_eeg.mat')
% load('epoched_eeg_off.mat')
load('f_onTrials.mat')

if exist('epoched_eeg_off', 'var')
    epoched_eeg = epoched_eeg_off;
    clear epoched_eeg_off
end



%% Prepare for Training and Decoding

% ################### Default Decoding Params ######################## %

% Basic params
direction     = -1; % 1 for forward, -1 for backward
sr            = 64; % sr = 64 Hz
duration      = 60;
tmin          = 0; % min time-lag(ms)
tmax          = 250; % max time-lag(ms)

% Time-lag (tau) params
tmin_index     = floor(tmin/1000*sr); % tmin2idx
tmax_index     = ceil(tmax/1000*sr);  % tmax2idx
time_lag_index = length(tmin_index:tmax_index); % time-lag points btw tmin and tmax
t              = sort( direction * ( linspace(tmin, tmax, time_lag_index) ) );

% # of eeg channels, # of trials
num_chans   = 60;
num_trials  = 30;


% ##################### Online Decoding Params ######################### %

% Sliding Window
timewin     = 15 * 1000; % in ms
timewin_idx = floor( timewin/1000*sr );

% num_frame
num_frame = duration - (timewin/1000) + 1; 

% Stride - Length as Variables
stride = 1 * 1000; % in ms
stride2idx = stride/1000;

% Lambda - The Order as Variables
lambda = 10;

% # of training trials & train_, test_idx
n_train_trial = 14;

% train_ & test_idx
all_trials            = 1:num_trials;
train_idx             = 1:14;
all_trials(train_idx) = [];
test_idx              = all_trials;
% test_idx              = 8:14;
trial_seq             = [train_idx, test_idx]; 
clear all_trials
   
% Direction
LR = {'L', 'R'};
s1_dir = [1 1 1 2 2 2, 2 1 1 2 2 1]; % 1 for L, 2 for R
s2_dir = [1 2 2 1]; % 1 for L, 2 for R
s2_jit = [31 32 34 28];
dir = [ones(1,7), ones(1,7)*2 s1_dir s2_dir]; % dir for the whole dataset


%%%%%% 16 Channels

% openbci
% openbci chans (17 chans): Fp1, Fp2, F7, F3, Fz, ...
%                           F4,  F8,  T7, C3, Cz, ...
%                           C4,  T8,  P7, P3, Pz, ...
%                           P4,  P8,  O1, O2
selec_chans = [1, 3, 6, 7, 10, ...
               12, 14, 24, 26, 28, ...
               30, 32, 44, 46, 48, ...
               50, 52, 61, 63];
num_chans   = length(selec_chans);

% from AAK
% AAK_tw15_s1_l10_tr14: F1,  F4,  FC3, Fz,  FC6,
%                       C3,  TP8, Pz,  FPz, FP2,
%                       AF3, FT7, FC5, C5,  C4,  T8
% selec_chans     = [9,  12, 17, 10, 22, ...
%                    26, 42, 48, 2,  3, ...
%                    4,  15, 16, 25, 30, 32];

% AADC_tw15_s1_l10_tr14: FC6, C5,  FC5, FC3, FC4,
%                        FT8, P7,  FP2, FCz, T7,
%                        C3,  C4,  C6,  CP3, P3, PO3
% selec_chans     = [22, 25, 16, 17, 21, ...
%                    23, 44, 3,  19, 24, ...
%                    26, 30, 31, 36, 46, 55];

% offline 16 channel selection: F5,  F3,  F1,  Fz,  F6,
%                               F8,  FT7, FC5, T7,  Pz,
%                               PO7, PO5, PO3, PO4, PO6, PO8
% selec_chans   = [7,  8,  9,  10, 13, ...
%                 14, 15, 16, 24, 48, ...
%                 53, 54, 55, 57, 58, 59];

% offline 8 channel selection: F1,  F6, F8,  FC5,
%                              PO7, PO5, PO3, PO4
% selec_chans     = [9, 13, 14, 16, ...
%                    53, 54, 55, 57];           

% 16 chans based on weights values
% selec_chans = [32, 42, 17, 16, 18, ...
%                31, 19, 7,  6,  8, ...
%                23, 24, 29, 52, 51, 53];

% % 16 chans based on weights norm values
% selec_chans = [32, 15, 42, 17, 16, ...
%                21, 18, 62, 52, 31, ...
%                35, 19, 7,  6,  26, 8];      
               
               
%% Prepare for Speech Data

% load speech data
load('AAK_right.mat'); % journey
att_speech_bank = AAK_right; 

load('AAK_left.mat'); % tweenty
unatt_speech_bank = AAK_left;

clear AAK_left; clear AAK_right;



%% Construct Band-pass filter

% ====================== Original Data Info ============================

% original sampling rate of eeg
orig_sr = 1000; % in Hz

% ===================== Filtering Params ==============================

% nyquist frequency
nyquist = orig_sr/2;

% lower & high filter bound ###############
lower_filter_bound = 0.5; % Hz
upper_filter_bound = 8; % Hz 

% transition width for filter construction
transition_width   = 0.25;

revfilt = 0; % low & band pass filter, 1 for high-pass

% ============ Create filter weights as in pop_eegfiltnew ================

edgeArray = sort([lower_filter_bound upper_filter_bound]);

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = nyquist - edgeArray(end);
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
df = min([max([edgeArray(1) * transition_width 2]) maxDf]);

% filter_order
filter_order = 3.3 / (lower_filter_bound / orig_sr); % Hamming window
filter_order = ceil(filter_order / 2) * 2; % Filter order must be even.

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;

% Window
winArray = windows('hamming', filter_order + 1);

% Filter coefficients
filterweights = firws(filter_order, cutoffArray / nyquist, winArray);


%% initialize matrix

% for left
w_l = zeros(num_chans, length(t));
b_l = zeros(1);

% for right
w_r = zeros(num_chans, length(t));
b_r = zeros(1);

% init empty vector for saving results
rL = []; rL_biased = [];
rR = []; rR_biased = [];
acc = []; acc_biased = [];

% saving pred
ReconEnv = zeros(length(test_idx), num_frame, sr*(timewin/1000));



%% Figure for Display

figure(123)
clf
set(gcf, 'color', 'w')

% Plot for BOTH Decoder
subplot(2,1,1)
h1=plot(0, 0, 'LineWidth', 1.3);
hold on
h2=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)],'ylim',[-0.5 0.5])

for k = 1:length(test_idx)
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR{dir(trial_seq(length(train_idx)+k))}, ...
            'HorizontalAlignment', 'center')
    else
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR{dir(trial_seq(length(train_idx)+k))}, ...
            'HorizontalAlignment', 'center')
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        if k > length(test_idx)-4
            h5=plot((k-1)*46+s2_jit(k-length(test_idx)+4)-14, 0, 'marker', 'h', 'MarkerSize', 10);
            h5.MarkerFaceColor = [0.9290 0.6940 0.1250];
            h5.MarkerEdgeColor = [0.9290 0.6940 0.1250];
        end
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title(' Both Decoder ', 'FontWeight', 'bold')

% Plot for Biased Decoder
subplot(2,1,2)
h3=plot(0, 0, 'LineWidth', 1.3);
hold on
h4=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)],'ylim',[-0.5 0.5])

for k = 1:length(test_idx)
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR{dir(trial_seq(length(train_idx)+k))}, ...
            'HorizontalAlignment', 'center')
    else
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR{dir(trial_seq(length(train_idx)+k))}, ...
            'HorizontalAlignment', 'center')
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        if k > length(test_idx)-4
            h5=plot((k-1)*46+s2_jit(k-length(test_idx)+4)-14, 0, 'marker', 'h', 'MarkerSize', 10);
            h5.MarkerFaceColor = [0.9290 0.6940 0.1250];
            h5.MarkerEdgeColor = [0.9290 0.6940 0.1250];
        end
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h3, h4], {'Audio From Left', 'Audio From Right'})
title(' Biased Decoder ', 'FontWeight', 'bold')



%% Training And Test

% loop over trials

for triali = 1:num_trials
    
    % loop over num_frame
    for framei = 1:46 
              
        
        % ======================= Bring EEG data ===========================
        
        eeg2use = squeeze(epoched_eeg(trial_seq(triali), :, (framei-1)*1000+1:(framei-1)*1000+15*1000));
        
        % ======================= Re-referencing ===========================
        
%         eeg2use([60 64 65 66 67 68],:) = []; % excl. CB1, CB2, HEO, VEO, EKG, EMG
        eeg2use = eeg2use(selec_chans,:); % only openbci channels
        eeg2use = bsxfun(@minus, eeg2use, mean(eeg2use, 1)); % average reference
        
        
        % ====================== Band-pass filtering ========================
        
%         eeg2use([33 43],:) = []; % excl. M1, M2
%         eeg2use = eeg2use(selec_chans, :);
        eeg2use = filtfilt(filterweights, 1, eeg2use'); % time_pnts by chans
        
        
        % ======== downsampling to 64 Hz, including an anti-alias ============
        
        eeg2use = resample(eeg2use, 64, 1000); % time by chans, if filter applied
        
        % normalize eeg data (ad hoc mean & std)
        eeg2use = bsxfun(@rdivide, bsxfun(@minus, eeg2use, mean(eeg2use, 1)), std(eeg2use, [], 1));
        
        
        
        % ===================== Speech preprocessing ==========================
        
        % speech_cursor
        speech_cursor = round((f_onTrials(trial_seq(triali))/1000)*64);
        if speech_cursor < 1
            speech_cursor = 1; end
%         speech_cursor = 1;
        
        % bring left speech data
        att_speech2use = att_speech_bank(trial_seq(triali),speech_cursor:end);  % speech clipping
        att_speech2use = [att_speech2use zeros(1,60*64-length(att_speech2use))]; % padding clipped speech part
        att_speech2use = att_speech2use((framei-1)*64+1:(framei-1)*64+15*64); % window size
        
        % bring right speech data
        unatt_speech2use = unatt_speech_bank(trial_seq(triali),speech_cursor:end); % speech clipping
        unatt_speech2use = [unatt_speech2use zeros(1,60*64-length(unatt_speech2use))]; % padding clipped speech part
        unatt_speech2use = unatt_speech2use((framei-1)*64+1:(framei-1)*64+15*64); % window size
        
        
        % 7. %%%%%%%%%%%%%%%%%%%% Training or Test %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ===================== Decoder Training ==========================
        
        if triali < length(train_idx)+1
            
            % When attended to LEFT speech
            if dir(trial_seq(triali)) == 1
                
                % assign attended speech
                attended_speech   = att_speech2use';
                
                
                % parameter training
                temp_model = mTRFtrain(attended_speech, eeg2use, ...
                    64, -1, 0, 250, lambda); % sr, dir, tmin, tmax, lambda
                
                
                % parameter updating
                w_l = w_l + temp_model.w;
                b_l = b_l + temp_model.b;
                
                
                % When attended to RIGHT speech
            elseif dir(trial_seq(triali)) == 2
                
                % assign attended speech
                attended_speech   = att_speech2use';
                
                % parameter training
                temp_model = mTRFtrain(attended_speech, eeg2use, ...
                    64, -1, 0, 250, lambda); % sr, dir, tmin, tmax, lambda
                
                % parameter updating
                w_r = w_r + temp_model.w;
                b_r = b_r + temp_model.b;
                
            end % end of dir2use if
            
            % display training phase
            disp([' Training Phase - Trial_num: ' num2str(trial_seq(triali)), ...
                  ' /  Frame_num: ' num2str(framei+14), ...
                  ' /  direction: ' LR{dir(trial_seq(triali))}])
            
            
            % construct decoder model
            if (triali == length(train_idx)) && (framei == 46)
               
                % BOTH - decoder params into structure
                model.w    = (w_l + w_r)/(46*length(train_idx)); % divide through Nframes by NtrainTrials
                model.b    = (b_l + b_r)/(46*length(train_idx)); % divide through Nframes by NtrainTrials
                model.t    = t;
                model.fs   = 64;
                model.dir  = -1;
                model.type = 'multi';
                
                % LEFT - decoder params into structure
                model_l.w    = w_l/(46*length(train_idx)/2); % divide through Nframes by NtrainTrials
                model_l.b    = b_l/(46*length(train_idx)/2); % divide through Nframes by NtrainTrials
                model_l.t    = t;
                model_l.fs   = 64;
                model_l.dir  = -1;
                model_l.type = 'multi';
                
                % RIGHT - decoder params into structure
                model_r.w    = w_r/(46*length(train_idx)/2); % divide through Nframes by NtrainTrials
                model_r.b    = b_r/(46*length(train_idx)/2); % divide through Nframes by NtrainTrials
                model_r.t    = t;
                model_r.fs   = 64;
                model_r.dir  = -1;
                model_r.type = 'multi';
                
                % display the end of training
                disp( ' End of Training!! ' )
                
            end
            
            
        % ====================== Decoder Test ===========================
        
        else
            
            % dir2use
            if triali < 27               
                dir2use = dir(trial_seq(triali));
            else
                if framei+14 < s2_jit(trial_seq(triali)-26)
                    dir2use = dir(trial_seq(triali));
                else
                    if dir(trial_seq(triali)) == 1
                        dir2use = 2;
                    elseif dir(trial_seq(triali)) == 2
                        dir2use = 1;
                    end
                end
            end % end dir2use if
            
            
            % ----------------- Direction-biased DECODER --------------------
            % When attended to LEFT speech
            if dir2use == 1
                
                % ---------------------- BOTH DECODER -------------------------
                % Decoding Left Speech
                [pred,stats] = mTRFpredict( att_speech2use', eeg2use, model );
                ReconEnv(triali-length(train_idx), framei, :) = pred;
                rL = [rL stats.acc];
                
                % Decoding Right Speech
                [~,stats] = mTRFpredict( unatt_speech2use', eeg2use, model );
                rR = [rR stats.acc];
                
                
                % ----------------- Direction-biased DECODER --------------------
                % Decoding Left Speech
                [~,stats] = mTRFpredict( att_speech2use', eeg2use, model_l );
                rL_biased = [rL_biased stats.acc];
                
                % Decoding Right Speech
                [~,stats] = mTRFpredict( unatt_speech2use', eeg2use, model_l );
                rR_biased = [rR_biased stats.acc];
                
                
                % ----------- Correctness of an Decoding Instance -------------
                
                if rL(end) > rR(end)
                    corr = 1;
                else
                    corr = 0;
                end
                acc = [acc corr];
                
                
                if rL_biased(end) > rR_biased(end)
                    corr_biased = 1;
                else
                    corr_biased = 0;
                end
                acc_biased = [acc_biased corr_biased];
                
                
                % When attended to RIGHT speech
            elseif dir2use == 2
                
                % ---------------------- BOTH DECODER -------------------------
                % Decoding Left Speech
                [~,stats] = mTRFpredict( unatt_speech2use', eeg2use, model );
                rL = [rL stats.acc];
                
                % Decoding Right Speech
                [pred,stats] = mTRFpredict( att_speech2use', eeg2use, model );
                ReconEnv(triali-length(train_idx), framei, :) = pred;
                rR = [rR stats.acc];
                
                % ----------------- Direction-biased DECODER --------------------
                
                % Decoding Left Speech
                [~,stats] = mTRFpredict( unatt_speech2use', eeg2use, model_r );
                rL_biased = [rL_biased stats.acc];
                
                % Decoding Right Speech
                [~,stats] = mTRFpredict( att_speech2use', eeg2use, model_r );
                rR_biased = [rR_biased stats.acc];
                
                
                % ----------- Correctness of an Decoding Instance -------------
                
                if rR(end) > rL(end)
                    corr = 1;
                else
                    corr = 0;
                end
                acc = [acc corr];
                
                
                if rR_biased(end) > rL_biased(end)
                    corr_biased = 1;
                else
                    corr_biased = 0;
                end
                acc_biased = [acc_biased corr_biased];
                
            end % end of dir2use if
            
            
            % Plot rL & rR
            set(h1,'XData',1:length(rL),'YData',rL)
            set(h2,'XData',1:length(rR),'YData',rR)
            
            % Plot rL_ & rR_biased
            set(h3,'XData',1:length(rL_biased),'YData',rL_biased)
            set(h4,'XData',1:length(rR_biased),'YData',rR_biased)
            
            pause(0.01)
            
            % display the latest corr. coef.
            disp(['  Test Phase - Trial_num: ' num2str(trial_seq(triali)), ...
                  ' /  Frame_num: ' num2str(framei+14), ...
                  ' /  direction: ' LR{dir2use}])
            
            disp(['rL: ' num2str(round(rL(end),2)), ...
                '  /  rR: ' num2str(round(rR(end),2)), ...
                '  /  corr:' num2str(acc(end)), ...
                '  /  corrAve: ' num2str(mean(acc))])
            
            disp(['rL_bias: ' num2str(round(rL_biased(end),2)), ...
                '  /  rR_bias: ' num2str(round(rR_biased(end),2)), ...
                '  /  corr_bias:' num2str(acc_biased(end)), ...
                '  /  corrAve_bias: ' num2str(mean(acc_biased))])
            
        end % end train or test
        
        
    end % end num_frame loop
    
    
end % end triali-loop

% save result
result = struct('model', [], 'model_l', [], 'model_r', [], ...
                'rL', [], 'rR', [], 'rL_biased', [], 'rR_biased', [], ...
                'acc', [], 'acc_biased', [], 'ReconEnv', []);
result.model = model; result.model_l = model_l; result.model_r = model_r;
result.rL = rL; result.rR = rR; result.rL_biased = rL_biased; result.rR_biased = rR_biased;
result.acc = acc; result.acc_biased = acc_biased; result.ReconEnv = ReconEnv;
            
save('result_w_clip_delta+theta_openbci.mat', 'result');
 

% mean(result.acc)
% mean(result.acc_biased)
% 
% mean(result.acc(1:368))
% mean(result.acc_biased(1:368))
% 
% mean(result.acc(369:end))
% mean(result.acc_biased(369:end))



%% Offline Baseline - preparing offline data

load('epoched_eeg_off.mat')

% init mtx for trial-wise eeg data
eeg_filt = zeros(30, 60, 60*64);

% loop over trial
for i = 1:30

    % ======================= Bring EEG data ===========================
    
    eeg2use = squeeze(epoched_eeg_off(i, :, :));
    
    % ======================= Re-referencing ===========================
    
    eeg2use([60 64 65 66 67 68],:) = []; % excl. CB1, CB2, HEO, VEO, EKG, EMG
    eeg2use = bsxfun(@minus, eeg2use, mean(eeg2use, 1)); % average reference
    
    
    % ====================== Band-pass filtering ========================
    
    eeg2use([33 43],:) = []; % excl. M1, M2
    eeg2use = filtfilt(filterweights, 1, eeg2use'); % time_pnts by chans
    
    
    % ======== downsampling to 64 Hz, including an anti-alias ============
    
    eeg2use = resample(eeg2use, 64, 1000); % time by chans, if filter applied
    
    % normalize eeg data
    eeg2use = bsxfun(@rdivide, bsxfun(@minus, eeg2use, mean(eeg2use, 1)), std(eeg2use, [], 1));

    % ======================= eeg_filt ===========================
    eeg_filt(i,:,:) = eeg2use';
    
end

save('eeg_filt', 'eeg_filt');



%% Offline Baseline

% load eeg data
load('eeg_filt.mat');

% load speech data
load('AAK_left.mat');
unattended_speech_bank = AAK_left;

load('AAK_right.mat');
attended_speech_bank = AAK_right;

clear AAK_left; clear AAK_right;


% Default Params settings for Decoder
direction     = -1; % 1 for forward, -1 for backward
sampling_rate = 64;
duration      = 60;
tmin          = 0; % min time-lag(ms)
tmax          = 250; % max time-lag(ms)
lambda        = 10;

% Time-lag (tau) params
tmin_index     = floor(tmin/1000*sampling_rate); % tmin2idx
tmax_index     = ceil(tmax/1000*sampling_rate);  % tmax2idx
time_lag_index = length(tmin_index:tmax_index); % time-lag points btw tmin and tmax
t              = sort( direction * ( linspace(tmin, tmax, time_lag_index) ) );

% # of eeg channels, # of trials
num_chans   = 60;
num_trials  = 26; % * note!!!: exclude switch data

% direction
s1_dir = [1 1 1 2 2 2, 2 1 1 2 2 1]; % 1 for L, 2 for R
dir = [ones(1,7), ones(1,7)*2 s1_dir]; % dir for the whole dataset


% initialize cell array
eeg_cell        = cell( num_trials, 1 ); % cell for eeg data
attended_cell   = cell( num_trials, 1 ); % cell for attended speech
unattended_cell = cell( num_trials, 1 ); % cell for unattended speech


% put each trial data into the corresponding cell
% loop over triali loop
for triali = 1:num_trials
    eeg_cell{triali}        = squeeze( eeg_filt(triali, :, :) )';
    attended_cell{triali}   = attended_speech_bank(triali,:)';
    unattended_cell{triali} = unattended_speech_bank(triali,:)';
end % end of triali-loop

 
% % Left trials
% L_trials = find(abs(dir-2));
% 
% % loop over triali loop
% for i = 1:length(L_trials)
%     eeg_cell{i}        = squeeze( eeg_filt(L_trials(i), :, :) )';
%     attended_cell{i}   = left_speech_bank(L_trials(i),:)';
%     unattended_cell{i} = right_speech_bank(L_trials(i),:)';
% end % end of triali-loop


% % Right trials
% R_trials = find(dir-1);
% 
% % loop over triali loop
% for i = 1:length(R_trials)
%     eeg_cell{i}        = squeeze( eeg_filt(R_trials(i), :, :) )';
%     attended_cell{i}   = right_speech_bank(R_trials(i),:)';
%     unattended_cell{i} = left_speech_bank(R_trials(i),:)';
% end % end of triali-loop


% ATTENDED DECODER
% cross-validation for optimizing the value of lambda
[~,stats_aa,stats_au,~] = mTRFattncrossval(attended_cell, unattended_cell, eeg_cell, ...
    sampling_rate, direction, tmin, tmax, lambda);

% save the results
acc = stats_aa.acc' > stats_au.acc';

r_aa = stats_aa.acc';
r_au = stats_au.acc';


% % save the results in the structure form
% offline_result.acc = acc;
% 
% offline_result.r_aa = r_aa;
% offline_result.r_au = r_au;
% 
% 
% filename = 'offline_result.mat';
% save(filename, 'offline_result')

disp([' Decoding Accuracy: ' num2str(mean(acc)*100) ' (%)'])

