%% Live Streaming - Chunked Data

%% Online Decoding Simulation
% Written by Seung-Cheol Baek
% Date: Thr. July, 16th. 2020

%% Load Data

% set working directory
cd 'C:\Users\Seung-Cheol\Desktop\OnlineDecodingSimulation_200715'

% load data, EEG & SPEECH normalized
load('DATA(3conditions).mat');

% % AAK
% ALLEEG    = DATA(1).EEG;
% ALLSPEECH = DATA(1).SPEECH;
% ALLINDEX  = DATA(1).INDEX;
% dataset   = 'AAK';
% clear DATA

% % AAKHL
% ALLEEG    = DATA(2).EEG;
% ALLSPEECH = DATA(2).SPEECH;
% ALLINDEX  = DATA(2).INDEX;
% clear DATA

% AADC
ALLEEG    = DATA(3).EEG;
ALLSPEECH = DATA(3).SPEECH;
ALLINDEX  = DATA(3).INDEX;
dataset   = 'AADC';
clear DATA


%% Prepare for Decoder Construction

% pick up one subject
sub_num  = 9;
if strcmpi(dataset, 'AAK')
    rSub_num = {char('01'), char('03'), char('04'), char('05'), char('06'), ...
        char('07'), char('08'), char('09'), char('10'), char('11')};
elseif strcmpi(dataset, 'AADC')
    rSub_num = {char('02'), char('03'), char('04'), char('05'), char('06'), ...
        char('07'), char('08'), char('09'), char('10')};
end

% ################### Default Decoding Params ######################## %

% Basic params
direction     = -1; % 1 for forward, -1 for backward
sr            = ALLEEG(sub_num).fs; % sr = 64 Hz
duration      = ALLEEG(sub_num).duration;
tmin          = 0; % min time-lag(ms)
tmax          = 250; % max time-lag(ms)

% Time-lag (tau) params
tmin_index     = floor(tmin/1000*sr); % tmin2idx
tmax_index     = ceil(tmax/1000*sr);  % tmax2idx
time_lag_index = length(tmin_index:tmax_index); % time-lag points btw tmin and tmax
t              = sort( direction * ( linspace(tmin, tmax, time_lag_index) ) );

% # of eeg channels, # of trials
num_chans   = size( ALLEEG(sub_num).data, 2 );
num_trials  = size( ALLEEG(sub_num).data, 1 );


% ##################### Online Decoding Params ######################### %

% Sliding Window
timewin     = 15 * 1000; % in ms
timewin_idx = floor( timewin/1000*sr );

% Stride - Length as Variables
stride = 1 * 1000; % in ms
stride2idx = stride/1000;

% Lambda - The Order as Variables
lambda = 10;

% # of training trials & train_, test_idx ( !!!!!!!!!!! Check)
n_train_trial = 14;

if strcmpi(dataset, 'AAK')
    % ################ for AAK #################
    % train_ & test_idx
    all_trials            = 1:num_trials;
    train_idx             = 1:n_train_trial;
    all_trials(train_idx) = [];
    test_idx              = all_trials;
    clear all_trials
    
elseif strcmpi(dataset, 'AADC')
    % ################ for AADC #################
    
    % left_idx
    left_idx = find(ALLINDEX(sub_num).left_is_1);
    left_idx_test = left_idx(8:end);
    left_idx = left_idx(1:7);
    
    % right_idx
    right_idx = find(1-ALLINDEX(sub_num).left_is_1);
    right_idx_test = right_idx(8:end);
    right_idx = right_idx(1:7);
    
    % both_idx & left_idx & right_idx & test_idx
    both_idx = sort([left_idx; right_idx]);
    
    % test_idx
    test_idx = [left_idx_test; right_idx_test];
    test_idx = sort(test_idx);
%     test_idx = [20 23 24 29 30]; % AADC06 for attention switching simul.
    test_idx'
end

% % ################ IncludeTrain #################
IncludeTrain = 0;


%% Prepare for Speech Data

% initialize list
att_speech_bank_long = [];
unatt_speech_bank_long = [];

if strcmpi(dataset, 'AAK')    
    % ########## One-sided Attention Task ###########
    if IncludeTrain % True
        % fetch Attended & Unattended speech
        attended_speech_bank    = squeeze( ALLSPEECH( ALLINDEX(sub_num).a, : ) );
        attended_speech_bank    = attended_speech_bank([train_idx test_idx], :);
        
        unattended_speech_bank  = squeeze( ALLSPEECH( ALLINDEX(sub_num).u, : ) );
        unattended_speech_bank  = unattended_speech_bank([train_idx test_idx], :);
    else % False
        % fetch Attended & Unattended speech
        attended_speech_bank    = squeeze( ALLSPEECH( ALLINDEX(sub_num).a, : ) );
        attended_speech_bank    = attended_speech_bank(test_idx, :);
        
        unattended_speech_bank  = squeeze( ALLSPEECH( ALLINDEX(sub_num).u, : ) );
        unattended_speech_bank  = unattended_speech_bank(test_idx, :);
    end 
elseif strcmpi(dataset, 'AADC')
    % ########## Attention Switching Task - AADC ###########
    % att. denotes left side & unatt. right side
    
    % init Attended & Unattended Speech Bank Matrices
    attended_speech_bank   = zeros(num_trials, sr*duration);
    unattended_speech_bank = zeros(num_trials, sr*duration);
    
    % fetch Attended & Unattended speech
    for triali = 1:num_trials
        % if left_is_1 is true
        if ALLINDEX(sub_num).left_is_1(triali) == 1
            % when subject attended to the left side
            attended_speech_bank(triali,:)   = ALLSPEECH( ALLINDEX(sub_num).a( triali ), : ); % left-left
            unattended_speech_bank(triali,:) = ALLSPEECH( ALLINDEX(sub_num).u( triali ), : ); % left-right
        else
            % when subject attended to the right side
            attended_speech_bank(triali,:)   = ALLSPEECH( ALLINDEX(sub_num).u( triali ), : ); % right-left
            unattended_speech_bank(triali,:) = ALLSPEECH( ALLINDEX(sub_num).a( triali ), : ); % right-right
        end
    end % end of triali-loop
    
    % Sort Attended & Unattended Speech Corresponding to the test trials
    attended_speech_bank    = attended_speech_bank(test_idx, :);
    unattended_speech_bank  = unattended_speech_bank(test_idx, :);
    
%     % for Attention Switching Simulation
%     switch_pnt = size(attended_speech_bank,2)/2; % at 30 second.
%     attended_speech_switch    = [attended_speech_bank(:,1:switch_pnt) unattended_speech_bank(:,switch_pnt+1:end)];
%     unattended_speech_switch  = [unattended_speech_bank(:,1:switch_pnt) attended_speech_bank(:,switch_pnt+1:end)];
%     attended_speech_bank   = attended_speech_switch;
%     unattended_speech_bank = unattended_speech_switch;
    
    % Left & Right
    LR = ['R', 'L'];
    LR_idx = ALLINDEX(sub_num).left_is_1+1;    
end

% concat data
for i = 1:size(attended_speech_bank, 1)
    att_speech_bank_long   = [att_speech_bank_long    attended_speech_bank(i,:)];
    unatt_speech_bank_long = [unatt_speech_bank_long  unattended_speech_bank(i,:)];
end

% delete EEG & Speech data from the workspace
clear ALL*
clear attended_speech_bank; clear unattended_speech_bank;


%% Decoder Model

if IncludeTrain % True
    % ################ When training model from scratch #####################
    % initialize decoder weight
    weight = zeros(num_chans, time_lag_index);
    bias   = 0;
else % False
    % ################### When using pre-trained params #####################
    if strcmpi(dataset, 'AAK')
        
        % *********** Revision in needed in this part ****************
        
    elseif strcmpi(dataset, 'AADC')
        %%%%%%%%%%%%%% AADC6 %%%%%%%%%%%%%%%%
        % load pretrained params
        load('AADC_params.mat')
        
        % ==================== BOTH DECODER ========================
        % Structure for ATTENDED decoder (use pretrained params)
        model.w    = squeeze(mean(mean(AADC_params(sub_num).weights(both_idx, :, :, :), 2), 1));
        model.b    = mean(mean(AADC_params(sub_num).biases(both_idx, :, :, :), 2));
        model.t    = t;
        model.fs   = sr;
        model.dir  = direction;
        model.type = 'multi';
        
        % ==================== LEFT DECODER ========================
        % Structure for ATTENDED decoder (use pretrained params)
        model_l.w    = squeeze(mean(mean(AADC_params(sub_num).weights(left_idx, :, :, :), 2), 1));
        model_l.b    = mean(mean(AADC_params(sub_num).biases(left_idx, :), 2));
        model_l.t    = t;
        model_l.fs   = sr;
        model_l.dir  = direction;
        model_l.type = 'multi';
        
        % ==================== RIGHT DECODER ========================
        % Structure for ATTENDED decoder (use pretrained params)
        model_r.w    = squeeze(mean(mean(AADC_params(sub_num).weights(right_idx, :, :, :), 2), 1));
        model_r.b    = mean(mean(AADC_params(sub_num).biases(right_idx, :), 2));
        model_r.t    = t;
        model_r.fs   = sr;
        model_r.dir  = direction;
        model_r.type = 'multi';
        
        % delete loaded params
        clear AADC_params
    end
end


%% Construct Band-pass filter

% ====================== Original Data Info ============================

% original sampling rate of eeg
orig_sr = 1000; % in Hz

% ===================== Filtering Params ==============================

% nyquist frequency
nyquist = orig_sr/2;

% lower & high filter bound
lower_filter_bound = 2; % Hz
upper_filter_bound = 8; % Hz

% transition width for filter construction
transition_width   = 0.25;

revfilt = 0; % low & band pass filter, 1 for high-pass

% ============ Create filter weights as in pop_eegnewfilt ================

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


%% Prepare for Streaming

% ===================== Streaming Params ==============================

% Downsampling Params
multiplier = floor(orig_sr/sr);

% Electrode indices 2 exclude from re-referencing
exc_reref = [60 64 65 66 67 68]; % CB1, CB2, HEO, VEO, EKG, EMG

% Electrode indices 2 exclude from AAD
exc_aad = [33 43]; % M1, M2

% ====================== initialization ================================

% initialize eeg ralated params
eeg_buffer        = [];  
time              = [];
eeg_data          = [];
eeg_cursor        = 0;

% initialize speech related params
attended_buffer   = [];
unattended_buffer = [];
speech_cursor     = 0;

% initialize list for saving correlation coefficients (both)
r_aa   = [];
r_au   = [];

% initialize list for saving correlation coefficients (biased)
if strcmpi(dataset, 'AADC')
    r_aa2    = [];
    r_au2    = [];
    
    r_att    = [];
    r_unatt  = [];
    
    r_att2   = [];
    r_unatt2 = [];
end

% initialize matrix for recording irregular data streaming
ir_time = [];

% condition for breaking while-loop
if IncludeTrain
    whole_data = num_trials*duration;
else
    whole_data = length(test_idx)*duration;
end

% ======================= figure 2 plot rs =============================

% figure to plot correlation coefficients
figure(1)
clf
set(gcf, 'Color', 'w');
h1=plot(0, 0, 'LineWidth', 1.3);
hold on
h2=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
for k = 1:length(test_idx)
    if k == 1
        if strcmpi(dataset, 'AAK')    
            text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
                'HorizontalAlignment', 'center')
        elseif strcmpi(dataset, 'AADC')    
            text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Num: ' num2str(test_idx(k)) '' LR(LR_idx(test_idx(k)))], ...
                'HorizontalAlignment', 'center')            
        end
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        if strcmpi(dataset, 'AAK')    
            text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
                'HorizontalAlignment', 'center')
        elseif strcmpi(dataset, 'AADC')    
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Num: ' num2str(test_idx(k)) '' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
        end
    end
end
xlabel('Time (sec.)'); ylabel('Correlation Coefficients (both) ')
if strcmpi(dataset, 'AAK')
    legend({'r\_aa', 'r\_au'})
elseif strcmpi(dataset, 'AADC')
    legend([h1, h2], {'Audio From Left', 'Audio From Right'})
end

% another figure
if strcmpi(dataset, 'AADC')
    % figure to plot correlation coefficients
    figure(2)
    clf
    set(gcf, 'Color', 'w');
    h3=plot(0, 0, 'LineWidth', 1.3);
    hold on
    h4=plot(0, 0, 'LineWidth', 1.3);
    grid on
    set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
    for k = 1:length(test_idx)
        if k == 1
            if strcmpi(dataset, 'AAK')
                text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
                    'HorizontalAlignment', 'center')
            elseif strcmpi(dataset, 'AADC')
                text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Num: ' num2str(test_idx(k)) '' LR(LR_idx(test_idx(k)))], ...
                    'HorizontalAlignment', 'center')
            end
        else
            plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
            if strcmpi(dataset, 'AAK')
                text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
                    'HorizontalAlignment', 'center')
            elseif strcmpi(dataset, 'AADC')
                text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Num: ' num2str(test_idx(k)) '' LR(LR_idx(test_idx(k)))], ...
                    'HorizontalAlignment', 'center')
            end
        end
    end
    xlabel('Time (sec.)'), ylabel('Correlation Coefficients (biased) ')
    legend([h3, h4], {'Audio From Left', 'Audio From Right'})
end

% ###################### normalize parameter ############################
% load('AAK01_mean.mat')
% z_means = AAK01_mean;
% clear AAK01_mean
% 
% load('AAK01_std.mat')
% z_stds  = AAK01_std;
% clear AAK01_std

%% Live Streaming

% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG'); end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});


%%%%%%%%%%%%%%%%% only for test (trials respectively) %%%%%%%%%%%%%%%%%%%%

% initialize count
i = -1; % considering blank phase (bc. time.sleep of the outlet)

% while loop
disp('Now receiving chunked data...');
while true
    
    % count
    i = i + 1;
    
    % ######################## Data Acquisition ##########################
    
    % measure transmission time
    tic; % TIC, pair 1
    
    % get chunk from the inlet
    [chunk, stamps] = inlet.pull_chunk();
    
    T1 = toc; % TOC, pair 1

    
    % measure data processing time
    t_process = tic; % TIC pair 2

    % ######################### Warm up ################################
    
    % save data when there are data points
    if i == 0 % no data is streamed into Matlab
        T2 = toc; % TOC pair 2
        % wait for the next data chunk
        if T2 <= 1
            pause(1-T2+T1) % considering transmission jitter
        end
        % display the process
        disp([ 'Phase Number: 0 , Processing Time: ' num2str(T2)])
    else
        
        % ################### Data Buffer & Cursor ##################### %
        
        % specify input data points
        new_pnts        = size(chunk, 2);
        % update eeg points
        eeg_cursor      = eeg_cursor + new_pnts/orig_sr;
        % record irregular transmission
        if ~(new_pnts == orig_sr)
           ir_time = [ir_time i];  
        end
        % save eeg data
        eeg_buffer = [eeg_buffer chunk];
        time       = [time stamps];
        
        % specify data points for speech
        new_pnts4speech = round((new_pnts/orig_sr)*sr);
        % save speech data, in new points are not zero
        if ~(new_pnts4speech == 0)
            if speech_cursor+new_pnts4speech <= length(att_speech_bank_long) % prevent index error
                attended_buffer   = [attended_buffer    att_speech_bank_long(speech_cursor+1:speech_cursor+new_pnts4speech)];
                unattended_buffer = [unattended_buffer  unatt_speech_bank_long(speech_cursor+1:speech_cursor+new_pnts4speech)];
            else
                attended_buffer   = [attended_buffer    att_speech_bank_long(speech_cursor+1:end)];
                unattended_buffer = [unattended_buffer  unatt_speech_bank_long(speech_cursor+1:end)];
            end
        end
        % updata speech_pnts
        speech_cursor = speech_cursor + new_pnts4speech;
        
        % Data Processing Condition
        if eeg_cursor == 0

            T2 = toc; % TOC pair 2
            
            % wait for the next data chunk
            if T2 <= 1
                pause(1-T2+T1+0.0163) % considering transmission jitter
            end
        
        elseif (mod(eeg_cursor, 60) > 0) && (mod(eeg_cursor, 60) < 15) % for the first 14 seconds in each trial
            
            T2 = toc; % TOC pair 2
            
            % display the process
            if IncludeTrain && (eeg_cursor < n_train_trial*duration)
                disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration)), ...
                    ', Train Phase: ' num2str(eeg_cursor), ...
                    ' Waiting for enough buffer size,  Processing Time: ' num2str(T2)])
            else
                if IncludeTrain
                    disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration)), ...
                        ', Test Phase: ' num2str(eeg_cursor), ...
                        ', r_aa: ' num2str(NaN) ', r_au: ' num2str(NaN), ...
                        ', Correctness: ' num2str(NaN) ', Processing Time: ' num2str(T2)])
                else
                    disp([' Trial Number: ' num2str(test_idx(ceil(eeg_cursor/duration))), ...
                        ', Test Phase: ' num2str(eeg_cursor), ...
                        ', r_aa: ' num2str(NaN) ', r_au: ' num2str(NaN), ...
                        ', Correctness: ' num2str(NaN) ', Processing Time: ' num2str(T2)])
                end
            end
            
            % wait for the next data chunk
            if T2 <= 1
                pause(1-T2+T1+0.0163) % considering transmission jitter
            end

        else     % after the first 14 seconds in each trial
            
            % ########## Update Buffer & Keep Buffer Size Uniform ###########
        
            % keep the size of the buffer no bigger than the pre-specified size
            if size(eeg_buffer, 2) > timewin
            
                % update eeg_buffer & time stamps
                eeg_buffer = eeg_buffer(:, new_pnts+1:end);
                time       = time(new_pnts+1:end);
                
                % update speech buffer
                attended_buffer   = attended_buffer(new_pnts4speech+1:end);
                unattended_buffer = unattended_buffer(new_pnts4speech+1:end);
            
            end % end size if
            
            % #################### Data Preprocessing ########################
            
            % =================== re-referencing ====================
            eeg2use = eeg_buffer;
            eeg2use(exc_reref,:) = [];
            
            if strcmpi(dataset, 'AAK')
                eeg2use = (eeg2use - repmat(mean(eeg2use,1), [62,1]));
            end
            eeg2use(exc_aad,:) = []; % result: re-referenced, 60 by time_pnts

            % ================== Band-pass filtering ==================
            if strcmpi(dataset, 'AAK')
                eeg2use = filtfilt(filterweights, 1, eeg2use'); % time_pnts by 60
            end
            
            % =================== downsample eeg data ===================
            % downsampling to 64 Hz, including an anti-alias
            if strcmpi(dataset, 'AAK')
                eeg_resamp = resample(eeg2use, sr, orig_sr); % time by chans
            elseif strcmpi(dataset, 'AADC')
                eeg_resamp = resample(eeg2use', sr, orig_sr); % time by chans, if no filter applied
            end           

            % ==================== Z-Scoring ===========================
            % normalize eeg data (used pre-specified mean & std)
%             eeg_data = bsxfun(@rdivide, bsxfun(@minus, eeg_data, z_means(ceil(eeg_cursor/duration)+n_train_trial,:)), ...
%                         z_stds(ceil(eeg_cursor/duration)+n_train_trial,:));
            
            % normalize eeg data (ad hoc mean & std)
            eeg_data = bsxfun(@rdivide, bsxfun(@minus, eeg_resamp, mean(eeg_resamp, 1)), std(eeg_resamp, [], 1));
            
            % ====================== ETC ===========================
            % copy speech data
            attended_speech   = attended_buffer;
            unattended_speech = unattended_buffer;
            
            % zero-padding: eeg_data
            if size(eeg_data, 1) < timewin_idx
                eeg_data = [eeg_data;  zeros(timewin_idx-size(eeg_data, 1), 60)];
            end
            
            % zero-padding: speech_data
            if length(attended_speech) < timewin_idx
                attended_speech   = [attended_speech   zeros(1, timewin_idx-length(attended_speech))];
                unattended_speech = [unattended_speech zeros(1, timewin_idx-length(unattended_speech))];
            end
            
            % to make sure eeg_data and speech_data have the same size
            if size(eeg_data, 1) ~= length(attended_speech)
                % consider two cases
                if size(eeg_data, 1) > length(attended_speech)
                    attended_speech   = [attended_speech    zeros(1, size(eeg_data,1)-length(attended_speech))];
                    unattended_speech = [unattended_speech  zeros(1, size(eeg_data,1)-length(attended_speech))];
                else
                    eeg_data = [eeg_data;  zeros(length(attended_speech)-size(eeg_data,1), 60)];
                end                
            end
            
            %%%%%%%%%%%%%%%%%%%%%% Train or Test %%%%%%%%%%%%%%%%%%%%%%%%
            if IncludeTrain && (eeg_cursor <= n_train_trial*duration)
               
                % ######################## Training #############################
                % train decoder model after buffer being full
                if  ~(new_pnts == 0)% 15 sec.
                    % parameter training
                    model = mTRFtrain( attended_speech',eeg_data, sr,direction,tmin,tmax,lambda); % lambda equals 10
                    
                    % parameter updating
                    weight = weight + model.w;
                    bias   = bias   + model.b;
                end % end if
                
                if abs(eeg_cursor-n_train_trial*duration) >= 0 && abs(eeg_cursor-n_train_trial*duration) < 1
                    % devide params through # of training & update model params
                    model.w  = weight / (duration-timewin/1000+1)*n_train_trial;
                    model.b  = bias   / (duration-timewin/1000+1)*n_train_trial;
                end    
                    
                % initialize eeg_data
                eeg_data = [];
                
                T2 = toc; % TOC pair 2
                
                % display the process
                disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration)), ...
                    ', Train Phase: ' num2str(eeg_cursor), ', Processing Time: ' num2str(T2)])
                if abs(eeg_cursor-n_train_trial*duration) >= 0 && abs(eeg_cursor-n_train_trial*duration) < 1
                    disp(' End of Training ')
                end    
                
            else
                % ######################## Inference #############################

                if  ~(new_pnts == 0) % to prevent redundant inferences
                    % --------------------- Both Decoder -------------------
                    % test data & save correlation coefficients
                    [~,stats_aa] = mTRFpredict( attended_speech', eeg_data, model );
                    [~,stats_au] = mTRFpredict( unattended_speech', eeg_data, model );

                    % -------- When you use LEFT & RIGHT Decoder Both --------
                    % Right Decoder
                    if  strcmpi(dataset, 'AADC') && LR_idx(test_idx(ceil(eeg_cursor/duration))) == 1
                        % test data & save correlation coefficients
                        [~,stats_aa2] = mTRFpredict( attended_speech', eeg_data, model_r );
                        [~,stats_au2] = mTRFpredict( unattended_speech', eeg_data, model_r );
                        % Left Decoder
                    elseif strcmpi(dataset, 'AADC') && LR_idx(test_idx(ceil(eeg_cursor/duration))) == 2
                        % test data & save correlation coefficients
                        [~,stats_aa2] = mTRFpredict( attended_speech', eeg_data, model_l );
                        [~,stats_au2] = mTRFpredict( unattended_speech', eeg_data, model_l );
                    end
                end
                
                % initialize eeg_data
                eeg_data = [];
                
                % ###################### Display Results ########################

                % save rs
                r_aa  = [r_aa stats_aa.acc];
                r_au  = [r_au stats_au.acc];
                
                r_aa2 = [r_aa2 stats_aa2.acc];
                r_au2 = [r_au2 stats_au2.acc];
                
                if  strcmpi(dataset, 'AADC') && LR_idx(test_idx(ceil(eeg_cursor/duration))) == 1
                    r_att    = [r_att stats_au.acc];
                    r_unatt  = [r_unatt stats_aa.acc];
                    
                    r_att2   = [r_att2 stats_au2.acc];
                    r_unatt2 = [r_unatt2 stats_aa2.acc];
                elseif strcmpi(dataset, 'AADC') && LR_idx(test_idx(ceil(eeg_cursor/duration))) == 2
                    r_att    = [r_att stats_aa.acc];
                    r_unatt  = [r_unatt stats_au.acc];
                    
                    r_att2   = [r_att2 stats_aa2.acc];
                    r_unatt2 = [r_unatt2 stats_au2.acc];
                end
                
                % plot rs
                set(h1,'XData',1:length(r_aa),'YData',r_aa)
                set(h2,'XData',1:length(r_au),'YData',r_au)
                if strcmpi(dataset, 'AADC')
                    set(h3,'XData',1:length(r_aa2),'YData',r_aa2)
                    set(h4,'XData',1:length(r_au2),'YData',r_au2)
                end
                 
                T2 = toc; % TOC pair 2
                
                % display the process
                if IncludeTrain
                    disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration))  ', Test Phase: ' num2str(eeg_cursor), ...
                        ', r_aa: ' num2str(round(r_aa(end),3)) ', r_au: ' num2str(round(r_au(end),3)), ...
                        ', Correctness: ' num2str(r_aa(end) > r_au(end)) ', Processing Time: ' num2str(T2)])
                else
                    if strcmpi(dataset, 'AADC')
                        disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration)+n_train_trial)  ', Test Phase: ' num2str(eeg_cursor), ...
                        ', r_aa: ' num2str(round(r_aa(end),2)) '; ' num2str(round(r_aa2(end),2)), ...
                        ', r_au: ' num2str(round(r_au(end),2)) '; ' num2str(round(r_au2(end),2)), ...
                        ', Correctness: ' num2str(r_aa(end) > r_au(end)) '; ' num2str(r_aa2(end) > r_au2(end)), ...
                        ', Processing Time: ' num2str(T2)])                       
                    else    
                        disp([' Trial Number: ' num2str(ceil(eeg_cursor/duration)+n_train_trial)  ', Test Phase: ' num2str(eeg_cursor), ...
                        ', r_aa: ' num2str(round(r_aa(end),3)) ', r_au: ' num2str(round(r_au(end),3)), ...
                        ', Correctness: ' num2str(r_aa(end) > r_au(end)) ', Processing Time: ' num2str(T2)])
                    end
                end
               
            end % end of trial& test-if
            
            % wait for the next data chunk
            if T2 <= 1
                pause(1-T2+T1+0.0163) % considering transmission jitter
            end
            
       end % end of data processing condition
       
        % display data-log
        if mod(i, 5) == 0
            disp( num2str( eeg_buffer(1, 1:5000:end) ) )
        end
        
        % ##### initialize eeg_buffer & speech_buffer for every trial ####
        if abs(eeg_cursor-ceil(eeg_cursor/duration)*60) >= 0 && abs(eeg_cursor-ceil(eeg_cursor/duration)*60) < 1
            eeg_buffer = [];
            attended_buffer   = [];
            unattended_buffer = [];
        end 
        
    end % end of 'i==0' if
    
    % escape while-loop when data transmission ends
    if eeg_cursor == whole_data
        break
    end
    
end % end while-loop

% save rs
if strcmpi(dataset, 'AAK')
    result = struct('w', [], 'b', [], 'r_aa', [], 'r_au', []);
    result.w = weight; result.b = bias;
    result.r_aa = r_aa; result.r_au = r_au;
elseif strcmpi(dataset, 'AADC')
    result2 = struct('r_aa', [], 'r_au', [], 'r_aa2', [], 'r_au2', [], ...
                     'r_att', [], 'r_unatt', [], 'r_att2', [], 'r_unatt2', []);
    result2.r_aa = r_aa;   result2.r_au = r_au;
    result2.r_aa2 = r_aa2; result2.r_au2 = r_au2;
    result2.r_att = r_unatt; result2.r_unatt = r_att;
    result2.r_att2 = r_unatt2; result2.r_unatt2 = r_att2;
end

filename = [dataset rSub_num{sub_num} '_result2.mat'];
save(filename, 'result2')

mean(result2.r_att < result2.r_unatt)
mean(result2.r_att2 > result2.r_unatt2)

%%
r_aa  = result.r_aa;
r_aa2 = result.r_aa2;
r_au  = result.r_au;
r_au2 = result.r_au2;

% LEFT & RIGHT Decoder
if length(r_aa) < length(test_idx)*46
    r_aa  = [r_aa zeros(1, length(test_idx)*46 - length(r_aa))];
    r_au  = [r_au zeros(1, length(test_idx)*46 - length(r_au))];
    
    r_aa2 = [r_aa2 zeros(1, length(test_idx)*46 - length(r_aa2))];
    r_au2 = [r_au2 zeros(1, length(test_idx)*46 - length(r_au2))];
end
acc = [];
acc2 = [];
for m = 1:length(test_idx)
    if LR_idx(test_idx(m)) == 2
        acc  = [acc  r_aa((m-1)*46+1:m*46) > r_au((m-1)*46+1:m*46)];
        acc2 = [acc2 r_aa2((m-1)*46+1:m*46) > r_au2((m-1)*46+1:m*46)];
    elseif LR_idx(test_idx(m)) == 1
        acc  = [acc  r_au((m-1)*46+1:m*46) > r_aa((m-1)*46+1:m*46)];
        acc2 = [acc2 r_au2((m-1)*46+1:m*46) > r_aa2((m-1)*46+1:m*46)];
    end
end

mean(acc)
mean(acc2)


% figure to plot correlation coefficients
figure(1)
clf
set(gcf, 'Color', 'w');
h1=plot(r_aa, 'LineWidth', 1.3);
hold on
h2=plot(r_au, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
for k = 1:length(test_idx)
    if k == 1
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'r\_aa', 'r\_au'})
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title( [' AADC06 - LEFT & RIGHT Decoder for LEFT & RIGHT trials respectively  / Accuracy:  ' num2str(mean(acc)*100) ' (%) ' ], ...
    'FontSize', 15, 'FontWeight', 'bold')


% Smooth corr. coef. within each trial
r_aa_3 = zeros(1,length(r_aa));
r_au_3 = zeros(1,length(r_au));

r_aa_5 = zeros(1,length(r_aa));
r_au_5 = zeros(1,length(r_au));

r_aa_7 = zeros(1,length(r_aa));
r_au_7 = zeros(1,length(r_au));

for m = 1:length(test_idx)
    % trial data
    temp_r_aa = r_aa2((m-1)*46+1:m*46);
    temp_r_au = r_au2((m-1)*46+1:m*46);
    
    for i = 1:length(temp_r_aa)
        if i < 3
            r_aa_3((m-1)*46+i) = mean(temp_r_aa(1:i));
            r_au_3((m-1)*46+i) = mean(temp_r_au(1:i));
        else
            r_aa_3((m-1)*46+i) = mean(temp_r_aa(i-2:i));
            r_au_3((m-1)*46+i) = mean(temp_r_au(i-2:i));
        end
    end
    
     for j = 1:length(temp_r_aa)
        if j < 5
            r_aa_5((m-1)*46+j) = mean(temp_r_aa(1:j));
            r_au_5((m-1)*46+j) = mean(temp_r_au(1:j));
        else
            r_aa_5((m-1)*46+j) = mean(temp_r_aa(j-4:j));
            r_au_5((m-1)*46+j) = mean(temp_r_au(j-4:j));
        end
    end
    
     for k = 1:length(temp_r_aa)
        if k < 7
            r_aa_7((m-1)*46+k) = mean(temp_r_aa(1:k));
            r_au_7((m-1)*46+k) = mean(temp_r_au(1:k));
        else
            r_aa_7((m-1)*46+k) = mean(temp_r_aa(k-6:k));
            r_au_7((m-1)*46+k) = mean(temp_r_au(k-6:k));
        end
    end    
end

acc_3 = [];
for m = 1:length(test_idx)
    if LR_idx(test_idx(m)) == 2
        acc_3 = [acc_3 r_aa_3((m-1)*46+1:m*46) > r_au_3((m-1)*46+1:m*46)];
    elseif LR_idx(test_idx(m)) == 1
        acc_3 = [acc_3 r_au_3((m-1)*46+1:m*46) > r_aa_3((m-1)*46+1:m*46)];
    end
end

acc_5 = [];
for m = 1:length(test_idx)
    if LR_idx(test_idx(m)) == 2
        acc_5 = [acc_5 r_aa_5((m-1)*46+1:m*46) > r_au_5((m-1)*46+1:m*46)];
    elseif LR_idx(test_idx(m)) == 1
        acc_5 = [acc_5 r_au_5((m-1)*46+1:m*46) > r_aa_5((m-1)*46+1:m*46)];
    end
end

acc_7 = [];
for m = 1:length(test_idx)
    if LR_idx(test_idx(m)) == 2
        acc_7 = [acc_7 r_aa_7((m-1)*46+1:m*46) > r_au_7((m-1)*46+1:m*46)];
    elseif LR_idx(test_idx(m)) == 1
        acc_7 = [acc_7 r_au_7((m-1)*46+1:m*46) > r_aa_7((m-1)*46+1:m*46)];
    end
end

mean(acc_3)
mean(acc_5)
mean(acc_7)



figure(2)
clf
set(gcf, 'Color', 'w');

subplot(3, 1, 1)
h1=plot(r_aa_3, 'LineWidth', 1.3);
hold on
h2=plot(r_au_3, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
for k = 1:length(test_idx)
    if k == 1
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'r\_aa\_3', 'r\_au\_3'})
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title([' 3 Point Moving Average Filter  /  Accuracy: ' num2str(100*mean(acc_3)), ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')

subplot(3, 1, 2)
h1=plot(r_aa_5, 'LineWidth', 1.3);
hold on
h2=plot(r_au_5, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
for k = 1:length(test_idx)
    if k == 1
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'r\_aa\_3', 'r\_au\_3'})
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title([' 5 Point Moving Average Filter  /  Accuracy: ' num2str(100*mean(acc_5)), ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')

subplot(3, 1, 3)
h1=plot(r_aa_7, 'LineWidth', 1.3);
hold on
h2=plot(r_au_7, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*length(test_idx)+1],'ylim',[-0.5 0.5])
for k = 1:length(test_idx)
    if k == 1
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
%             'HorizontalAlignment', 'center')
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
            'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'r\_aa\_3', 'r\_au\_3'})
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title([' 7 Point Moving Average Filter  /  Accuracy: ' num2str(100*mean(acc_7)), ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')

sgtitle({' AADC06 - LEFT & RIGHT Decoder for LEFT & RIGHT trials respectively', ...
    ' Moving Average Filter Applied '}, 'FontSize', 15, 'FontWeight', 'bold')


%%

trans_pnts = zeros(1, length(test_idx));
% attention switch simulation plot
for m=1:length(test_idx)
   temp_r_aa = r_aa((m-1)*46+1:m*46);
   temp_r_au = r_au((m-1)*46+1:m*46);
   [~,idx] = min(abs(temp_r_aa(1:30) - temp_r_au(1:30)));
   trans_pnts(m) = idx;
end
trans_pnts = trans_pnts + 14;

pad = zeros(1, 60-46);
rs_aa = zeros(5, 60);
rs_au = zeros(5, 60);
for m=1:length(test_idx)
   rs_aa(m,:) = [pad r_aa((m-1)*46+1:m*46)];
   rs_au(m,:) = [pad r_au((m-1)*46+1:m*46)];
end


for i = 1:length(test_idx)
    figure(i)
    clf
    set(gcf, 'Color', 'w');
    h1=plot(rs_aa(i,:), 'LineWidth', 1.3);
    hold on
    h2=plot(rs_au(i,:), 'LineWidth', 1.3);
    h3=plot(trans_pnts(i), rs_aa(i,trans_pnts(i)), '*', 'MarkerSize', 15);
    set(gca, 'xlim', [1 60], 'ylim', [-0.3 0.3], 'FontSize', 13)
    h4=plot( [15, 15], get(gca,'ylim'), 'k--', 'LineWidth', 0.5);
    grid on; hold off;
    xlabel('Time (sec.)', 'FontSize', 18, 'FontWeight', 'bold')
    ylabel('Correlation Coefficients', 'FontSize', 18, 'FontWeight', 'bold')
    legend([h1, h2, h3, h4], {'Audio From Left', 'Audio From Right', 'Transition Point', ...
        ' Start Inference '})
    title({' AADC06 - Attention Switch Simulation, Attention switched at 30 sec. ', ...
        [' Trial Number: ' num2str(test_idx(i)) ' /  Transition Point: ' num2str(trans_pnts(i)) ' (sec.) ']}, ...
        'FontSize', 18, 'FontWeight', 'bold')
end



%%

% figure(1)
% title({' AAK01 - Online AAD Simulation w/ Data Streaming using LSL ( no Preproc & online Training )', ...
%     [' Accuracy (16 test): ' num2str(100*mean(r_aa > r_au)), ...
%     ' (%)  /  Accuracy (15 test): ' num2str(100*mean(r_aa(47:end) > r_au(47:end))) ' (%)']}, ...
%     'FontSize', 15, 'FontWeight', 'bold')

r_aa = [];
r_au = [];

for k = 1:15
    r_aa = [r_aa  all_results.results.rs_aa(k,:)];
    r_au = [r_au  all_results.results.rs_au(k,:)];
end

% figure to plot correlation coefficients
figure(1)
clf
set(gcf, 'Color', 'w');
h1=plot(r_aa, 'LineWidth', 1.3);
hold on
h2=plot(r_au, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*15+1],'ylim',[-0.5 0.5])
for k = 1:15
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend({'r\_aa', 'r\_au'})
% legend([h1, h2], {'Audio From Left', 'Audio From Right'})
% title({' AAK01 - Online AAD Simulation w/ Data Streaming using LSL ( no Preproc & pretrained )', ...
%     [' Accuracy (15 test): ' num2str(100*mean(r_aa > r_au)) ' (%)']}, ...
%     'FontSize', 15, 'FontWeight', 'bold')
title({' AAK01 - Online AAD Simulation w/o Live Streaming', ...
    [' Accuracy (15 test): ' num2str(100*mean(r_aa > r_au)) ' (%)']}, ...
    'FontSize', 15, 'FontWeight', 'bold')

% =========================== 3 filter ===================================
r_aa_3 = zeros(1,length(r_aa));
r_au_3 = zeros(1,length(r_au));

for i=1:length(r_aa)
    if i < 3
        r_aa_3(i) = mean(r_aa(1:i));
        r_au_3(i) = mean(r_au(1:i));
    else
        r_aa_3(i) = mean(r_aa(i-2:i));
        r_au_3(i) = mean(r_au(i-2:i));
    end
end

r_aa_5 = zeros(1,length(r_aa));
r_au_5 = zeros(1,length(r_au));

for i=1:length(r_aa)
    if i < 5
        r_aa_5(i) = mean(r_aa(1:i));
        r_au_5(i) = mean(r_au(1:i));
    else
        r_aa_5(i) = mean(r_aa(i-4:i));
        r_au_5(i) = mean(r_au(i-4:i));
    end
end

r_aa_7 = zeros(1,length(r_aa));
r_au_7 = zeros(1,length(r_au));

for i=1:length(r_aa)
    if i < 7
        r_aa_7(i) = mean(r_aa(1:i));
        r_au_7(i) = mean(r_au(1:i));
    else
        r_aa_7(i) = mean(r_aa(i-6:i));
        r_au_7(i) = mean(r_au(i-6:i));
    end
end

mean(r_aa > r_au)
mean(r_aa_3 > r_au_3)
mean(r_aa_5 > r_au_5)
mean(r_aa_7 > r_au_7)

% ====================== 3 filter - 15 Tr ==============================
r_aa_2 = r_aa(47:end);
r_au_2 = r_au(47:end);

r_aa_3_2 = zeros(1,length(r_aa_2));
r_au_3_2 = zeros(1,length(r_au_2));

for i=1:length(r_aa_2)
    if i < 3
        r_aa_3_2(i) = mean(r_aa_2(1:i));
        r_au_3_2(i) = mean(r_au_2(1:i));
    else
        r_aa_3_2(i) = mean(r_aa_2(i-2:i));
        r_au_3_2(i) = mean(r_au_2(i-2:i));
    end
end

r_aa_5_2 = zeros(1,length(r_aa_2));
r_au_5_2 = zeros(1,length(r_au_2));

for i=1:length(r_aa_2)
    if i < 5
        r_aa_5_2(i) = mean(r_aa_2(1:i));
        r_au_5_2(i) = mean(r_au_2(1:i));
    else
        r_aa_5_2(i) = mean(r_aa_2(i-4:i));
        r_au_5_2(i) = mean(r_au_2(i-4:i));
    end
end

r_aa_7_2 = zeros(1,length(r_aa_2));
r_au_7_2 = zeros(1,length(r_au_2));

for i=1:length(r_aa_2)
    if i < 7
        r_aa_7_2(i) = mean(r_aa_2(1:i));
        r_au_7_2(i) = mean(r_au_2(1:i));
    else
        r_aa_7_2(i) = mean(r_aa_2(i-6:i));
        r_au_7_2(i) = mean(r_au_2(i-6:i));
    end
end


figure(2)
clf
set(gcf, 'Color', 'w');

subplot(3, 1, 1)
h1=plot(r_aa_3, 'LineWidth', 1.3);
hold on
h2=plot(r_au_3, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*15+1],'ylim',[-0.5 0.5])
for k = 1:15
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend({'r\_aa\_3', 'r\_au\_3'})
% legend([h1, h2], {'Audio From Left', 'Audio From Right'})
% title([' 7 Point Moving Average Filter  /  Accuracy (16 tests): ' num2str(100*mean(r_aa_3 > r_au_3)), ...
%     ' (%)  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_3_2 > r_au_3_2)) ' (%) '], ...
%     'FontSize', 13, 'FontWeight', 'bold')
title([' 3 Point Moving Average Filter  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_3 > r_au_3)), ...
    ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')


subplot(3, 1, 2)
h1=plot(r_aa_5, 'LineWidth', 1.3);
hold on
h2=plot(r_au_5, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*15+1],'ylim',[-0.5 0.5])
for k = 1:15
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend({'r\_aa\_5', 'r\_au\_5'})
% legend([h1, h2], {'Audio From Left', 'Audio From Right'})
% title([' 5 Point Moving Average Filter  /  Accuracy (16 tests): ' num2str(100*mean(r_aa_5 > r_au_5)), ...
%     ' (%)  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_5_2 > r_au_5_2)) ' (%) '], ...
%     'FontSize', 13, 'FontWeight', 'bold')
title([' 5 Point Moving Average Filter  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_5 > r_au_5)), ...
    ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')

subplot(3, 1, 3)
h1=plot(r_aa_7, 'LineWidth', 1.3);
hold on
h2=plot(r_au_7, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*15+1],'ylim',[-0.5 0.5])
for k = 1:15
    if k == 1
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    else 
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
        text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Trial ' num2str(test_idx(k+1))], ...
            'HorizontalAlignment', 'center')
%         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, ['Tr Num: ' num2str(test_idx(k)) ' ' LR(LR_idx(test_idx(k)))], ...
%             'HorizontalAlignment', 'center')
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend({'r\_aa\_7', 'r\_au\_7'})
% legend([h1, h2], {'Audio From Left', 'Audio From Right'})
% title([' 7 Point Moving Average Filter  /  Accuracy (16 tests): ' num2str(100*mean(r_aa_7 > r_au_7)), ...
%     ' (%)  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_7_2 > r_au_7_2)) ' (%) '], ...
%     'FontSize', 13, 'FontWeight', 'bold')
title([' 7 Point Moving Average Filter  /  Accuracy (15 tests): ' num2str(100*mean(r_aa_7 > r_au_7)), ...
    ' (%) '], ...
    'FontSize', 13, 'FontWeight', 'bold')

sgtitle({' AAK01 - Online AAD Simulation w/o Live Streaming', ...
    ' Moving Average Filter Applied '}, 'FontSize', 15, 'FontWeight', 'bold')




% figure(1)
% title([' AAK01 - Online AAD Simulation w/ Data Streaming using LSL / Accuracy: ' num2str(100*mean(r_aa > r_au)) ' (%)'], ...
%     'FontSize', 15, 'FontWeight', 'bold')
% % title([' AADC6 - Online AAD Simulation w/ Data Streaming using LSL (tracking Left-Right) / Accuracy: ' num2str(100*mean(acc)) ' (%)'], ...
% %     'FontSize', 15, 'FontWeight', 'bold')
% 
% % r_aa = [r_aa 0];
% % r_au = [r_au 0];
% % 
% % acc = zeros(1, 690);
% % 
% % for k = 1:15
% %     if (LR_idx(k+15)-1) == 1
% %         acc((k-1)*46+1:k*46) = r_aa((k-1)*46+1:k*46) > r_au((k-1)*46+1:k*46);
% %     else    
% %         acc((k-1)*46+1:k*46) = r_au((k-1)*46+1:k*46) > r_aa((k-1)*46+1:k*46);
% %     end
% % end 
% % 
% % mean(acc)
% 
% r_aa = [];
% r_au = [];
% 
% for k = 1:46
%     r_aa = [r_aa  all_results.results.rs_aa(:,k)'];
%     r_au = [r_au  all_results.results.rs_au(:,k)'];
% end
% 
% mean(r_aa > r_au)
% 
% 
% 
% 
% windowSize = 3; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s3_r_aa = filter(b, a, r_aa);
% s3_r_au = filter(b, a, r_au);
% mean(s3_r_aa > s3_r_au)*100 % 0.7377
% 
% windowSize = 5; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s5_r_aa = filter(b, a, r_aa);
% s5_r_au = filter(b, a, r_au);
% mean(s5_r_aa > s5_r_au)*100 % 0.7551
% 
% windowSize = 7; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s7_r_aa = filter(b, a, r_aa);
% s7_r_au = filter(b, a, r_au);
% mean(s7_r_aa > s7_r_au)*100 % 0.7667
% 
% 
% % figure to plot correlation coefficients
% figure(2)
% clf
% set(gcf, 'Color', 'w');
% 
% subplot(3, 1, 1)
% plot(s3_r_aa, 'LineWidth', 1.3, 'Color', [0.3010 0.7450 0.9330]);
% hold on
% plot(s3_r_au, 'LineWidth', 1.3, 'Color', [0.9290 0.6940 0.1250]);
% grid on
% set(gca,'xlim',[1 700],'ylim',[-0.5 0.5])
% for k = 1:46
%     if ~(k == 1)
%         plot( [(k-1)*15+1 (k-1)*15+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
% %         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR(LR_idx(k+15)))
%     end
% end
% xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'s3\_r\_aa', 's3\_r\_au'})
% 
% title([' AAK01 ( w/o Streaming, 3 samples ) / Accuracy (wrong): ' num2str(100*mean(s3_r_aa > s3_r_au)) ' (%)'], ...
%     'FontSize', 15, 'FontWeight', 'bold')
% 
% subplot(3, 1, 2)
% plot(s5_r_aa, 'LineWidth', 1.3, 'Color', [0.3010 0.7450 0.9330]);
% hold on
% plot(s5_r_au, 'LineWidth', 1.3, 'Color', [0.9290 0.6940 0.1250]);
% grid on
% set(gca,'xlim',[1 700],'ylim',[-0.5 0.5])
% for k = 1:46
%     if ~(k == 1)
%         plot( [(k-1)*15+1 (k-1)*15+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
% %         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR(LR_idx(k+15)))
%     end
% end
% xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'s5\_r\_aa', 's5\_r\_au'})
% 
% title([' AAK01 ( w/o Streaming, 5 samples ) / Accuracy (wrong): ' num2str(100*mean(s5_r_aa > s5_r_au)) ' (%)'], ...
%     'FontSize', 15, 'FontWeight', 'bold')
% 
% subplot(3, 1, 3)
% plot(s7_r_aa, 'LineWidth', 1.3, 'Color', [0.3010 0.7450 0.9330]);
% hold on
% plot(s7_r_au, 'LineWidth', 1.3, 'Color', [0.9290 0.6940 0.1250]);
% grid on
% set(gca,'xlim',[1 700],'ylim',[-0.5 0.5])
% for k = 1:46
%     if ~(k == 1)
%         plot( [(k-1)*15+1 (k-1)*15+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
% %         text( (((k-1)*46+1) + (k*46+1))/2, -0.45, LR(LR_idx(k+15)))
%     end
% end
% xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
% legend({'s7\_r\_aa', 's7\_r\_au'})
% 
% title([' AAK01 ( w/o Streaming, 7 samples ) / Accuracy (wrong): ' num2str(100*mean(s7_r_aa > s7_r_au)) ' (%)'], ...
%     'FontSize', 15, 'FontWeight', 'bold')
% 
% 
% rs_aa = all_results.results.rs_aa;
% rs_au = all_results.results.rs_au;
% 
% windowSize = 3; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s3_r_aa = filter(b, a, rs_aa);
% s3_r_au = filter(b, a, rs_au);
% mean(s3_r_aa > s3_r_au, 'all')*100 % 0.7377
% 
% windowSize = 5; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s5_r_aa = filter(b, a, rs_aa);
% s5_r_au = filter(b, a, rs_au);
% mean(s5_r_aa > s5_r_au, 'all')*100 % 0.7551
% 
% windowSize = 7; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% s7_r_aa = filter(b, a, rs_aa);
% s7_r_au = filter(b, a, rs_au);
% mean(s7_r_aa > s7_r_au, 'all')*100 % 0.7667



% acc3 = s3_r_aa > s3_r_au
% acc3_t1 = acc3(1:46)
% 
% acc3_re = reshape(acc3, [], 15)
% acc3_re = acc3_re'
% 
% mean(mean(acc3_re, 2))


