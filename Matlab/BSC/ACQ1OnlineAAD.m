%% Online AAD
% written by Seung-Cheol Baek
% Update Date: Sat. 15th August, 2020
%

%% Things needed to be solved
%
% 1. Polish Codes - done
% 2. Correct Jitter - done
% 3. check if the revised code works - done
% 3. Re-conduct pilot test w/ Participants
%

%% Parameters for Simulation in KIST 
%  - should be changed into quote when run this code in a real experiment
%  situation

% % =================== Simulation Settings ==========================
% training_marker = 2;
% test_one_marker = 3;
% test_two_marker = 4;
% end_marker      = 5;
% 
% trial_marker    = 1200001; % make use of photo sensor
% dir_marker      = [6, 7]; % 6 for LEFT / 7 for RIGTH
% insampleratehz  = 1000;
% 
% 
% % ================= LSL Library Load ============================
% global eeg_inlet
% [indat, ~] = eeg_inlet.pull_chunk();
% 
% global marker_inlet
% [marker, marker_time] = marker_inlet.pull_chunk();
% inevents = [marker_time; marker];
 


%% Real-time AAD Part

%%%%%%%%%%%%%%%%%%% 1. Global Variable Settings %%%%%%%%%%%%%%%%%%%%%%%%

% ============== Load or Initialize Global Variables ====================

% Initialize EEG_buffer
if isempty(whos('global','eeg_buffer'))
    initGlobal_eeg_buffer; end

% Load Pre-processed Speech Data
if isempty(whos('global','left_speech_bank')) || isempty(whos('global','right_speech_bank'))
    setGlobal_speech_bank; end


% ============== Initialize Toggles for Each Session ====================

% Initialize IsTrain
if isempty(whos('global','IsTrain'))
    initGlobal_IsTrain; end
global IsTrain

% Initialize IsTestOne
if isempty(whos('global','IsTestOne'))
    initGlobal_IsTestOne; end
global IsTestOne

% Initialize IsTestTwo
if isempty(whos('global','IsTestTwo'))
    initGlobal_IsTestTwo; end
global IsTestTwo

% Initialize Trial Count
if isempty(whos('global','trial_count'))
    initGlobal_trial_count; end

% Initialize dir
if isempty(whos('global','dir'))
    initGlobal_dir; end



%%%%%%%%%%%%% 2. On & Off Train, Test, and Trial Toggles  %%%%%%%%%%%%%%%

% ===================== On & Off IsTrain Toggels =======================

% initialize global 'IsTrain' when there is training_marker in the chunk
if any(ismember(inevents(2,:), 2)) && (IsTrain == 0) % training_marker = 2
    setGlobal_IsTrain; % init weights & bias
elseif any(ismember(inevents(2,:), 5)) && (IsTrain == 1) % end_marker = 5
    % off IsTrain when training_marker come out one more time
    offGlobal_IsTrain; % save decoder model
end
IsTrain2use = IsTrain; % copy IsTrain


% ===================== On & Off IsTestOne Toggels =========================

% initialize global 'IsTestOne' when there is test_one_marker in the chunk
if any(ismember(inevents(2,:), 3)) && (IsTestOne == 0)  % test_one_marker = 3
    setGlobal_IsTestOne; % init rL & rR; rL_biased & rR_biased
elseif any(ismember(inevents(2,:), 5)) && (IsTestOne == 1) % end_marker = 5
    % off IsTest when test_two_marker come out one more time
    offGlobal_IsTestOne; % save r_left & r_right
end
IsTestOne2use = IsTestOne;  % copy IsTest


% ===================== On & Off IsTestTwo Toggels =========================

% initialize global 'IsTestOne' when there is test_two_marker in the chunk
if any(ismember(inevents(2,:), 4)) && (IsTestTwo == 0)  % test_two_marker = 4
    setGlobal_IsTestTwo; % init r_left & r_right
elseif any(ismember(inevents(2,:), 5)) && (IsTestTwo == 1) % end_marker = 5
    % off IsTest when test_two_marker come out one more time
    offGlobal_IsTestTwo; % save r_left & r_right
end
IsTestTwo2use = IsTestTwo;  % copy IsTest



%%%%%%%%%%%%%%%%%%%%%%%% 3. Update EEG buffer %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update eeg_buffer
global eeg_buffer
eeg_buffer = [eeg_buffer indat];

% ====== when the size of eeg_buffer is bigger than time window =======

if size(eeg_buffer,2) > 15*insampleratehz
    eeg_buffer = eeg_buffer(:,end-15*insampleratehz+1:end); % trim eeg_buffer
    eeg2use = eeg_buffer; % copy eeg_buffer for processing
end


    
%%%%%%%%%%%%%%%%%%%%%%%% 4. Trial Markers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================== Direction Marker ============================

% initialize global LeftIsOne when there is direction-related marker in the chunk
% dir_marker(1) = left_marker = 6 / dir_marker(2) = right_marker = 7
if any(ismember(inevents(2,:), [6,7]))
    setGlobal_dir(inevents, [6,7]); end
global dir;
dir2use = dir;

% direction marker
LR = {'NaN', 'L', 'R'};


% ======================== Trial Onset Marker ============================

if (IsTrain2use || IsTestOne2use || IsTestTwo2use)
    if (dir2use ~= 0) && isempty(whos('global','onTrial')) && any(ismember(inevents(2,:), 1000001)) % audio_marker = 1000001
        % initialize global onTrial when there audio stimulus starts
        setGlobal_onTrial(inevents, insampleratehz);
    end
end

% update marker info when there is onTrial
if ~isempty(whos('global','onTrial'))
    % onTrial Var.
    global onTrial
    onTrial     = onTrial + 1; % update process
    onTrial2use = onTrial; % copy onTrial to use in this script
    
    % flag onTrial
    global f_onTrial
    f_onTrial2use = f_onTrial; % copy flag to use in this script
    
    % trial_count Var.
    global trial_count
    trial_count2use = trial_count;
end


%%%%%%%%%%%%%%%%%%%%%%%% 5. Display Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(0,0)
set(gca, 'xlim', [-12 12], 'ylim', [-3 3])
text(0,2, ['Training Phase: ' num2str(IsTrain2use), ...
    '  /  Test 1 Phase: ' num2str(IsTestOne2use), ...
    '  /  Test 2 Phase: ' num2str(IsTestTwo2use)], ...
    'HorizontalAlignment', 'center')
if ~isempty(whos('onTrial2use'))
    text(0,1, ['Train Trial\_num: ' num2str(trial_count2use), ...
        '  /  Frame\_num: ' num2str(onTrial2use), ...
        '  /  Direction: ' num2str(LR{dir2use+1})], ...
        'HorizontalAlignment', 'center')
end
if ~isempty(inevents) % to check marker input
    text(0,0, ['Input Markers: ' num2str(inevents(2,:)), ...
        ' / f\_onTrial: ' num2str(f_onTrial2use)], ...
        'HorizontalAlignment', 'center')
end



% %%%%%%%%%%%%%%%%%%%% 6. Data Preprocessing Part %%%%%%%%%%%%%%%%%%%%%%%%%

% ================= conditions for data preprocessing =====================

if exist('onTrial2use', 'var') && onTrial2use >= 15
    
    
    % ======================= EEG preprocessing ===========================
   
    % copy eeg_buffer for use
    eeg2use = double(eeg2use);

    
    % ======================= Re-referencing ===========================
    
    eeg2use([60 64 65 66 67 68 69],:) = []; % excl. CB1, CB2, HEO, VEO, EKG, EMG, Trigger 
    eeg2use = bsxfun(@minus, eeg2use, mean(eeg2use, 1)); % average reference
   
    
    % ====================== Band-pass filtering ========================
    
    eeg2use([33 43],:) = []; % excl. M1, M2
    
    global filterweights;
    eeg2use = filtfilt(filterweights, 1, eeg2use'); % time_pnts by chans
    
    
    % ======== downsampling to 64 Hz, including an anti-alias ============
    
    eeg2use = resample(eeg2use, 64, insampleratehz); % time by chans, if filter applied
%     eeg2use = resample(eeg2use', 64, insampleratehz); % time by chans, if no filter applied
    
    % normalize eeg data (ad hoc mean & std)
    eeg2use = bsxfun(@rdivide, bsxfun(@minus, eeg2use, mean(eeg2use, 1)), std(eeg2use, [], 1));
    
    
    % ===================== Speech preprocessing ==========================
    
    % speech_cursor
    speech_cursor = round(f_onTrial2use*64);
    
    % bring left speech data
    global left_speech_bank
    left_speech2use = left_speech_bank(trial_count2use,speech_cursor+1:end);  % speech clipping
    left_speech2use = [left_speech2use zeros(1,60*64-length(left_speech2use))]; % padding clipped speech part
    left_speech2use = left_speech2use((onTrial2use-15)*64+1:(onTrial2use-15)*64+15*64); % window size
    
    % bring right speech data
    global right_speech_bank
    right_speech2use = right_speech_bank(trial_count2use,speech_cursor+1:end); % speech clipping
    right_speech2use = [right_speech2use zeros(1,60*64-length(right_speech2use))]; % padding clipped speech part
    right_speech2use = right_speech2use((onTrial2use-15)*64+1:(onTrial2use-15)*64+15*64); % window size
        
    
    % 7. %%%%%%%%%%%%%%%%%%%% Training or Test %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ===================== Decoder Training ==========================
    
    if IsTrain2use
        
        % When attended to LEFT speech 
        if dir2use == 1
            
            % recall global vars.
            global w_l; global b_l;
                        
            % assign attended speech
            attended_speech   = left_speech2use';
            
            % parameter training
            temp_model = mTRFtrain(attended_speech, eeg2use, ...
                64, -1, 0, 250, 10); % sr, dir, tmin, tmax, lambda
            
            % parameter updating
            w_l = w_l + temp_model.w;
            b_l = b_l + temp_model.b;
            
        % When attended to RIGHT speech    
        elseif dir2use == 2
            
            % recall global vars.
            global w_r; global b_r;
            
            % assign attended speech
            attended_speech   = right_speech2use';
            
            % parameter training
            temp_model = mTRFtrain(attended_speech, eeg2use, ...
                64, -1, 0, 250, 10); % sr, dir, tmin, tmax, lambda
            
            % parameter updating
            w_r = w_r + temp_model.w;
            b_r = b_r + temp_model.b;
            
        end % end of dir2use if
                
    end % end IsTrain2use if
    
    
    % ====================== Decoder Test ===========================
    
    if (IsTestOne2use || IsTestTwo2use)
        
        % call decoder structure
        global model; global model_l; global model_r
        model2use = model;
        model_l2use = model_l;
        model_r2use = model_r;
        
        % call global r_left & r_right
        global rL; global rL_biased;
        global rR; global rR_biased;
        
        % call accs
        global acc; global acc_biased;
                
        
        % ----------------- Direction-biased DECODER --------------------
        % When attended to LEFT speech
        if dir2use == 1
            
            % ---------------------- BOTH DECODER -------------------------
            % Decoding Left Speech
            [~,stats] = mTRFpredict( left_speech2use', eeg2use, model2use );
            rL = [rL stats.acc];
            
            % Decoding Right Speech
            [~,stats] = mTRFpredict( right_speech2use', eeg2use, model2use );
            rR = [rR stats.acc];
   
            
            % ----------------- Direction-biased DECODER --------------------
            % Decoding Left Speech
            [~,stats] = mTRFpredict( left_speech2use', eeg2use, model_l2use );
            rL_biased = [rL_biased stats.acc];
            
            % Decoding Right Speech
            [~,stats] = mTRFpredict( right_speech2use', eeg2use, model_l2use );
            rR_biased = [rR_biased stats.acc];
            
            
            % ----------- Correctness of an Decoding Instance -------------
            % ?? modification needed for corr. smoothing
            
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
            [~,stats] = mTRFpredict( left_speech2use', eeg2use, model2use );
            rL = [rL stats.acc];
            
            % Decoding Right Speech
            [~,stats] = mTRFpredict( right_speech2use', eeg2use, model2use );
            rR = [rR stats.acc];
            
            % ----------------- Direction-biased DECODER --------------------
            
            % Decoding Left Speech
            [~,stats] = mTRFpredict( left_speech2use', eeg2use, model_r2use );
            rL_biased = [rL_biased stats.acc];
            
            % Decoding Right Speech
            [~,stats] = mTRFpredict( right_speech2use', eeg2use, model_r2use );
            rR_biased = [rR_biased stats.acc];
            
            
            % ----------- Correctness of an Decoding Instance -------------
            % ?????????? modification needed for corr. smoothing
            
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
        global h1; global h2;
        set(h1,'XData',1:length(rL),'YData',rL)
        set(h2,'XData',1:length(rR),'YData',rR)

        % Plot rL_ & rR_biased
        global h3; global h4;
        set(h3,'XData',1:length(rL_biased),'YData',rL_biased)
        set(h4,'XData',1:length(rR_biased),'YData',rR_biased)

        
        % display the latest corr. coef.
        figure(2)
        text(0,-1, ['rL: ' num2str(round(rL(end),2)), ...
            '  /  rR: ' num2str(round(rR(end),2)), ...
            '  /  corr:' num2str(acc(end)), ...
            '  /  corrAve: ' num2str(mean(acc))], 'HorizontalAlignment', 'center')
            
        text(0,-2, ['rL\_bias: ' num2str(round(rL_biased(end),2)), ...
            '  /  rR\_bias: ' num2str(round(rR_biased(end),2)), ...
            '  /  corr\_bias:' num2str(acc_biased(end)), ...
            '  /  corrAve\_bias: ' num2str(mean(acc_biased))], 'HorizontalAlignment', 'center')

        
        % ###################################################################
        % #                                                                 #
        % #                                                                 #
        % ######### A Place for LSL Outlet!: Send Corr. Coef. in NFB ########
        % #                                                                 #
        % #                                                                 #
        % ###################################################################
        
    end % end IsTest
    
    % ======================== After Trial Offset ============================
    
    if  exist('onTrial2use', 'var') && (onTrial2use >= 60) % when onTrial is over the trial length
        initGlobal_eeg_buffer   % initialize eeg_buffer
        
        clear global onTrial;   % terminate onTrial
        clear onTrial2use;      % terminate onTrial2use
        
        clear global f_onTrial; % terminate f_onTrial
        clear f_onTrial2use;    % terminate f_onTrial2use
        
        initGlobal_dir;         % initialize dir
        dir2use = 0;            % initialize dir2use
    end
    
end % end of data processing conditon-if




