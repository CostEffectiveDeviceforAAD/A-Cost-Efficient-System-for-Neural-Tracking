%% AADC_offline decoder

clear
%% load data_AAK
clear
%path = 'C:\Users\user\Desktop\hy-kist\OpenBCI\Recording data\'
path = 'C:\Users\LeeJiWon\Desktop\hykist\AAD\OnlineAAD\Recording data\'

sub = string(2)

file = strcat(path, sub, '\E_', sub, '.mat')
file2 = strcat(path, sub, '\A_', sub, '.mat')

% Resp load
load(file)      % EEG ;1by30 cell / timeby channel
load(file2)

% Stim load
load 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\Allspeech.mat' 
%load 'C:\Users\user\Desktop\hy-kist\OpenBCI\AAD\Matlab\Allspeech.mat'

%% Struct stim, resp 

stim_R = Allspeech(1:30,:);
stim_L = Allspeech(31:end,:);

% Stim cell array
StimC_L = cell(1,30);
for i=1:30
    sl = stim_L(i,:)';
    StimC_L{1,i} = sl;
end

StimC_R = cell(1,30);
for i=1:30
    sr = stim_R(i,:)';
    StimC_R{1,i} = sr;
end

%% Resp preprocessing

%path = 'C:\Users\user\Desktop\hy-kist\OpenBCI\Recording data\'
path = 'C:\Users\LeeJiWon\Desktop\hykist\AAD\OnlineAAD\Recording data\'

% Stim load
load 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\Allspeech.mat' 
%load 'C:\Users\user\Desktop\hy-kist\OpenBCI\AAD\Matlab\Allspeech.mat'
i=1
for s = [2,5,6,7,9,10,11]
    
    sub_eeg = [];
    
    sub = num2str(s);

    file = strcat(path, sub, '\E_', sub, '.mat');
    file2 = strcat(path, sub, '\A_', sub, '.mat');

    % Resp load
    load(file);      % EEG ;1by30 cell / timeby channel
    load(file2);

    for tr = 1:30
        raw = EEG{tr};
        aux = AUX{tr};

        % Delete Fp1
        raw(:,8) =[];

        % find stimuli onset
        tri = diff(aux(:,2));
        index = find(tri>0);
        speech = index(end)+(125*3);

        eeg = raw(speech+1:speech+(125*60),:);

        % rereference
        for t = 1:15
            m = mean(eeg(:,t));
            eeg(:,t) = eeg(:,t) - m;
        end

        % filter
        design = designfilt('bandpassfir', 'FilterOrder',601, ...
            'CutoffFrequency1', 0.5, 'CutoffFrequency2', 8,...
            'SampleRate', 125);

        eeg = filtfilt(design, eeg);

        % downsample
        eeg = resample(eeg,64,125);

        % zscore
        eeg = zscore(eeg);

        sub_eeg = cat(3, sub_eeg, eeg);
    end
    RespC{i} = eeg;
    DATA{i} = sub_eeg;
    i = i+1
    
end

%% Searching Optimal lambda

%------------------------
% Model parameters
Dir = -1;             
tmin = 0;         
tmax = 250;            
%lambda = 2.^(0:2:20); 
lambda = 10.^(-6:1:4)
fs = 64              


%%
% Atteded decoder
% Search optimal lambda
cv = mTRFcrossval(StimC_L,RespCC,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',1);
[rmax,idx] = max(mean(cv.r));
optmal = lambda(idx)

%% Online


for idx = 1:11
    
    model_tr = zeros(15,17);
    bias_tr = 0;

    % Train
    for tr = 1:14   % trials

        eeg = RespC{tr};
        model_w = zeros(15,17);
        bias_w = 0;

        for i = 1:46      % window number

            resp = eeg(fs*(i-1)+1 : fs*(i+14), :);
            stim_a = stim_L(tr, fs*(i-1)+1 : fs*(i+14))'; 

            % Preprocessing
            % rereference
            for t = 1:15
                m = mean(resp(:,t));
                resp(:,t) = resp(:,t) - m;
            end
            % filter
            resp = filtfilt(design, resp');
            % downsample
            resp = resample(resp,64,125);
            % zscore
            resp = zscore(resp);

            % Train
            model = mTRFtrain(stim_a,resp,fs,Dir,tmin,tmax,10,'zeropad',1); 

            % sum model
            model_w = model_w + model.w;
            bias_w = bias_w + model.b; 

        end     %  a trial end

        model_tr = model_tr + model_w;
        bias_tr = bias_tr + bias_w;

    end    

    % Average model
    model.w = model_tr/(14*46)
    model.b = bias_tr/(14*46)


    % Test
    for tr = 15:30   % trials

        eeg = RespC{tr};

        for i = 1:46      % window number

            resp = eeg(125*(i-1)+1 : 125*(i+14), :);
            stim_a = stim_L(tr, fs*(i-1)+1 : fs*(i+14))';
            stim_u = stim_R(tr, fs*(i-1)+1 : fs*(i+14))';

            % Preprocessing      
            % rereference
            for t = 1:15
                m = mean(resp(:,t));
                resp(:,t) = resp(:,t) - m;
            end
            % filt
            resp = filtfilt(design, resp);
            % downsample
            resp = resample(resp,64,125);
            % zscore
            resp = zscore(resp);

            % reconstruction
            pred_at = mTRFpredict(stim_a, resp, model);     % 960 by 1

            % Evaluate
            r_at(i) = mTRFevaluate(stim_a, pred_at, 'corr', 'pearson', 'error', 'mse');
            r_un(i) = mTRFevaluate(stim_u, pred_at, 'corr', 'pearson', 'error', 'mse'); 

            if r_at(i) > r_un(i)
                acc(i) = 1;
            else
                acc(i) = 0;
            end

        end

        Acc(tr) = mean(acc)*100;

    end     %  a trial end

    % Estimate Accuracy
    Accuracy(idx) = mean(Acc(15:end));   
    acc_fix(idx) = mean(Acc(15:26));
    acc_swi(idx) = mean(Acc(27:end));

end
 
[Acc_max, Acc_idx] = max(Accuracy);
Acc_fix = acc_fix(Acc_idx);
Acc_swi = acc_swi(Acc_idx);
Acc_lambda = lambda(Acc_idx);

