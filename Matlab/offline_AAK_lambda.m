%% AADC_offline decoder

clear
%% load data_AAK

sub = '0806_LJH'
file = strcat('C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\', sub, '\E_', sub, '.mat')
file2 = strcat('C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\', sub, '\A_', sub, '.mat')
% Stim & Resp load
load(file)      % EEG ;1by30 cell / timeby channel
load(file2)
load 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\Allspeech.mat' 

%% Decoder

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

%------------------------------------------
%% eeg restruct
for i = 1:30
    eeg = EEG{i};
    aux = AUX{i};

    tri = diff(aux(:,2));
    index = find(tri>0);
    speech = index+(125*3)+1;

    eeg = eeg(speech:speech+(125*60)-1,:);

    % rereference
    for t = 1:16
        m = mean(eeg(t, :));
        eeg(t,:) = eeg(t, :) - m;
    end

    % filter
    design = designfilt('bandpassfir', 'FilterOrder',601, ...
        'CutoffFrequency1', 0.5, 'CutoffFrequency2', 30,...
        'SampleRate', 125);

    eeg = filtfilt(design, eeg);

    % downsample
    eeg = resample(eeg,64,125);

    % zscore
    eeg = zscore(eeg);

    RespC{i} = eeg;
end

%%

%------------------------
% Model parameters
Dir = -1;             
tmin = 0;         
tmax = 250;            
%lambda = 2.^(0:2:20); 
lambda = 10.^(-6:2:6)
fs = 64              


% Atteded decoder
% Search optimal lambda
cv = mTRFcrossval(StimC_L,RespC,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',1);
[rmax,idx] = max(mean(cv.r));
index = idx;

% Train
model = mTRFtrain(StimC_L,RespC,fs,Dir,tmin,tmax,lambda(idx),'zeropad',1);  %lambda(idx)

% Predicted envelope
pred_at = mTRFpredict(StimC_L,RespC,model);

% Estimate correlatioin coefficient
pred = pred_at';

for i = 1:30
    yl = stim_L(i,:)';
    yr = stim_R(i,:)';
    pi = pred{i};

    r_L(i) = mTRFevaluate(yl,pi, 'corr', 'pearson', 'error', 'mse');
    r_R(i) = mTRFevaluate(yr,pi, 'corr', 'pearson', 'error', 'mse');
        
    % estimate accuracy

    if r_L(i) > r_R(i)

        acc(i) = 1;

    else acc(i) = 0;

    end
    r_Unat(1,i) = r_R(i);
    r_Att(1,i) = r_L(i);
end
    
Acc_at = mean(acc);




%% plot
% Entire accuracy of each subject

figure
plot(Acc_at*100, '-ob', 'LineWidth', 1); %hold on 
%plot(Acc_un*100, '--^r'); hold off
mu = mean(Acc_at*100);
hline = refline([0 mu]);
hline.Color = 'b';
hline.LineWidth = 0.5;
set(gcf, 'color', 'white')
set(gca, 'color', 'white')
axis([1 10 0 100]); grid on
xlabel('Subject')
ylabel('Accuracy (%)')
title('Accuracy - Decording_ AAK(19chan))','FontSize',10)
%legend('Attended decorder','Unatteded decorder')
 
    
    

