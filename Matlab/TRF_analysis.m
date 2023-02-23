%% Temporal Response Functions (TRF)
%% DATA format - 0.5 ~ 25 Hz
% Data analysis
clear
clc
% load
subject = '0928_kms';
load('C:\Users\LeeJiWon\Desktop\hykist\AAD\Speech\2022\2022_Speech\Allsource_snr.mat');
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\'+string(subject)+'\RAW_'+string(subject)+'.mat');
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\'+string(subject)+'\AUX_'+string(subject)+'.mat');

% searching onset 
temp = diff(AUX(2,:));
start = find(temp > 0)+1;
% start = [1, start];
try
    for i = length(start):-1:1

        if abs(start(i)-start(i-1)) == 1

            start(i)=[];
        end
    end
end

% preprocessing
EEG = {};
srate = 125;

% filter
design = designfilt('bandpassfir', 'FilterOrder',701, ...
    'CutoffFrequency1', 0.5, 'CutoffFrequency2', 25,...
    'SampleRate', 125);
% fvtool(design);
    
for j = 1:length(start)
    eeg = RAW([1:7,9:16],(start(j)+(srate*3)):(start(j)+(srate*63))-1); 

    % rereference
    for ch = 1:size(eeg,1)
        m = mean(eeg(ch,:));
        eeg(ch,:) = eeg(ch,:) - m;
    end
    
    eeg = filtfilt(design, eeg');

    % downsample
    eeg = resample(eeg,64,srate);

    % zscore
    eeg = zscore(eeg);

    EEG{j} = eeg;
end

% speech envelope format
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\data\SAT_envelope.mat');  % sound env
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\data\AAK_envelope.mat'); 
att_env={};
utt_env={};
for i = 1:length(EEG)
    if i < 15
        att_env{i} = Allspeech_JT(i,:);
        utt_env{i} = Allspeech_JT(i+30,:);
    else              
        att_env{i} = Allspeech_SAT(i-14,:);
        utt_env{i} = Allspeech_SAT(i+16,:);
    end 
end

% all data format
load ('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\Rawdata_TRF_SPL.mat') % RAWDATA_SPL
l = length(RAWDATA_TRF_SPL)+1;
RAWDATA_TRF_SPL(l).subject = subject;
RAWDATA_TRF_SPL(l).data = EEG;
RAWDATA_TRF_SPL(l).att_stim = att_env;
RAWDATA_TRF_SPL(l).utt_stim = utt_env;

% save('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\Rawdata_TRF_SPL.mat', 'RAWDATA_TRF_SPL')

%%
% ddd = EEG;
% att = att_env;
% utt = utt_env;

ori = 28+14
pre = 28+14;

EEG(ori) = ddd(pre);
att_env(ori) = att(pre);
utt_env(ori) = utt(pre);

%% Offline - TRF - Encoding
clear 
clc
load ('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\Rawdata_TRF_SPL.mat');
% load ('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\Rawdata_SPL.mat');
% RAWDATA_TRF_SPL = RAWDATA_SPL;

trials_mcl = [15,22,24,27,33,38,42];
trials_20 = [17,20,25,28,32,36,41];
trials_90 = [16,19,23,30,34,37,40];
trials_srt = [18,21,26,29,31,35,39];
condition = {'mcl', '20', '90', 'srt'};
list = {'MCL', '20', '90', 'SRT'};
front = [1,2,3,4,8,9,10,11];
% front = [1:15];

fs = 64;
Dir    = 1;
tmin    = -150;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 500;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );
lambda = 10^4;
weight = {};

for c = 1:length(condition)
    eval(['Ave_weight_',condition{c},' = zeros(timeLag, length(front));'])   % timelag index by channel
end
Ave_weight = zeros(timeLag, length(front));

% design = designfilt('bandpassfir', 'FilterOrder',601, ...
%     'CutoffFrequency1', 2, 'CutoffFrequency2', 8,...
%     'SampleRate', 125);

for SubIdx = 1:length(RAWDATA_TRF_SPL)
    
    att_env = RAWDATA_TRF_SPL(SubIdx).att_stim;
    utt_env = RAWDATA_TRF_SPL(SubIdx).utt_stim;
    EEG = RAWDATA_TRF_SPL(SubIdx).data;
    
    % reject channel
    for i = 1:length(EEG)
        EEG{i} = EEG{i}(:,front);
    end   
    model = {};
    
    % freq band 변경
%     for tr = 1:length(EEG)
%         eeg = EEG{tr};
%         EEG{tr} = filtfilt(design, eeg);
%     end
    
    % encoder train
    model = mTRFtrain(att_env(1:14), EEG(1:14), fs, Dir, tmin, tmax, lambda, 'verbose', 0);
    weight{SubIdx} = squeeze(model.w);
    
    % train model for trf
    for c = 1:length(condition)
        eval(['model_trf = mTRFtrain(att_env(trials_',condition{c},'), EEG(trials_',condition{c},'), fs, Dir, tmin, tmax, lambda, ''verbose'', 0);'])    
        eval(['weight_',condition{c},'{SubIdx} = squeeze(model_trf.w);'])
    end
 
    for li = 1:length(list)
        eval(['Offcorr_att_',list{li},'= [];'])
        eval(['Offcorr_utt_',list{li},'= [];']); 
        eval(['Pred_att_',condition{li},'= {};'])
        eval(['Pred_utt_',condition{li},'= {};']);       
    end

    for i = 15:length(EEG)
%         j = i+14;
        if length(find(trials_mcl == i)) == 1;
            [pred_att_mcl, stats_att_mcl] = mTRFpredict(att_env{i}, EEG{i}, model, 'verbose', 0);
            [pred_utt_mcl, stats_utt_mcl] = mTRFpredict(utt_env{i}, EEG{i}, model, 'verbose', 0);
            
            Pred_att_mcl = [Pred_att_mcl, pred_att_mcl];
            Pred_utt_mcl = [Pred_utt_mcl, pred_utt_mcl];
            Offcorr_att_MCL = [Offcorr_att_MCL; stats_att_mcl.r];
            Offcorr_utt_MCL = [Offcorr_utt_MCL; stats_utt_mcl.r];

        elseif length(find(trials_20 == i)) == 1;
            [pred_att_20, stats_att_20] = mTRFpredict(att_env{i}, EEG{i}, model, 'verbose', 0);
            [pred_utt_20, stats_utt_20] = mTRFpredict(utt_env{i}, EEG{i}, model, 'verbose', 0);
            
            Pred_att_20 = [Pred_att_20, pred_att_20];
            Pred_utt_20 = [Pred_utt_20, pred_utt_20];
            Offcorr_att_20 = [Offcorr_att_20; stats_att_20.r];
            Offcorr_utt_20 = [Offcorr_utt_20; stats_utt_20.r];

        elseif length(find(trials_90 == i)) == 1;
            [pred_att_90, stats_att_90] = mTRFpredict(att_env{i}, EEG{i}, model, 'verbose', 0);
            [pred_utt_90, stats_utt_90] = mTRFpredict(utt_env{i}, EEG{i}, model, 'verbose', 0);
            
            Pred_att_90 = [Pred_att_90, pred_att_90];
            Pred_utt_90 = [Pred_utt_90, pred_utt_90];
            Offcorr_att_90 = [Offcorr_att_90; stats_att_90.r];
            Offcorr_utt_90 = [Offcorr_utt_90; stats_utt_90.r];

        elseif length(find(trials_srt == i)) == 1;
            [pred_att_srt, stats_att_srt] = mTRFpredict(att_env{i}, EEG{i}, model, 'verbose', 0);
            [pred_utt_srt, stats_utt_srt] = mTRFpredict(utt_env{i}, EEG{i}, model, 'verbose', 0);
            
            Pred_att_srt = [Pred_att_srt, pred_att_srt];
            Pred_utt_srt = [Pred_utt_srt, pred_utt_srt];
            Offcorr_att_SRT = [Offcorr_att_SRT; stats_att_srt.r];
            Offcorr_utt_SRT = [Offcorr_utt_SRT; stats_utt_srt.r];
        end                
    end 
    
    % weight 정리
    Ave_weight = Ave_weight + weight{SubIdx};
    for c = 1:length(condition)
        eval(['Ave_weight_',condition{c},' = Ave_weight_',condition{c},' + weight_',condition{c},'{SubIdx};'])
    end

    %
    TRF_DATA(SubIdx).subject = RAWDATA_TRF_SPL(SubIdx).subject;
    for con = 1:length(condition)
        eval(['TRF_DATA(SubIdx).pred_att_',condition{con},' = Pred_att_',condition{con},';'])
        eval(['TRF_DATA(SubIdx).pred_utt_',condition{con},' = Pred_utt_',condition{con},';'])
        
        eval(['TRF_DATA(SubIdx).corr_att_',condition{con},' = Offcorr_att_',list{con},';'])
        eval(['TRF_DATA(SubIdx).corr_utt_',condition{con},' = Offcorr_utt_',list{con},';'])
    end
%     TRF_DATA_2(SubIdx).subject = RAWDATA_TRF_SPL(SubIdx).subject;
%     TRF_DATA_2(SubIdx).all_acc = mean([OffAcc_mcl,OffAcc_20, OffAcc_90,OffAcc_srt])*100;
% 
%     for i = 1:4
% %         eval(['OffDATA(SubIdx).corr_att_',condition{i},'= Offcorr_att_',list{i},';']);
%         eval(['OffDATA_2(SubIdx).corr_att_',condition{i},'= Offcorr_att_',list{i},';']);
%     end
%     for i = 1:4
% %         eval(['OffDATA(SubIdx).corr_utt_',condition{i},'= Offcorr_utt_',list{i},';']);
%         eval(['OffDATA_2(SubIdx).corr_utt_',condition{i},'= Offcorr_utt_',list{i},';']);
%     end
    
    disp(['subject',num2str(SubIdx),' finished!'])
end

Ave_weight = Ave_weight/SubIdx;
for c = 1:length(condition)
    eval(['Ave_weight_',condition{c},'= Ave_weight_',condition{c},'/SubIdx;'])
end
    
% save('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\TRF_DATA.mat', 'TRF_DATA')

%% Encoding plot
% load ('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\TRF_DATA.mat');
ch_name = {'Fz'; 'Cz'; 'C3'; 'C4'; 'P7'; 'P8'; 'Pz'; 'F7'; 'F8'; 'F3'; 'F4'; 'T7'; 'T8'; 'P3'; 'P4'};
condition = {'mcl', '20', '90', 'srt'};
list = {'MCL', '20', '90', 'SRT'};


trf_eeg = {};
ttemp = zeros(size(TRF_DATA(1).pred_att_mcl{1},1), size(TRF_DATA(1).pred_att_mcl{1},2));
for SubIdx = 1:length(TRF_DATA)
    pred_att_mcl = TRF_DATA(SubIdx).pred_att_20;   % 1 by 7 cell
    temp = zeros(size(pred_att_mcl{1},1), size(pred_att_mcl{1},2));

    for tr = 1:7
        temp = temp + pred_att_mcl{tr}; % condition 합치기 / 60*64 by channel
    end

    ttemp = ttemp + (temp/7);   % 피험자별 합치기
    
end
% 역전된 뇌파 바꿔주기
for i = [5,6,7,14,15]   % temporal 들
    ttemp(:,i) = -ttemp(:,i);
end

trf_eeg_mcl = ttemp/(length(TRF_DATA));
trf_one_mcl = mean(trf_eeg_mcl,2);

%
sec = 0.5;
fs = 64;

figure(1)
clf
x = linspace(-0.1,sec*1000,fs*sec);
[a,inx] = max(trf_one_mcl)
% peak = x(inx)
% plot(x, -trf_eeg_mcl(1:fs*sec,:))
plot(x, -trf_one_mcl(1:fs*sec,:))
% ylim([-0.04 0.06])
% ylim([-0.01 0.04])
title('SI 50%')

figure(2)
clf
eeg = ttemp/length(TRF_DATA);
for ch = 1 : size(pred_att_mcl{1},2)
    subplot(4,4,ch)
    plot(x, -trf_eeg_mcl(1:fs*sec,ch))
%     ylim([-0.04 0.06])
    title(ch_name{ch})
end

figure(3)
clf
plot(x, -trf_eeg_mcl(1:fs*sec,:))

%% TRF plotting
ch_name = {'Fz'; 'Cz'; 'C3'; 'C4'; 'P7'; 'P8'; 'Pz'; 'F7'; 'F8'; 'F3'; 'F4'; 'T7'; 'T8'; 'P3'; 'P4'};
condition = {'mcl', '20', '90', 'srt'};
list = {'MCL', '20', '90', 'SRT'};

trf_all_one = [];
for c = 1:length(condition)
    eval(['trf = Ave_weight_',condition{c},';'])
    % 역전된 뇌파 바꿔주기
%     for i = [5,6,7,12,13,14,15]   % p7,8,z,3,4, t7,8
%         trf(:,i) = -trf(:,i);
%     end

    trf_one = mean(trf, 2);
    trf_all_one(:,c) = trf_one; 

    figure(c)
    clf
    plot(t, trf); hold on
    title(condition{c})
    ylim([-0.02 0.02])
    xline(0, '--k')
%     legend(ch_name)
    saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\Ave\TRF_allchan_' +string(condition{c})+ '.jpg')

    figure(c+length(condition))
    clf
    plot(t, trf_one, 'Linewidth', 2)
    title(condition{c})
    ylim([-0.05 0.05])
    xline(0, '--k')
    saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\Ave\TRF_mean_' +string(condition{c})+ '.jpg')
    %peak 찾기
    [a,inx] = max(trf_one(16:30));
    peak(c) = t(inx+16)   
end

% train trials 들 결과값
trf = Ave_weight;
% 역전된 뇌파 바꿔주기
% for i = [5,6,7,12,13,14,15]   % p7,8,z,3,4, t7,8
%     trf(:,i) = -trf(:,i);
% end

trf_one = mean(trf, 2);
trf_all_one(:,5) = trf_one;

figure(9)
clf
plot(t, trf)
title('Train trials')
ylim([-0.25 0.25])
xline(0, '--k')
%     legend(ch_name)
saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\Ave\TRF_allchan_Train trials.jpg')

figure(10)
clf
plot(t, trf_one, 'Linewidth', 2)
title('Train trials')
ylim([-0.08 0.08])
xline(0, '--k')
saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\Ave\TRF_mean_Train trials.jpg')

% peak 찾기
[a,inx] = max(trf_one(16:30));
peak(c+1) = t(inx+16)   

% 모두 다같이
figure(11)
clf
plot(t, trf_all_one, 'Linewidth' ,2);
legend('MCL', 'MCL-20dBA', 'SI 90%', 'SI 50%', 'Train', 'AutoUpdate','off')
xline(0, '--k')
saveas(gcf, 'C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\Ave\TRF_All conditions.jpg')



%% TRF - previouse data
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\DATA(3conditions)_rev210330.mat');
load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\AllEEGdata_AAK_25.mat')
% load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\AllEEG_AAK.mat') % AllEEG_AAK 

% channel names
chanlocs = DATA(1).EEG.chanlocs;
chanlocs = struct2cell(chanlocs);
chanlocs = squeeze(chanlocs(1,1,:));
RoI = {'F5','F3','FZ', 'F4', 'F6','FC5', 'FC3', 'FCZ', 'FC4', 'FC6'};
channum = [];
for r = 1:length(RoI)    
    channum(r) = find(strcmp(chanlocs, RoI{r}));    
end
% 
% data = AllEEG_AAK;   % struct-subject / trial by channel by time
% data = AllEEGdata_AAk
data = DATA.EEG;
Stim = DATA.SPEECH;      % trial by time
%
att_env = {};
utt_env = {};
for i = 1:30
    att_env{i} = Stim(i,:)';
    utt_env{i} = Stim(i+30,:)';
end

fs = 64;
Dir = 1;
% Time-lag (tau) params
tmin    = -125;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 500;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );
lambda = 10^3;

%filter
% design = designfilt('bandpassfir', 'FilterOrder',601, ...
%     'CutoffFrequency1', 2, 'CutoffFrequency2', 8,...
%     'SampleRate', 125);

Ave_weight = zeros(length(channum), length(t));
for SubIdx = 1:length(data)
    EEG = {};
    model = {};
    
    % cell 만들기
    for tr = 1:30
%         EEG{tr} = filtfilt(design, double(squeeze(data(SubIdx).data(tr,channum,:))'));    % trial by channel by time   
        EEG{tr} = double(squeeze(data(SubIdx).data(tr,channum,:))');    % trial by channel by time
    end
  
    % encoder train
    model = mTRFtrain(att_env, EEG, fs, Dir, tmin, tmax, lambda, 'verbose', 0);
    weight{SubIdx} = squeeze(model.w);
    
    Ave_weight = Ave_weight + weight{SubIdx}';
    disp(['subject',num2str(SubIdx),' finished!'])
end

Ave_weight = Ave_weight/SubIdx;
one_weight = mean(Ave_weight,1);
% plot
figure(1)
clf
plot(t, Ave_weight);

figure(2)
clf
plot(t, one_weight, 'Linewidth' ,2);
xline(0, '--k')










