%% Online Simulation
%% Subject
clear
clc

Mov = 1;
Width = 15;

fs = 64;
Dir    = -1;
tmin    = 0;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 250;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
lambda = 10;
model_cell = {};
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );

%filter
design = designfilt('bandpassfir', 'FilterOrder',601, ...
    'CutoffFrequency1', 0.5, 'CutoffFrequency2', 8,...
    'SampleRate', 125);
% fvtool(design);

Allspeech = readNPY('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\Allspeech.npy');
unatt_stim_cell = Allspeech(1:30,:);     % 1 by trials cell - 1 by time
att_stim_cell = Allspeech(31:end,:);   % 1 by trials cell - 1 by time
for SubIdx = 1:9
                
    load('C:\Users\LeeJiWon\Desktop\hykist\인수인계\Recording data_CEDA\'+string(SubIdx)+'\RAW_'+string(SubIdx)+'.mat');
    load('C:\Users\LeeJiWon\Desktop\hykist\인수인계\Recording data_CEDA\'+string(SubIdx)+'\AUX_'+string(SubIdx)+'.mat');

    % searching onset 
    temp = diff(TRIGGER(2,:));
    start = find(temp > 0)+1;
    % start = [1, start];
    try
        for i = length(start):-1:1
            if abs(start(i)-start(i-1)) == 1
                start(i)=[];
            end
        end
    end

    %preprocessing
    EEG = {};
    srate = 125;
    % 채널 하나 제거
    for j = 1:length(start)
        eeg = RAW([1:7,9:16],(start(j)+(srate*3)):(start(j)+(srate*63))-1); 
        EEG{j} = eeg';
    end

    eeg_cell = EEG;%RAWDATA_SPL(SubIdx).data;              % 1 by trials cell - time by channel

    m = zeros(15,length(t));
    b = 0;
    % Train model
    for tr = 1:14
        eeg = double(eeg_cell{tr});                % time by channel
        att_stim = att_stim_cell(tr,:);      % 1 by time
        unatt_stim = unatt_stim_cell(tr,:);  % 1 by time

        for w = 1:46
            Mov_eeg = eeg( (1+(Mov*(w-1))*srate) : (Width+(Mov*(w-1)))*srate,:);
            Mov_att_stim = att_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';

            % rereference
            for ch = 1:size(Mov_eeg,2)
                re = mean(Mov_eeg(:,ch));
                Mov_eeg(:,ch) = Mov_eeg(:,ch) - re;
            end
            % filter
                Mov_eeg = filtfilt(design, Mov_eeg);
%             eeg2use = filtfilt(filterweights, 1, Mov_eeg); % time_pnts by chans
            % downsample
            Mov_eeg = resample(Mov_eeg,fs,srate);
            % zscore
            Mov_eeg = zscore(Mov_eeg);            

            model = mTRFtrain(Mov_att_stim, Mov_eeg, fs, Dir, tmin, tmax, lambda, 'verbose', 0);

            m = m + model.w;
            b = b + model.b;
        end
    end
    disp(['Subject' num2str(SubIdx) ' Train end!']) 
    %%%%%%%%%%%%%%%%%%%%%%%
    % Averaging Decoder model
    model.w = m/(46*14);
    model.b = b/(46*14);

    % Test model
    for tr = 15:length(eeg_cell)        
        eeg = double(eeg_cell{tr});                % time by channel
        att_stim = att_stim_cell(tr,:);      % 1 by time
        unatt_stim = unatt_stim_cell(tr,:);  % 1 by time       

        for w = 1:46
            Mov_eeg = eeg( (1+(Mov*(w-1))*srate) : (Width+(Mov*(w-1)))*srate,:);
            Mov_att_stim = att_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';
            Mov_unatt_stim = unatt_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';

            % rereference
            for ch = 1:size(Mov_eeg,2)
                re = mean(Mov_eeg(:,ch));
                Mov_eeg(:,ch) = Mov_eeg(:,ch) - re;
            end
            Mov_eeg = filtfilt(design, Mov_eeg);
            % downsample
            Mov_eeg = resample(Mov_eeg,fs,srate);
            % zscore
            Mov_eeg = zscore(Mov_eeg); 

            [predict_att, stats_att] = mTRFpredict(Mov_att_stim, Mov_eeg, model, 'verbose', 0);
            [predict_unatt, stats_unatt] = mTRFpredict(Mov_unatt_stim, Mov_eeg, model, 'verbose', 0);

            if stats_att.r > stats_unatt.r
                acc(w) = 1;
            else acc(w) = 0; end
        end     % end 46 windows

        AccTrial(tr-14) = mean(acc);
    end

    Ave_Acc(SubIdx) = mean(AccTrial);
    disp(['Accuracy = ' num2str(mean(AccTrial))])
    disp(['Subject' num2str(SubIdx) ' finished!'])    
end



%% Subject - Bad group
clear
clc
trials_mcl = [15,22,24,27,33,38,42];
trials_20 = [17,20,25,28,32,36,41];
trials_90 = [16,19,23,30,34,37,40];
trials_srt = [18,21,26,29,31,35,39];
condition = {'mcl', '20', '90', 'srt'};
list = {'MCL', '20', '90', 'SRT'};

load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\DATA_withSNR.mat');
load ('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\data_withSNR\Rawdata_SPL_all.mat');
RAWDATA_SPL = RAWDATA_SPL_all;
Mov = 1;
Width = 15;
fs = 64;
Dir    = -1;
tmin    = 0;                       % min time-lag(ms)
tminIdx = floor(tmin/1000*fs);     % tmin2idx
tmax    = 250;                     % max time-lag(ms)
tmaxIdx = ceil(tmax/1000*fs);      % tmax2idx
lambda = 10;
model_cell = {};
timeLag = length(tminIdx:tmaxIdx); % time-lag points between tmin and tmax
t       = sort( 1 * ( linspace(tmin, tmax, timeLag) ) );

% filter
design = designfilt('bandpassfir', 'FilterOrder',601, ...
    'CutoffFrequency1', 0.5, 'CutoffFrequency2', 8,...
    'SampleRate', 125);
%fvtool(design);

tl=1;
% decoder
for SubIdx = 2:length(DATA_withSNR) % lde 빼고,,
    subject = DATA_withSNR(SubIdx).subject;

    load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\'+string(subject)+'\RAW_'+string(subject)+'.mat');
    load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Python\save_data\'+string(subject)+'\AUX_'+string(subject)+'.mat');

    % searching onset 
    temp = diff(AUX(2,:));
    start = find(temp > 0)+1;
    try
        for i = length(start):-1:1
            if abs(start(i)-start(i-1)) == 1
                start(i)=[];
            end
        end
    end
    if length(start) == 41
        start = [1, start];
    elseif length(start) < 5
        temp = diff(AUX(3,:));
        start = find(temp > 0)+1;
        try
            for i = length(start):-1:1
                if abs(start(i)-start(i-1)) == 1
                    start(i)=[];
                end
            end
        end
    end
    %preprocessing
    EEG = {};
    srate = 125;
    % 채널 하나 제거
    for j = 1:length(start)
        eeg = RAW([1:7,9:16],(start(j)+(srate*3)):(start(j)+(srate*63))-1); 
        EEG{j} = eeg';
    end

    eeg_cell = EEG;%RAWDATA_SPL(SubIdx).data;              % 1 by trials cell - time by channel
    att_stim_cell = RAWDATA_SPL(SubIdx).att_stim;     % 1 by trials cell - 1 by time
    unatt_stim_cell = RAWDATA_SPL(SubIdx).utt_stim;   % 1 by trials cell - 1 by time

    m = zeros(15,length(t));
    b = 0;
    % Train model
    for tr = 1:14
        eeg = double(eeg_cell{tr});                % time by channel
        att_stim = double(att_stim_cell{tr});      % 1 by time
        unatt_stim = double(unatt_stim_cell{tr});  % 1 by time

        for w = 1:46
            Mov_eeg = eeg( (1+(Mov*(w-1))*srate) : (Width+(Mov*(w-1)))*srate,:);
            Mov_att_stim = att_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';

            % rereference
            for ch = 1:size(Mov_eeg,2)
                re = mean(Mov_eeg(:,ch));
                Mov_eeg(:,ch) = Mov_eeg(:,ch) - re;
            end
            % filter
            Mov_eeg = filtfilt(design, Mov_eeg);
            % downsample
            Mov_eeg = resample(Mov_eeg,fs,srate);
            % zscore
            Mov_eeg = zscore(Mov_eeg);            

            model = mTRFtrain(Mov_att_stim, Mov_eeg, fs, Dir, tmin, tmax, lambda, 'verbose', 0);

            m = m + model.w;
            b = b + model.b;
        end
    end
    disp(['Subject' num2str(SubIdx) ' Train end!']) 
    %%%%%%%%%%%%%%%%%%%%%%%
    % Averaging Decoder model
    model.w = m/(46*14);
    model.b = b/(46*14);

    % Test model
    for tr = 15:length(eeg_cell)        
        eeg = double(eeg_cell{tr});                % time by channel
        att_stim = double(att_stim_cell{tr});      % 1 by time
        unatt_stim = double(unatt_stim_cell{tr});  % 1 by time        

        for w = 1:46
            Mov_eeg = eeg( (1+(Mov*(w-1))*srate) : (Width+(Mov*(w-1)))*srate,:);
            Mov_att_stim = att_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';
            Mov_unatt_stim = unatt_stim(1, (1+(Mov*(w-1))*fs) : (Width+(Mov*(w-1)))*fs)';

            % rereference
            for ch = 1:size(Mov_eeg,2)
                re = mean(Mov_eeg(:,ch));
                Mov_eeg(:,ch) = Mov_eeg(:,ch) - re;
            end
            Mov_eeg = filtfilt(design, Mov_eeg);
            % downsample
            Mov_eeg = resample(Mov_eeg,fs,srate);
            % zscore
            Mov_eeg = zscore(Mov_eeg); 

            [predict_att, stats_att] = mTRFpredict(Mov_att_stim, Mov_eeg, model, 'verbose', 0);
            [predict_unatt, stats_unatt] = mTRFpredict(Mov_unatt_stim, Mov_eeg, model, 'verbose', 0);

            if stats_att.r > stats_unatt.r
                acc(w) = 1;
            else acc(w) = 0; end
        end     % end 46 windows

        AccTrial(tr-14) = mean(acc);
    end

    Ave_Acc(SubIdx,tl) = mean(AccTrial);
    Acc_mcl(SubIdx,:,tl) = AccTrial([trials_mcl]-14);
    Ave_Acc_mcl(SubIdx,tl) = mean(Acc_mcl(SubIdx,:,tl));

    Acc_20(SubIdx,:,tl) = AccTrial([trials_20]-14);
    Ave_Acc_20(SubIdx,tl) = mean(Acc_20(SubIdx,:,tl));

    Acc_90(SubIdx,:,tl) = AccTrial([trials_90]-14);
    Ave_Acc_90(SubIdx,tl) = mean(Acc_90(SubIdx,:,tl));

    Acc_srt(SubIdx,:,tl) = AccTrial([trials_srt]-14);
    Ave_Acc_srt(SubIdx,tl) = mean(Acc_srt(SubIdx,:,tl));

    disp(['Subject' num2str(SubIdx) ' finished!'])    
end

%% 
for s = 2:length(DATA_withSNR)
    AveAcc_simulation(s).subject = DATA_withSNR(s).subject;
    AveAcc_simulation(s).all_acc = Ave_Acc(s)*100;
    for c =1:length(condition)
        eval(['AveAcc_simulation(s).acc_',condition{c},' = Ave_Acc_',condition{c},'(s)*100;']);
    end
end
    
    

%%

% acc-condi
Xa = categorical({'All','MCL','MCL-20', 'SI90', 'SRT'});
Xa = reordercats(Xa,{'All','MCL','MCL-20', 'SI90', 'SRT'});
Y = [mean(Ave_Acc,1); mean(Ave_Acc_mcl,1); mean(Ave_Acc_20,1); mean(Ave_Acc_90,1); mean(Ave_Acc_srt,1)]*100;

figure(24)
clf
b = bar([1,2,3,4,5], Y, 'FaceColor', 'flat');  hold on
plot([1-0.15, 1+0.15], Ave_Acc*100, '--ok')
for i = 1:4
    plot([(i+1)-0.15,(i+1)+0.15],eval(['Ave_Acc_',condition{i},'*100']), '--ok')
end


set(gcf, 'color', 'white')
ylim([0 100])
ylabel('Decoder Accuracy (%)')
% refline([0, chance]);
title('Online simulation - Mean Accuracy - Conditions')

%% 비교 bar 250/500

t250 = [[0.625000000000000,0.741459627329193,0.715838509316770,0.465838509316770,0.624223602484472,0.790372670807453,0.583850931677019,0.508540372670808]]*100;
t500 = [[0.651397515527950,0.765527950310559,0.778726708074534,0.482142857142857,0.609472049689441,0.758540372670807,0.628105590062112,0.518633540372671]]*100;       

figure(1)
x = categorical({'0-250', '0-500'});
x = reordercats(x,{'0-250', '0-500'});
y = [mean(t250); mean(t500)];
b = bar(x, y, 'FaceColor', 'flat');
b.CData(1,:) = [0.4660 0.6740 0.1880];
b.CData(2,:) = [0 0.4470 0.7410];
hold on
plot(x, [t250; t500], '--ok')
ylim([30 90])

%% Online exp vs simulation vs off
condition = {'mcl', '20', '90', 'srt'};
list = {'MCL', '20', '90', 'SRT'};

for sub = 1:length(DATA_withSNR) % lde 지운 후 / lji 도
    eval(['Ave_Acc(sub,1) = DATA_withSNR(sub).All_Acc;']);
    eval(['Ave_Acc(sub,2) = AveAcc_simulation(sub).all_acc;']);
    eval(['Ave_Acc(sub,3) = OffDATA(sub).all_acc;']);
    for c = 1:length(list)
        eval(['Ave_Acc_',condition{c},'(sub,1) = DATA_withSNR(sub).Acc_',list{c},';']);   
        eval(['Ave_Acc_',condition{c},'(sub,2) = AveAcc_simulation(sub).acc_',condition{c},';']); 
        eval(['Ave_Acc_',condition{c},'(sub,3) = OffDATA(sub).acc_',condition{c},';']);  
    end
end

% plot
% l = [12,16,17,18,20];
l = [1,6,8,9,12,14,16,17,19,20];
Y_exp = [Ave_Acc(l,1), Ave_Acc_mcl(l,1), Ave_Acc_20(l,1), Ave_Acc_90(l,1), Ave_Acc_srt(l,1)];
Y_simul = [Ave_Acc(l,2), Ave_Acc_mcl(l,2), Ave_Acc_20(l,2), Ave_Acc_90(l,2), Ave_Acc_srt(l,2)];
Y_off = [Ave_Acc(l,3), Ave_Acc_mcl(l,3), Ave_Acc_20(l,3), Ave_Acc_90(l,3), Ave_Acc_srt(l,3)];

figure(27)
clf
y = [mean(Y_simul); mean(Y_off)]';
ngroups = size(y, 1);
nbars = size(y, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5)); 
X = categorical({'All','MCL','MCL-20dB', 'SI 90%', 'SI 50%'});
X = reordercats(X,{'All','MCL','MCL-20dB', 'SI 90%', 'SI 50%'});

b = bar([1,2,3,4,5], y, 'FaceColor', 'flat');  hold on
for i = 1:5
    p = plot([i-0.15,i+0.15],[[Y_simul(:,i),Y_off(:,i)]'], '--o', 'Linewidth', 1.5);
%     p(1).Color = [0 0.4470 0.7410];      % blue
%     p(2).Color = [0.8500 0.3250 0.0980]; % red
%     p(3).Color = [0.9290 0.6940 0.1250]; % yellow
%     p(4).Color = [0.4940 0.1840 0.5560]; % purple
%     p(5).Color = [0.4660 0.6740 0.1880]; % green
%     legend()
end







    