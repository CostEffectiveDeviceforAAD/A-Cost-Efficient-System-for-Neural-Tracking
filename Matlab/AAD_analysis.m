%% Real-time AAD analysis
%clear

sub = '_0812_LKS'
subject = 7

% sub = '_0811_JJH'
% subject = 5

%%
all = mean(Acc);
fix = mean(Acc(1:12));
swi = mean(Acc(13:end));

%% Behavior
file = strcat('Behavior',sub, '.mat');
load(file)

beh_at = Behavior(:,3:4);
beh_un = Behavior(:,1:2);

correct_at = find(beh_at=='T');
correct_un = find(beh_un=='T');

accb_at = (length(correct_at)/60)*100;
accb_un = (length(correct_un)/60)*100;  

%% Bar plot - behavior

X = categorical({'Attended','Unattended'});

for sub = 1:size(C,1)
    be_at(sub) = table2array(C(sub, 5));
    be_un(sub) = table2array(C(sub, 6));
end

y = [be_at; be_un];
Y = [mean(be_at), mean(be_un)];

% plot

b = bar(X, Y); hold on
plot(X, y, '--ok');
grid on
ylim([0 100])
ylabel('Accuracy(%)')
set(gcf, 'color', 'white')
title('Behavior Result')

%% Accuracy Plot
%% Accuracy

file = strcat('Accuracy',sub, '.mat');
load(file)  % Acc

data = readtable('C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\result.xlsx');

overall = mean(Acc*100);
fixed = mean(Acc(1:12)*100);
switching = mean(Acc(13:16)*100);

%data = table(subject, overall, fixed, switching, accb_at, accb_un)  % first suvject
data = [data; table(subject, overall, fixed, switching, accb_at, accb_un)]

chance = 52.99

%% Write to Excel file 
writetable(data, 'C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\result.xlsx');

%% Read previous AAD result

pre_data = readtable('C:\Users\LeeJiWon\Desktop\OpenBCI\OnlineAAD_Performance.xlsx');


%% subject_total
%C = readtable('result_sample.xlsx');

X = reordercats(X,{'Overall','Fixed','Switching'});

y = table2array(C(:, 2:4));
Y = [mean(table2array(C(:,2))), mean(table2array(C(:,3))),...
    mean(table2array(C(:,4)))];


b = bar(X, Y);  hold on
plot(X, y, '--ok');
grid on
set(gcf, 'color', 'white')
ylim([50 80])
ylabel('Accuracy(%)')
% refline([0, chance]);
title('Total')

%% multi subject

%C = readtable('result_sample.xlsx');
Y = []
X = []

for sub = 1:length(data.subject)
    name = strcat('Sub', num2str(sub));
    X = [X, categorical({name})];
    Y = [Y; table2array(data(sub,2:4))];
    
end

b = bar(X, Y);
yline(52.99, '-.k')
hline.Color = 'r';
grid on
legend('overall','fixed','switching');
ylim([0 80])
set(gcf, 'color', 'white')
ylabel('Accuracy(%)')
xlabel('Subject')

%%
clear

load 'Accuracy_chloc.mat' 
train_trial_ori = mean(Acc(1:14));
test_ori = mean(Acc(15:30));
train_mean_ori = mean(Acc(31:end));

%% channel search
clear
%load 'C:\Users\user\Desktop\hy-kist\OpenBCI\save_data\Accuracy_chloc.mat' 
load 'C:\Users\LeeJiWon\Desktop\OpenBCI\save_data\Accuracy_chloc.mat' 

for i = 1:size(Acc,1)
    check_Acc(i) = mean(Acc(i,1:14))*100;
    test_Acc(i)  = mean(Acc(i,15:30))*100;
end

maxval = max(test_Acc);
index = find(test_Acc==maxval);

%% channel graph

%ch = ['T7' 'Fz' 'T8' 'P4' 'P7' 'C3' 'O2'...
%                'Fp1' 'O1' 'Pz' 'F3' 'Cz' 'C4' 'P3' 'P8' 'F4'];                   
            
%ac = [61.41, 65.76, 66.85, 67.39, 67.39, 66.17, 65.76...
%        63.72, 62.77, 61.96, 60.60, 59.78, 58.42, 56.79, 54.89, 54.62];
     
%ac = [71.47, 76.90, 76.22, 75.00, 72.83, 72.28, 70.38, 68.34, 67.39...
%        68.75, 69.43, 70.38, 71.88, 69.57]  % bsc

%ac = [62.36, 69.16, 69.57, 71.33, 73.10, 74.59, 76.36, 76.36, 76.77, 76.49, 76.77, 76.90, 75.54, 74.18, 71.60];     %kkm

ac = [68.75, 73.37, 74.05, 72.55, 75.68, 75.82, 76.36,...
       76.22, 75.27, 74.86, 73.64, 72.69, 73.37, 72.01, 64.95]


x = 1:15 

% plot
figure;
plot(x, ac,'-or');
set(gcf, 'color', 'white')
set(gca, 'xtick', [1:15])
%set(gca, 'xticklabel', char('T7', 'Fz', 'T8', 'P4', 'P7', 'C3', 'O2', 'Fp1'...
%                                'O1', 'Pz', 'F3', 'Cz', 'C4', 'P3', 'P8', 'F4'))   %ljw
%set(gca, 'xticklabel', char('Fp1', 'P7', 'T7', 'F4', 'Pz/P3', 'F7', 'P4',...
%                                'O2', 'Fz', 'C4', 'C3/P8', 'F3', 'T8', 'Cz'))       % bsc
%set(gca, 'xticklabel', char('T8', 'Cz', 'Fz','F3', 'F7', 'C4', 'C3','P3', 'Pz',...
%                                'P7', 'P4', 'P8', 'F4', 'F8', 'T7'))       % kkm                                
set(gca, 'xticklabel', char('P8', 'C3', 'F4','C4', 'Fz', 'P7','P3', 'T8', 'Cz', 'P4',...
                                'Pz', 'T7', 'F7', 'F3', 'F8'))       % ctm  

ylim([60 80])
grid on 
xlabel('Channel')
ylabel('Accuracy(%)')


%% Envelope
load 'Predict_L.mat'
load 'Allspeech.mat'
i = 20;
L = 960;
x = linspace(0, L/64, L);
%%
for i = 0
    figure()
    subplot(311)
    plot(x, Allspeech(31+i,1:960)')
    title('Speech-att')
    subplot(312)
    plot(x, Allspeech(1+i,1:960)')
    title('Speech-unatt')
    subplot(313)
    plot(x, pre(1,:)')
    title('Predicted')
end

%% 비교-envelope
tr = 24
pre_l = squeeze(Pre_L(tr,:,:));


for i = 16:20
    figure()
    subplot(211)
    plot(x, Allspeech(31+i, 64*(i)+1:64*(15+i))'); hold on
    plot(x, pre_l(i+1,:)*5')
    title('Speech-att')
    subplot(212)
    plot(x, Allspeech(1+i, 64*(i)+1:64*(15+i))'); hold on
    plot(x, pre_l(i+1,:)*5')
    title('Speech-unatt')
end


%% 겹치기-att
for (i = 9) %&& (i = 24:25)
    figure()
    plot(x, Allspeech(31+i,1:64*(15))', 'r'); hold on
    plot(x, (Pre_L{1+i}(1,:))', 'b')
    legend('Speech', 'Predicted')
    title(strcat('Trial - ', num2str(i)))
end
%%
for (i = 16:20)
    figure()
    plot(x, Allspeech(1+i,1:960)', 'r'); hold on
    plot(x, Pre_L{1+i}(1,:)', 'b')
    legend('speech', 'predicted')
end
%16,17,18,24,25,

%% bar and spot
Y=[]
X=[]
X = categorical({'Overall','Fixed','Switching'});
X = reordercats(X,{'Overall','Fixed','Switching'});

for sub = 1: length(data.subject)
    Y = [Y; table2array(data(sub,2:end))];
end

Ym = mean(Y(:,1:3), 1)

b = bar(X, Ym);  hold on
plot(X, Y(:,1:3)', '--o');
grid on
ylim([0 100])
ylabel('Accuracy(%)')
% refline([0, chance]);
title('Total')

%% bar and spot - comparision
Y=[];
X = categorical({'Overall','Fixed','Switching'});
X = reordercats(X,{'Overall','Fixed','Switching'});

for sub = 1: size(data,1)
    Y = [Y; table2array(data(sub,2:4))];
end

% previour result
pre_Y = table2array(pre_data(:,2:4));

Ym = [mean(Y, 1); mean(pre_Y, 1)];

b = bar(X, Ym);  hold on
%plot(X, Y', '-o');
grid on
ylim([0 100])
ylabel('Accuracy(%)')
set(gcf, 'color', 'white')
legend('OpenBCI', 'Neuroscan');
title('Comparison')

%% result plot

% real-time three
X = categorical({'Overall','Fixed','Switching'});
X = reordercats(X,{'Overall','Fixed','Switching'});

all_sub = [67.39 67.53 70.11 88.32 63.18 60.05 73.51];
all = 70.01;
fix_sub = [ 68.12 69.02 72.83 91.85 63.59 62.50 76.09];
fix = 72.00
swi_sub = [ 65.22 63.04 61.96 77.72 61.96 52.72 65.76];
swi = 64.05

barp = [all; fix; swi];
barp_x = [all_sub; fix_sub; swi_sub];
%% 
figure
b = bar(X, barp, 'FaceColor', 'flat'); hold on
b.CData(1,:) = [0 0.4470 0.7410]
b.CData(2,:) = [0.8500 0.3250 0.0980]
b.CData(3,:) = [0.9290 0.6940 0.1250]
plot(X, barp_x, '--ok', 'MarkerFaceColor', [1 1 1], 'LineWidth', 0.8);
ylim([0 100])
ylabel('Accuracy(%)')
set(gcf, 'color', 'white')
title('Real-time accuracy')
%box('off')

%% conditions
X = categorical({'sub02','sub05','sub06','sub07','sub09','sub10','sub11'});
X = reordercats(X,{'sub02','sub05','sub06','sub07','sub09','sub10','sub11'});

real = [67.39 67.53 70.11 88.32 63.18 73.51 60.05];
sma = [70.11	69.16	70.79	91.3	64.81	74.32	61.14];
wma = [68.48	70.24	70.52	91.03	63.86	74.59	61.28];
ema = [75 69.7 70.38 95.52 71.06 73.51 60.33];

thre = [real; sma; wma; ema];

b = bar(X, thre);
ylim([0 100])
ylabel('Accuracy(%)')
set(gcf, 'color', 'white')
legend('Real-time', 'SMA', 'WMA', 'EMA')
%box('off')

%% fre 2-8

all_28 = [ 58.56 60.05 72.69 70.92 56.11 56.39 51.09,53.13];
fix_28 = [ 57.79 63.95 74.09 71.74 54.35 61.23 48.37,50.91];
swi_28 = [ 60.87 48.37 68.48 68.48 61.41 41.85 59.24,59.78];

mean_all_28 = mean(all_28);
mean_fix_28 = mean(fix_28);
mean_swi_28 = mean(swi_28);

std_all_28 = std(all_28)/sqrt(length(all_28));
std_fix_28 = std(fix_28)/sqrt(length(fix_28));
std_swi_28 = std(swi_28)/sqrt(length(swi_28));

%% 3 condition
% real-time
real_all = [67.39 67.53 70.11 88.32 63.18 60.05 73.51];
real_mean_all = 70.01;
% fre 2-8
low_all = [ 58.56 60.05 72.69 70.92 56.11 56.39 51.09];
low_mean_all = 60.83;
% EMA
ema_all = [75 69.7 70.38 95.52 71.06 60.33 73.51];
ema_mean_all = 73.64;

X = categorical({'Real-time','High-freq','EMA'});
X = reordercats(X,{'Real-time','High-freq','EMA'});

thre = [real_mean_all; low_mean_all; ema_mean_all];
y = [real_all; low_all; ema_all];

b=zeros(3,1);

figure
b(1) = bar(X(1), real_mean_all,'FaceColor', 'flat');  hold on
b(2) = bar(X(2), low_mean_all, 'FaceColor', 'flat'); 
b(3) = bar(X(3), ema_mean_all, 'FaceColor', 'flat'); 
plot(X, y, '--ok', 'MarkerFaceColor', [1 1 1] ,'LineWidth', 1);
ylim([0 100]);
ylabel('Accuracy(%)');
set(gcf, 'color', 'white'); hold off
hleg = legend(b,'Real-time', 'High-freq', 'EMA', 'Location', 'southeast');
set(hleg, 'FontSize', 12);


%% all/fix /swi - real/ma

% real data
real_all = [67.39 67.53 70.11 88.32 63.18 60.05 73.51, 76.09];
real_all_m = mean(real_all);
real_all_std = std(real_all);
real_all_se = std(real_all)/sqrt(length(real_all));

real_fix = [ 68.12 69.02 72.83 91.85 63.59 62.50 76.09,77.54];
real_fix_m = mean(real_fix);
real_fix_std = std(real_fix);
real_fix_se = std(real_fix)/sqrt(length(real_fix));

real_swi = [ 65.22 63.04 61.96 77.72 61.96 52.72 65.76,71.74];
real_swi_m = mean(real_swi);
real_swi_std = std(real_swi);
real_swi_se = std(real_swi)/sqrt(length(real_swi));

% ema
ema_all = [75 69.7 70.38 95.52 71.06 60.33 73.51,86.55];
ema_all_m = mean(ema_all);
ema_all_std = std(ema_all);
ema_all_se = std(ema_all)/sqrt(length(ema_all));

ema_fix = [75.54 68.12 75.36 96.2 71.92 63.04 74.64,88.95];
ema_fix_m = mean(ema_fix);
ema_fix_std = std(ema_fix);
ema_fix_se = std(ema_fix)/sqrt(length(ema_fix));

ema_swi = [75.37 74.46 55.43 93.48 68.48 52.17 70.11,79.35];
ema_swi_m = mean(ema_swi);
ema_swi_std = std(ema_swi);
ema_swi_se = std(ema_swi)/sqrt(length(ema_swi));

% high
all_28 = [ 58.56 60.05 72.69 70.92 56.11 56.39 51.09];
fix_28 = [ 57.79 63.95 74.09 71.74 54.35 61.23 48.37];
swi_28 = [ 60.87 48.37 68.48 68.48 61.41 41.85 59.24];

all_28_m = mean(all_28);
fix_28_m = mean(fix_28);
swi_28_m = mean(swi_28);

s_all_28 = std(all_28);
s_fix_28 = std(fix_28);
s_swi_28 = std(swi_28);

%% all/fix /swi - real/ma - plot
X = categorical({'Overall','Fixed','Switching'});
X = reordercats(X,{'Overall','Fixed','Switching'});

% X = categorical({'High freq','Low freq'});
% X = reordercats(X,{'High freq','Low freq'});

re_all = [real_all_m; ema_all_m];
re_fix = [real_fix_m; ema_fix_m];
re_swi = [real_swi_m; ema_swi_m;];
re_se_all = [real_all_se,ema_all_se; real_fix_se,ema_fix_se; real_swi_se,ema_swi_se];

% re_all_std_high = [real_all_m + real_all_std; ema_all_m + ema_all_std];
% re_fix_std_high = [real_fix_m + real_fix_std; ema_fix_m + ema_fix_std]; 
% re_swi_std_high = [real_swi_m + real_swi_std; ema_swi_m + ema_swi_std]; 
% re_all_std_low = [real_all_m - real_all_std; ema_all_m - ema_all_std];
% re_fix_std_low = [real_fix_m - real_fix_std; ema_fix_m - ema_fix_std]; 
% re_swi_std_low = [real_swi_m - real_swi_std; ema_swi_m - ema_swi_std]; 
% std_all_high = [re_all_std_high, re_fix_std_high, re_swi_std_high];
% std_all_low = [re_all_std_low, re_fix_std_low, re_swi_std_low];

barp = [re_all, re_fix, re_swi];

figure
b = bar(barp', 'grouped', 'FaceColor', 'flat'); hold on
b(1).CData(1,:) = [0 0.4470 0.7410]
b(2).CData(1,:) = [0.9290 0.6940 0.1250]
b(2).CData(3,:) = [0.9290 0.6940 0.1250]
b(2).CData(2,:) = [0.9290 0.6940 0.1250]

[ngroups, nbars] = size(barp');

x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% b(1) = bar(X, re_all ,'FaceColor', 'flat');  hold on
% b(2) = bar(X, re_fix ,'FaceColor', 'flat');
% b(3) = bar(X(3), re_swi ,'FaceColor', 'flat');
%plot(X, fre_all, '--ok', 'MarkerFaceColor', [1 1 1] ,'LineWidth', 1);
errorbar(x', barp', re_se_all, 'k','linestyle','none', 'lineWidth', 1);
ylim([0 100]);
ylabel('Accuracy(%)');
set(gcf, 'color', 'white'); hold off
hleg = legend(b,'Real-time', 'EMA', 'Location', 'northeast');
set(hleg, 'FontSize', 12);
set(gca, 'xticklabel', X);







