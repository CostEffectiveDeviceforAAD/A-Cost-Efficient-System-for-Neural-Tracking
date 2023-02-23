
%%

[sound,fs] = audioread('001_L1.wav');

sound = resample(sound, 125,fs);

%% Sound figure

plot(sound(125*31:125*32,2),'r')
axis off
set(gca, 'box', 'off','color', 'none', 'Xtick', [], 'Ytick', [])
set(gcf,'InvertHardcopy','off')
set(gcf, 'color', 'none')

% saveas(d, 'sound.eps')

%% accuracy figure_4

%-- Real-Time --%
% all
for i = 1:length(DATA)
    RealT_all(i) = DATA(i).RealTime_all;
end
RealT_all_mean = mean(RealT_all);
RealT_all_std = std(RealT_all);

% fix
for i = 1:length(DATA)
    RealT_fix(i) = DATA(i).RealTime_fix;
end
RealT_fix_mean = mean(RealT_fix);
RealT_fix_std = std(RealT_fix);

% switch
for i = 1:length(DATA)
    RealT_swi(i) = DATA(i).RealTime_switch;
end
RealT_swi_mean = mean(RealT_swi);
RealT_swi_std = std(RealT_swi);

%-- EMA --%
% all
for i = 1:length(DATA)
    EMA_all(i) = DATA(i).EMA_all;
end
EMA_all_mean = mean(EMA_all);
EMA_all_std = std(EMA_all);

% fix
for i = 1:length(DATA)
    EMA_fix(i) = DATA(i).EMA_fix;
end
EMA_fix_mean = mean(EMA_fix);
EMA_fix_std = std(EMA_fix);

% switch
for i = 1:length(DATA)
    EMA_swi(i) = DATA(i).EMA_switch;
end
EMA_swi_mean = mean(EMA_swi);
EMA_swi_std = std(EMA_swi);

%-- Unattended --%
RealT_unatt = [44.88,42.38,50.43,48.84,56.81,46.81,64.78,45.07,46.52];
RealT_unatt_mean = mean(RealT_unatt);
RealT_unatt_std = std(RealT_unatt);


%% save
for i = 1:length(RealT_all)
    DATA(i).RealTime_all = RealT_all(i);
    DATA(i).RealTime_fix = RealT_fix(i);
    DATA(i).RealTime_switch = RealT_swi(i);
    DATA(i).EMA_all = EMA_all(i);
    DATA(i).EMA_fix = EMA_fix(i);
    DATA(i).EMA_switch = EMA_swi(i);
    DATA(i).Unattended = RealT_unatt(i);
end
   
%%
% xlable
X = categorical({'sub02','sub05','sub06','sub07','sub09','sub10','sub11','sub12','sub13'});
X = reordercats(X,{'sub02','sub05','sub06','sub07','sub09','sub10','sub11''sub12','sub13'});

% load
for i = 1:9
    RealT_all(i) = DATA(i).RealTime_all;
    EMA_all(i) = DATA(i).EMA_all;
end

RealT_all = sort(RealT_all,'descend');
EMA_all = sort(EMA_all, 'descend');

% plot
clf
figure(1)
bar(X, [RealT_all; EMA_all])
set(gca,'color','none', 'box', 'off', 'Xtick', [], 'Ytick', [40 60 80 100]);
set(gcf, 'color','none')
legend('Real-time', 'EMA')
yline(52, '--', 'lineWidth', 1);
ylim([30 100])

% plot_mean
acc = [RealT_all_mean, RealT_fix_mean, RealT_swi_mean; EMA_all_mean, EMA_fix_mean, EMA_swi_mean]';
std = [RealT_all_std, RealT_fix_std, RealT_swi_std; EMA_all_std, EMA_fix_std, EMA_swi_std]';

%
figure
b = bar(acc, 'FaceColor', 'flat'); hold on

[ngroups, nbars] = size(acc);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

errorbar(x', acc, std, 'k','linestyle','none', 'lineWidth', 1);
% b.CData(1,:) = [0 0.4470 0.7410]
% b.CData(2,:) = [0.8500 0.3250 0.0980]
set(gca,'color','none', 'box', 'off', 'Xtick', [], 'Ytick',[40 60 80 100]);
set(gcf, 'color','none')
ylim([30 100])


%% figure.4 _ ver.2

%% only real-time - 3 types

figure(1)
clf

% data set
% data = [RealT_all_mean; RealT_fix_mean; RealT_swi_mean];

data = [RealT_all_mean; EMA_all_mean];

% bar plot
b = bar( data ); hold on
b.FaceColor = 'flat';
b.CData(1,:) = [0.8500 0.3250 0.0980];
b.CData(2,:) = [0.9290 0.6940 0.1250];
set(gcf, 'color', 'none')
set(gca,'color','none', 'Ytick',[50 60 70 80  90]);
ylim([50 90])
yline(52, '--', 'lineWidth', 1);

% errorbar
% plot(1*ones(1,2), RealT(1)+[-RealT_all_std,RealT_all_std], 'color', 'k', 'LineWidth', 5)
% plot(2*ones(1,2), RealT(2)+[-RealT_fix_std,RealT_fix_std], 'color', 'k', 'LineWidth', 5)
% plot(3*ones(1,2), RealT(3)+[-RealT_swi_std,RealT_swi_std], 'color', 'k', 'LineWidth', 5)

plot(1*ones(1,2), data(1)+[-RealT_all_std, RealT_all_std], 'color', 'k', 'LineWidth', 5)
plot(2*ones(1,2), data(2)+[-EMA_all_std, EMA_all_std], 'color', 'k', 'LineWidth', 5)



%% individual subject at 3 type

% sub_type = [RealT_all; RealT_fix; RealT_swi]';
% sub_type = sort(sub_type,'descend');

% bar plot
b = bar( sub_type, 'FaceColor', 'flat' ); hold on

% b.CData(2,:) = [0.8500 0.3250 0.0980];
% b.CData(3,:) = [0.9290 0.6940 0.1250];
set(gcf, 'color', 'none')
set(gca,'color','none', 'Ytick',[50 60 70 80 90 100]);
ylim([50 100])
yline(52, '--', 'lineWidth', 1);

%% EMA / Real-time
figure(2)
clf

plot( sort(EMA_all, 'descend'), '-o', 'color', 'k', 'LineWidth', 3); hold on
plot( sort(RealT_all, 'descend'), '-o', 'color', 'b', 'LineWidth', 3)
set(gcf, 'color', 'none')
set(gca, 'color', 'none', 'Ytick',[50 60 70 80 90 100]);
xlim([0 10])
ylim([50 100])

%% EMA 결과로 condition 별, individual
%% trial type 별
% load('C:\Users\LeeJiWon\Desktop\OpenBCI\AAD\Matlab\DATA(OpenBCI)_sub9.mat')

all_ema = [];
fix_ema = [];
swi_ema = [];
for sub = 1:length(DATA)
    all_ema = [all_ema, DATA(sub).EMA_all];
    fix_ema = [fix_ema, DATA(sub).EMA_fix];
    swi_ema = [swi_ema, DATA(sub).EMA_switch];
end
all_std = std(all_ema)/sqrt(length(all_ema));
fix_std = std(fix_ema)/sqrt(length(fix_ema));
swi_std = std(swi_ema)/sqrt(length(swi_ema));
acc = [all_ema; fix_ema; swi_ema];
std = [all_std; fix_std; swi_std];

X = categorical({'All','Fixed', 'Switching'});
X = reordercats(X,{'All','Fixed', 'Switching'});

figure(40)
clf
b = bar(X, mean(acc,2), 'Facecolor', 'flat'); hold on
b.CData(1,:) = [0 0.4470 0.7410]
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.9290 0.6940 0.1250];
ylim([40 100])

% error bar
er = errorbar(X, mean(acc,2), -std, std);  
er.Color = [0 0 0];  
er.LineStyle = 'none';  
er.LineWidth = 1;  
yline(52.99, ':')
ylabel('Decoder accuracy (%)')
xlabel('Trial type')

%% individual
[acc_d,q] = sort(acc(1,:), 'descend');

for i = 1:length(q)
    acc_d(2,i) = acc(2,q(i));
    acc_d(3,i) = acc(3,q(i));
end

figure(41)
clf
bar(acc_d')
yline(52.99, ':')
ylabel('Decoder accuracy (%)')
ylim([40 100])
legend(X)




