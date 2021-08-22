%%
clear

all_acc=[]
all_false=[]
%%
acc = Acc(:,15:30);
mean_acc = mean(acc,2);

tlag = flip(-[-406.25 , -390.625, -375, -359.375, -343.75 , -328.125, ...
       -312.5  , -296.875, -281.25 , -265.625, -250, -234.375, ...
       -218.75 , -203.125, -187.5  , -171.875, -156.25 , -140.625, ...
       -125, -109.375,  -93.75 ,  -78.125,  -62.5  ,  -46.875, ...
        -31.25 ,  -15.625,    0])';  
    
%%
    
all_acc = [all_acc, mean_acc];

%%
all_mean = mean(all_acc,2);

%%
all_false = [all_false, mean_acc];

%%
all_mean_f = mean(all_false,2);

%% one

figure
plot(tlag, mean_acc*100, '-ok');  %'LineWidth', 2);
ylim([30,80])
ylabel('Accuracy (%)')
xlabel('Time-lags (ms)')
title('Accuracy at individual time-lags')
set(gcf, 'color', 'white')

%%
figure
plot(tlag, all_mean.*100, 'k', 'LineWidth', 3);
ylim([30,80])
set(gcf, 'color', 'white')
xlim([0 407])
ylabel('Accuracy (%)')
xlabel('Time-lags (ms)')
title('Accuracy at individual time-lags')
box('off')

%%
figure
plot(tlag, all_mean_f*100,'LineWidth', 3);
ylim([30,80])
set(gcf, 'color', 'white')
ylabel('Accuracy (%)')
xlabel('Time-lags (ms)')
colororder([0.7 0.7 0.7])
title('Accuracy at individual time-lags (False)')

%% two
figure
plot(tlag, all_mean*100, 'LineWidth', 3); hold on
plot(tlag, all_mean_un*100, 'LineWidth', 3);
ylim([30,80])
set(gcf, 'color', 'white')
xlim([0 407])
ylabel('Accuracy (%)')
xlabel('Time-lags (ms)')
colororder([0 0 0; 0.7 0.7 0.7])
legend('Attended', 'Unattended')
%title('Accuracy at individual time-lags')
box('off')

%%

for i = 1:7
    subplot(4,2,i)
    plot(tlag, all_acc(:,i)*100, '-ok')
    xlim([0 407])
    ylim([30,90])
    title(num2str(i))
end

