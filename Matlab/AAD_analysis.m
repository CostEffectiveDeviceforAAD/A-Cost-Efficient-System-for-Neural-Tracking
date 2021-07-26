%% Real-time AAD analysis
clear

sub = '_01LJW'
subject = 3

%% Behavior
file = strcat('Behavior',sub, '.mat');
load(file)

beh = Behavior(:,3:4);

correct = find(beh=='T');

accb = (length(correct)/60)*100;
  

%% Accuracy Plot
%% Accuracy

file = strcat('Accuracy',sub, '.mat');
load(file)  % Acc

C = readtable('C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\result.xlsx');

overall = mean(Acc*100);
fixed = mean(Acc(1:12)*100);
switching = mean(Acc(13:16)*100);

data = [C; table(subject,overall, fixed, switching)]

chance = 52.99

%% Write to Excel file 
writetable(data, 'C:\Users\LeeJiWon\Desktop\OpenBCI\Recording data\result.xlsx');

%% 1 subject
%C = readtable('result_sample.xlsx');

X = categorical({'Overall','Fixed','Switching'});
Y = table2array(C(subject, 2:end));

b = bar(X, Y);
grid on
ylim([0 80])
ylabel('Accuracy(%)')
% refline([0, chance]);
title('Sub1')

%% multi subject

%C = readtable('result_sample.xlsx');
Y = []
X = []

for sub = 1:length(data.subject)
    name = strcat('Sub', num2str(sub));
    X = [X, categorical({name})];
    Y = [Y; table2array(C(sub,2:end))];
    
end

b = bar(X, Y);
grid on
legend('overall','fixed','switching');
ylim([0 80])
ylabel('Accuracy(%)')
xlabel('Subject')
%%
clear

load 'Accuracy.mat' 
train_trial_ori = mean(Acc(1:14));
test_ori = mean(Acc(15:30));
train_mean_ori = mean(Acc(31:end));

%% channel search
clear
load 'Accuracy_chloc.mat' 

for i = 1:size(Acc,1)
    check_Acc(i) = mean(Acc(i,1:14))*100;
    test_Acc(i)  = mean(Acc(i,15:30))*100;
end

maxval = max(test_Acc);
index = find(test_Acc==maxval);

%% channel graph

%ch = ['T7' 'Fz' 'T8' 'P4' 'P7' 'C3' 'O2'...
%                'Fp1' 'O1' 'Pz' 'F3' 'Cz' 'C4' 'P3' 'P8' 'F4'];        
x = 1:14            
            
%ac = [61.41, 65.76, 66.85, 67.39, 67.39, 66.17, 65.76...
%        63.72, 62.77, 61.96, 60.60, 59.78, 58.42, 56.79, 54.89, 54.62];
     
ac = [71.47, 76.90, 76.22, 75.00, 72.83, 72.28, 70.38, 68.34, 67.39...
        68.75, 69.43, 70.38, 71.88, 69.57]
    
% plot
figure;
plot(x, ac,'-or');
set(gcf, 'color', 'white')
set(gca, 'xtick', [1:14])
%set(gca, 'xticklabel', char('T7', 'Fz', 'T8', 'P4', 'P7', 'C3', 'O2', 'Fp1'...
%                                'O1', 'Pz', 'F3', 'Cz', 'C4', 'P3', 'P8', 'F4'))   %ljw
set(gca, 'xticklabel', char('Fp1', 'P7', 'T7', 'F4', 'Pz/P3', 'F7', 'P4',...
                                'O2', 'Fz', 'C4', 'C3/P8', 'F3', 'T8', 'Cz'))
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

%% 비교
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

X = reordercats(X,{'Overall','Fixed','Switching'});

for sub = 1: length(data.subject)
    Y = [Y; table2array(C(sub,2:end))];
end

Ym = mean(Y, 1)

b = bar(X, Ym);  hold on
plot(X, Y', '-o');
grid on
ylim([0 80])
ylabel('Accuracy(%)')
% refline([0, chance]);
title('Total')

