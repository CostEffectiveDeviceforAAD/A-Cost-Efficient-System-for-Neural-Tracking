%% Real-time AAD analysis

sub = '_01LJW'
subject = 1

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

C = readtable('result_sample.xlsx');

overall = mean(Acc*100);
fixed = mean(Acc(1:12)*100);
switching = mean(Acc(13:16)*100);

data = [C; table(subject,overall, fixed, switching)]

chance = 52.99

%%
writetable(data, 'result_sample.xlsx');

%% 1 subject
C = readtable('result_sample.xlsx');

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

%% Envelope
load 'Predict_L.mat'
load 'Allspeech.mat'
i = 20;
L = 960;
x = linspace(0, L/64, L);
%%
for i = 1:5
    figure()
    subplot(311)
    plot(x, Allspeech(31+i,1:960)')
    title('Speech-att')
    subplot(312)
    plot(x, Allspeech(1+i,1:960)')
    title('Speech-unatt')
    subplot(313)
    plot(x, Pre_L{1+i}(1,:)')
    title('Predicted')
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
