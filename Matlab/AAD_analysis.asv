%% Real-time AAD analysis

%% Accuracy Plot
%% Accuracy
sub = '_01_LJW'
subject = 5
file = strcat('Accuracy',sub, '.mat');
load(file)  % Acc

C = readtable('result_sample.xlsx');

overall = mean(Acc*100);
fixed = mean(Acc(1:12)*100);
switching = mean(Acc(13:16)*100);

data = [C; table(subject,overall, fixed, switching)]
%%
writetable(data, 'result_sample.xlsx');

%% 1 subject
C = readtable('result_sample.xlsx');

X = categorical({'Overall','Fixed','Switching'});
Y = table2array(C(subject, 2:end));

b = bar(X, Y);
grid on
ylim([0 100])
ylabel('Accuracy(%)')
title('Sub1')

%% multi subject

C = readtable('result_sample.xlsx');
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
ylim([0 100])
ylabel('Accuracy(%)')
xlabel('Subject')


