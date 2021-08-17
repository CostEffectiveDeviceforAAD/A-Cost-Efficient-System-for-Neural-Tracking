
clear

load('C:\Users\user\Desktop\hy-kist\OpenBCI\AAD\Matlab\all_acc.mat')
%%
% all acc
mean_acc = mean(accuracy,2);
% moving average
mov_1 = movmean(accuracy, 1);
mov_3 = movmean(accuracy, 3);
mov_5 = movmean(accuracy, 5);
mov_7 = movmean(accuracy, 7);

mean_mov_1 = mean(mov_1,2);
mean_mov_3 = mean(mov_3,2);
mean_mov_5 = mean(mov_5,2);
mean_mov_7 = mean(mov_7,2);

% sum
com = [mean_acc, mean_mov_1, mean_mov_3, mean_mov_5, mean_mov_7];
% sub by var

% mean each
all = mean(com,1);
