%%
clear

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

%%

figure
plot(tlag, mean_acc*100, '-ok');
ylim([30,100])
grid on
set(gcf, 'color', 'white')

%%
figure
plot(tlag, all_mean*100, '-ok');
ylim([50,100])
grid on
set(gcf, 'color', 'white')
title('all')

%%
figure
plot(tlag, all_mean_f*100, '-ok');
ylim([50,100])
grid on
set(gcf, 'color', 'white')
title('flase')

