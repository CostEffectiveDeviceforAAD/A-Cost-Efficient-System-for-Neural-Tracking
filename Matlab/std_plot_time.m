%%

st_at = std(all_acc')'*100;     % one by individual time-lags
tmp_at = mean(all_acc,2)*100;

%%

st_un = std(all_acc')'*100;     % one by individual time-lags
tmp_un = mean(all_acc,2)*100;

%%
tlag = flip(-[-406.25 , -390.625, -375, -359.375, -343.75 , -328.125, ...
       -312.5  , -296.875, -281.25 , -265.625, -250, -234.375, ...
       -218.75 , -203.125, -187.5  , -171.875, -156.25 , -140.625, ...
       -125, -109.375,  -93.75 ,  -78.125,  -62.5  ,  -46.875, ...
        -31.25 ,  -15.625,    0])';  

    
x = [tlag(1:end,1); tlag(end:-1:1,1)];
y_at = [tmp_at(1:end)+st_at(1:end); tmp_at(end:-1:1)-st_at(end:-1:1)];
y_un = [tmp_un(1:end)+st_un(1:end); tmp_un(end:-1:1)-st_un(end:-1:1)];
%%

figure
p = fill(x' ,y_at', 'red');
p.FaceColor = [0.6350 0.0780 0.1840];      
p.EdgeColor = 'none'; 
p.FaceAlpha = 0.1

hold on
p = fill(x' ,y_un', 'blue');
p.FaceColor = [0 0.4470 0.7410];
p.EdgeColor = 'none'; 
p.FaceAlpha = 0.1

plot(tlag, tmp_at, '-r', 'LineWidth', 2.5);
plot(tlag, tmp_un, '-b', 'LineWidth', 2.5);

legend('','','Attended', 'Unattended')
ylim([30,90])
set(gcf, 'color', 'white')
xlim([0 407])
ylabel('Accuracy (%)')
xlabel('Time-lags (ms)')
title('Accuracy at individual time-lags')
box('off')



