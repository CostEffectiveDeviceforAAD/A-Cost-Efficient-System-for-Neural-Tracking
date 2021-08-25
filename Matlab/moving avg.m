
clear

load('C:\Users\user\Desktop\hy-kist\OpenBCI\AAD\Matlab\all_acc.mat')
%% moving avg - overall

% moving average
mov_l_1 = movmean(corr_l, 1);
mov_r_1 = movmean(corr_r, 1);

mov_l_3 = movmean(corr_l, 3);
mov_r_3 = movmean(corr_r, 3);

mov_l_5 = movmean(corr_l, 5);
mov_r_5 = movmean(corr_r, 5);

mov_l_7 = movmean(corr_l, 7);
mov_r_7 = movmean(corr_r, 7);

%%
% extimate acc
for tr = 1:16
    for i = 1:46
        if mov_l_1(tr,i) > mov_r_1(tr,i)
            acc_1(tr,i) = 1;
        else acc_1(tr,i) = 0;
        end
    end
end
for tr = 1:16
    for i = 1:46
        if mov_l_3(tr,i) > mov_r_3(tr,i)
            acc_3(tr,i) = 1;
        else acc_3(tr,i) = 0;
        end
    end
end
for tr = 1:16
    for i = 1:46
        if mov_l_5(tr,i) > mov_r_5(tr,i)
            acc_5(tr,i) = 1;
        else acc_5(tr,i) = 0;
        end
    end
end
for tr = 1:16
    for i = 1:46
        if mov_l_7(tr,i) > mov_r_7(tr,i)
            acc_7(tr,i) = 1;
        else acc_7(tr,i) = 0;
        end
    end
end

%%
mean_acc_1 = mean(acc_1,2);
mean_acc_1_all = mean(mean_acc_1);
mean_acc_1_fix = mean(mean_acc_1(1:12));
mean_acc_1_swi = mean(mean_acc_1(13:end));

mean_acc_3 = mean(acc_3,2);
mean_acc_3_all = mean(mean_acc_3);
mean_acc_3_fix = mean(mean_acc_3(1:12));
mean_acc_3_swi = mean(mean_acc_3(13:end));

mean_acc_5 = mean(acc_5,2);
mean_acc_5_all = mean(mean_acc_5);
mean_acc_5_fix = mean(mean_acc_5(1:12));
mean_acc_5_swi = mean(mean_acc_5(13:end));

mean_acc_7 = mean(acc_7,2);
mean_acc_7_all = mean(mean_acc_7);
mean_acc_7_fix = mean(mean_acc_7(1:12));
mean_acc_7_swi = mean(mean_acc_7(13:end));

%%

% sum
com_1_all = [com_1_all; mean_acc_1_all];
com_1_fix = [com_1_fix; mean_acc_1_fix];
com_1_swi = [com_1_swi; mean_acc_1_swi];

com_3_all = [com_3_all; mean_acc_3_all];
com_3_fix = [com_3_fix; mean_acc_3_fix];
com_3_swi = [com_3_swi; mean_acc_3_swi];

com_5_all = [com_5_all; mean_acc_5_all];
com_5_fix = [com_5_fix; mean_acc_5_fix];
com_5_swi = [com_5_swi; mean_acc_5_swi];

com_7_all = [com_7_all; mean_acc_7_all];
com_7_fix = [com_7_fix; mean_acc_7_fix];
com_7_swi = [com_7_swi; mean_acc_7_swi];

%%
% sub by var

% mean each
all = mean(com,1);

s = std(com*100);

%%
for tr = 1:16
    
    for i = 1:46

        if mov_l_7(tr,i) > mov_r_7(tr,i)
            acc(tr,i) = 1;
        else acc(tr,i) = 0;
        end

    end
end

mov_ac = mean(acc, 2);

%% moving - online
mov_corr_att = []
mov_corr_un = []
mov_acc_all = []
mov_acc_trial = []

for tr = 1:16
    
    % moving average filter
    for i = 1:2:7
        
        for w = 1:46

            if w <= i 
                mov_cor_att(w) = mean(corr_l(tr, 1:w));
                mov_cor_un(w) = mean(corr_r(tr, 1:w));
            else     
                mov_cor_att(w) = mean(corr_l(tr, w-i+1:w));
                mov_cor_un(w) = mean(corr_r(tr, w-i+1:w));
            end
            
            % estimate
            if mov_cor_att(w) > mov_cor_un(w)
                mov_acc(w) = 1;
            else mov_acc(w) = 0;
            
            end
        end     % one sample
        
        mov_corr_att = [mov_corr_att; mov_cor_att];     % 4 by 46
        mov_corr_un = [mov_corr_un; mov_cor_un];
        
        mov_acc_trial = [mov_acc_trial; mean(mov_acc)];     % sample number by one
    end     % end one tiral for 4sample
        mov_acc_all = [mov_acc_all, mov_acc_trial];     % samplw by trials
        mov_acc_trial = [];
end
    
    
mov_acc_overall_1 = mean(mov_acc_all(1,:));
mov_acc_overall_3 = mean(mov_acc_all(2,:));
mov_acc_overall_5 = mean(mov_acc_all(3,:));
mov_acc_overall_7 = mean(mov_acc_all(4,:));
    
mov_acc_fix_1 = mean(mov_acc_all(1,1:12));
mov_acc_fix_3 = mean(mov_acc_all(2,1:12));
mov_acc_fix_5 = mean(mov_acc_all(3,1:12));
mov_acc_fix_7 = mean(mov_acc_all(4,1:12));   
    
mov_acc_swi_1 = mean(mov_acc_all(1,13:end));
mov_acc_swi_3 = mean(mov_acc_all(2,13:end));
mov_acc_swi_5 = mean(mov_acc_all(3,13:end));
mov_acc_swi_7 = mean(mov_acc_all(4,13:end));   

%%

total_acc_all_1 = [67.39,67.53,70.11,88.32,63.18,73.51,60.05];
total_acc_all_3 = [67.39,69.02,69.96,90.49,63.04,74.05,60.73];
total_acc_all_5 = [68.75,70.11,69.70,91.03,63.45,73.64,60.46];
total_acc_all_7 = [70.11,69.16,70.79,91.30,64.81,74.32,61.14];

total_acc_fix_1 = [68.12,69.20,72.83,91.85,63.59,76.09,62.50];
total_acc_fix_3 = [67.93,71.20,74.64,93.30,63.41,76.27,63.77];
total_acc_fix_5 = [68.66,72.49,74.28,93.66,63.95,76.06,63.41];
total_acc_fix_7 = [69.57,71.20,75.54,93.66,65.40,76.27,63.77];

total_acc_swi_1 = [65.22,63.04,61.96,77.72,61.96,65.76,52.72];
total_acc_swi_3 = [65.76,62.50,55.98,82.07,61.96,67.39,51.63];
total_acc_swi_5 = [69.02,63.04,55.98,83.15,61.96,66.30,51.63];
total_acc_swi_7 = [71.74,63.04,56.52,84.24,63.04,68.48,53.26];

total_mean_all_1 = mean(total_acc_all_1);
std_all_1 = std(total_acc_all_1);
total_mean_all_3 = mean(total_acc_all_3);
std_all_3 = std(total_acc_all_3);
total_mean_all_5 = mean(total_acc_all_5);
std_all_5 = std(total_acc_all_5);
total_mean_all_7 = mean(total_acc_all_7);
std_all_7 = std(total_acc_all_7);
    
total_mean_fix_1 = mean(total_acc_fix_1);
std_fix_1 = std(total_acc_fix_1);
total_mean_fix_3 = mean(total_acc_fix_3);
std_fix_3 = std(total_acc_fix_3);
total_mean_fix_5 = mean(total_acc_fix_5);
std_fix_5 = std(total_acc_fix_5);
total_mean_fix_7 = mean(total_acc_fix_7);
std_fix_7 = std(total_acc_fix_7);

total_mean_swi_1 = mean(total_acc_swi_1);
std_swi_1 = std(total_acc_swi_1);
total_mean_swi_3 = mean(total_acc_swi_3);
std_swi_3 = std(total_acc_swi_3);
total_mean_swi_5 = mean(total_acc_swi_5);
std_swi_5 = std(total_acc_swi_5);
total_mean_swi_7 = mean(total_acc_swi_7);
std_swi_7 = std(total_acc_swi_7);

%% weighted moving average

n=7

for tr = 1:16
    for i = 1:46
         if i >= n
             for k = 0:n-1
                % att
                wv_at(k+1) = Corr_att(tr,i-(k))*(n-k);        
                % un
                wv_un(k+1) = Corr_un(tr,i-(k))*(n-k);        
             end
             
             wma_at(tr,i) = sum(wv_at)/(n*(n+1)/2);
             wma_un(tr,i) = sum(wv_un)/(n*(n+1)/2);
             
         else 
             for k = 0:i-1
                 % att
                wv_at(k+1) = Corr_att(tr,i-(k))*(i-k);        
                % un
                wv_un(k+1) = Corr_un(tr,i-(k))*(i-k);        
             end
             
             wma_at(tr,i) = sum(wv_at)/(i*(i+1)/2);
             wma_un(tr,i) = sum(wv_un)/(i*(i+1)/2);
             
         end
         wv_at = [];
         wv_un = [];
        % estimate
        if wma_at(tr,i) > wma_un(tr,i)
            acc(tr,i) = 1;
        else acc(tr,i) = 0;
        end
    end
end

mean_acc = mean(acc,2);

mean_all = mean(mean_acc);
mean_fix = mean(mean_acc(1:12));
mean_swi = mean(mean_acc(13:16));

%% weight

a_3 = std([66.98	68.75	69.7	90.22	63.86	74.05	60.46]);
a_5 = std([67.53	69.7	69.84	90.9	63.86	74.05	60.46]);
a_7 = std([68.48	70.24	70.52	91.03	63.86	74.59	61.28]);

f_3 = std([67.39	70.29	74.46	93.3	64.49	76.09	63.41]);
f_5 = std([67.93	71.92	74.46	93.66	64.31	76.27	63.77]);
f_7 = std([68.3	72.64	75.18	93.66	64.67	76.45	64.49]);

s_3 = std([65.76	64.13	55.43	80.98	61.96	67.93	51.63]);
s_5 = std([66.3	63.04	55.98	82.61	62.5	67.39	50.54]);
s_7 = std([69.02	63.04	56.52	83.15	61.41	69.02	51.63]);

%% offline

a_3 = std([76.36	79.35	82.07	99.73	70.52	92.8	66.71]);
a_5 = std([82.22	85.33	86.68	100	74.18	97.15	72.83]);
a_7 = std([86.14	89.67	88.32	100	74.46	99.18	75]);

f_3 = std([76.81	79.53	85.14	100	69.2	95.29	69.02]);
f_5 = std([84.06	84.6	91.3	100	69.75	99.28	76.45]);
f_7 = std([89.31	88.04	94.02	100	71.74	99.82	78.8]);

s_3 = std([75	78.8	72.83	98.91	74.46	85.33	59.78]);
s_5 = std([76.63	87.5	72.83	100	87.5	90.76	61.69]);
s_7 = std([76.63	94.57	71.2	100	82.61	97.28	63.59]);

%% Exponential MA

ema_at = []
ema_un = []

for tr = 1:16
    
    ema_at(tr,1) = corr_l(tr,1);
    ema_un(tr,1) = corr_r(tr,1);

    for i = 2:46

        a = 2/(i+1);
        ema_at(tr,i) = (a*corr_l(tr,i)) + ((1-a)*ema_at(tr,i-1));
        ema_un(tr,i) = (a*corr_r(tr,i)) + ((1-a)*ema_un(tr,i-1));
        
        if ema_at(tr,i) > ema_un(tr,i)
            ema_corr(tr,i) = 1;
        else ema_corr(tr,i) = 0;end
    end
    
end

ema_all = mean(mean(ema_corr,2))*100;
ema_fix = mean(mean(ema_corr(1:12,:),2))*100;
ema_swi = mean(mean(ema_corr(13:16,:),2))*100;
            
%%

all = [75 69.7 70.38 95.52 71.06 60.33 73.51];
fix = [75.54 68.12 75.36 96.2 71.92 63.04 74.64];
swi = [75.37 74.46 55.43 93.48 68.48 52.17 70.11];

mean_all = mean(all);
std_all = std(all);
mean_fix = mean(fix);
std_fix = std(fix);
mean_swi = mean(swi);
std_swi = std(swi);





