% Sensitivity to switching auditory attention 

%% Load correlation values
EMAcorr_att = readNPY('C:\Users\LeeJiWon\Desktop\hykist\AAD\OnlineAAD\Recording data\save_data\Allcorr_att_EMA_ceda.npy');
EMAcorr_unatt = readNPY('C:\Users\LeeJiWon\Desktop\hykist\AAD\OnlineAAD\Recording data\save_data\Allcorr_unatt_EMA_ceda.npy');


%% EMA - sensitivity 확인 
% 전화점 앞뒤 5초 accuracy 
Ptrial = [30, 31, 33, 27]; % trial 별 switching point 
tr = [27, 28, 29, 30];

width = 10;

Sensit_att_fr = {};  Sensit_utt_fr = {};
Sensit_att_bc = {};  Sensit_utt_bc = {};
Sacc_tr_fr = [];     Sacc_tr_bc = [];
acc_fr = [];     acc_bc = [];
Compar_bc = {};
for i = 1:length(Ptrial)
    point = Ptrial(i)-14;    
%     point = 1;
    Sensit_att_fr{i} = EMAcorr_att(:,tr(i)-14, 1:point-1);
    Sensit_utt_fr{i} = EMAcorr_unatt(:,tr(i)-14, 1:point-1);
    Sensit_att_bc{i} = EMAcorr_att(:,tr(i)-14, point:end);
    Sensit_utt_bc{i} = EMAcorr_unatt(:,tr(i)-14, point:end);
end
    
for sub = 1:size(EMAcorr_att,1)
    for t = 1:size(Sensit_att_fr,2)
        acc_fr =[]; acc_bc = [];
        for wf = 1:size(Sensit_att_fr{t},3)
            att_fr(wf) = Sensit_att_fr{t}(sub,1,wf);
            utt_fr(wf) = Sensit_utt_fr{t}(sub,1,wf);
                                  
            if att_fr(wf) > utt_fr(wf)
                acc_fr(wf) = 1;
            else acc_fr(wf) = 0;             
            end
            
        end
        for wb = 1:size(Sensit_att_bc{t},3)
            att_bc(wb) = Sensit_att_bc{t}(sub,1,wb);
            utt_bc(wb) = Sensit_utt_bc{t}(sub,1,wb);
            
            if att_bc(wb) > utt_bc(wb)
                acc_bc(wb) = 1;
            else acc_bc(wb) = 0;             
            end
        end          
        Indi_Acc_swi(sub,t) = mean([acc_fr, acc_bc])*100;
        Sacc_tr_fr(sub,t) = mean(acc_fr)*100;
        Sacc_tr_bc(sub,t) = mean(acc_bc)*100;
        Compar_bc{sub,t} = acc_bc;
    end
end

latency = [];
latc = [];
f_one = {};
f = [];
f_ix = [];
for s = 1:size(Compar_bc,1)
    for t = 1:size(Compar_bc,2)
        f_one{s,t} = find(Compar_bc{s,t} == 1);
        try f_one{s,t}(5)-f_one{s,t}(1) == 4 && length(f_one{s,t}) > 1
            latency(s,t) = f_one{s,t}(1);
        catch %length(f_one) < 1
            latency(s,t) = NaN(1);
            f = [f,1];
        end
    end
    % standard
    idx = find(Sacc_tr_fr(s,:) > 55.99); 
%     idx = find(Indi_Acc_swi(s,:) > 55.99);
    f_ix = [f_ix, (4-length(idx))];
    latc(s) = nanmean(latency(s,idx));
    indi_latency{s} = latency(s,idx);
end

disp('-------------------------')
nanmean(latc)
sensitivi_sem = std(latc)/sqrt(length(latc))
mean(mean(Sacc_tr_fr))
mean(mean(Sacc_tr_bc))

