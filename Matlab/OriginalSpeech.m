%%
% speech 원본 가공;;



%% Load
s_T = {};
s_J = {};
% for j = 1:30
%     [s_T{j},fs] = audioread('해저'+string(j)+'.wav'); % fs = 44100
%     if j <10
%         [s_J{j},fs] = audioread('0'+string(j)+'_Journey.wav'); 
%     else
%         s_J{j} = audioread(string(j)+'_Journey.wav'); end
%         %fs = 48000
% end


for j = 1:30
%     [s_T{j},fs] = audioread('해저'+string(j)+'.wav'); % fs = 44100
    if j <10
        [s_J{j},fs] = audioread('음원31-'+string(j-1)+'.mp3'); 
    else
        s_J{j} = audioread('음원31-'+string(j-1)+'.mp3'); end
        %fs = 48000
end
%% Adjust length of sound
ss_J = {};
for j = 1:30
    H_len = fs*60;
    
    if length(s_J{j}) > H_len
        ss_J{j} = s_J{j}(1:H_len);
    else
        ss_J{j} = cat(1,s_J{j}(:,1), zeros(H_len-length(s_J{j}),1));
    end
end

% Adjust RMS of all speech to H_rms 
H_rms = 0.02;
for j = 1:length(ss_J)
    ss_J{j} = (ss_J{j}./rms(ss_J{j}).*H_rms); end

% save each speech
for j = 1:length(ss_J)

    %audiowrite(string(j)+'_Journey.wav', ss_J{j},fs);
    audiowrite(string(j)+'_Twenty.wav', ss_J{j},fs);
end

%% sat sound
% load
for j = 1:30
    [s2(j,:),fs] = audioread('C:\Users\LeeJiWon\Desktop\hykist\AAD\Speech\2022\2022_수능31-60\반향제거\진짜\'+string(j+30)+'_SAT.wav');
    [s1(j,:),fs] = audioread('C:\Users\LeeJiWon\Desktop\hykist\AAD\Speech\2022\2022_수능1-30\'+string(j)+'_SAT.wav');
end

for j = 1:30
    [s1(j,:),fs] = audioread('C:\Users\LeeJiWon\Desktop\hykist\AAD\Speech\2022\2022_Journey\'+string(j)+'_Journey.wav');
    [s2(j,:),fs] = audioread('C:\Users\LeeJiWon\Desktop\hykist\AAD\Speech\2022\2022_Twenty\'+string(j)+'_Twenty.wav');
end

% merge
sat = {};
for j = 1:30
    sat{j} = [s1(j,:); s2(j,:)];
    if j <10
        audiowrite('00'+string(j)+'_SAT.wav',sat{j}',fs);
    else
        audiowrite('0'+string(j)+'_SAT.wav',sat{j}',fs);
    end
end
save('SAT_speech.mat','sat')

aak = {};
for j = 1:30
    aak{j} = [s2(j,:); s1(j,:)];
    
    if j <10
        audiowrite('00'+string(j)+'_AAK.wav', aak{j}',fs);
    else
        audiowrite('0'+string(j)+'_AAK.wav', aak{j}',fs);
    end
end


%%

% sampling rate 44100 으로 맞춤
for j = 1:30
    s_J{j} = resample(s_J{j},44100,48000);end
fs = 44100;

% Adjust RMS of all speech to H_rms 
H_rms = 0.02;
for j = 1:length(s_T)
    Rs_T{j} = (s_T{j}./rms(s_T{j}).*H_rms); end
Rs_T{length(s_T)+1} = fs;

% Terminate silence 0.5s

diff = diff(abs(Rspeech{i}));
diffIdx = find(abs(diff) < 0.001);

%
for j = 1:length(Rs_T)-1
    a = abs(Rs_T{j});
    aa = a;
    t = length(a);
    m = mean(a);
    while t > 1+(fs*0.5)

        if a(t) < m   % 낮은 지점 찾기
            l = find(a(t-fs*0.5:t) < m);    % 그 부분부터 0.5초 길이동안 다 낮으면 silence로 인식

            if length(l) >= fs*0.5
                disp(['!!']) 
                Rs_T{j}(t-fs*0.4:t-fs*0.2) = []; % 0.2초 길이 묵음 제거
                t = t-fs*0.5;
            end
        end
        t = t-1;
    end
end
save('adj_speech0.4_0.2.mat','Rs_T');



%%
ss_J={};
lek= {};
H_len = fs*60;

lek{1} = H_len - length(s_J{1});

i=i+1;

lek{i} = H_len - length(s_J{i}(lek{i-1}-fs*3:end));
ss_J{i} = cat(1,s_J{i}(lek{i-1}-fs*3:end),s_J{i+1}(1:lek{i}));

% lek 넘친후
ss_J{i} = cat(1,s_J{i}(lek{i-1}-fs*3:end),s_J{i+1});
lek2 = H_len - length(ss_J{i});
ss_J{i} = cat(1,ss_J{i}, s_J{i+2}(1:lek2));

% 한번하고 그 다음
ss_J{i} = s_J{i+1}(lek2-fs*3:end);
ss_J{i} = cat(1,ss_J{i},s_J{i+2}(1:H_len-length(ss_J{i})));


sound(ss_J{i},fs);







