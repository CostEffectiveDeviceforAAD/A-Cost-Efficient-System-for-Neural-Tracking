%% Figure 2.  Results of Frequency Search
%
% working directory should be set "~/Supplementary"

pwd = 'C:\Users\LeeJiWon\Desktop\hykist\AAD\BSC\Supplementary-20221123T035500Z-001\Supplementary'
%% Data Preprocessing for Plotting

% ================================ AAK ====================================
% subject number
subNum = [1 3:11];

% init mtx
AAK_Att = zeros(length(subNum), 19, 19);

for subi = 1:length(subNum)

    % load result
    load([pwd '/AAK/AAK' num2str(subNum(subi)) '_ResultsOverFreq.mat'])

    % calculate decoding accuracy
    for i = 1:size(ResultsOverFreq, 1)
        for j = 1:size(ResultsOverFreq, 2)
            if ~(j < i)
                AAK_Att(subi,i,j)=mean(ResultsOverFreq(i,j).result.acc);
            end
        end
    end
end   


% =============================== AADC ====================================
% subject number
subNum = 2:11;

% init mtx
AADC_Att = zeros(length(subNum), 19, 19);

for subi = 1:length(subNum)
    
    % load result
    load([pwd '/AADC/AADC' num2str(subNum(subi)) '_ResultsOverFreq.mat'])

    % calculate decoding accuracy
    for i = 1:size(ResultsOverFreq, 1)
        for j = 1:size(ResultsOverFreq, 2)
            if ~(j < i)
                AADC_Att(subi,i,j)=mean(ResultsOverFreq(i,j).result.acc);
            end
        end
    end
end    


% ================================ EXP ====================================
% subject number
subNum = 1:6;

% init mtx
Exp_Att = zeros(length(subNum), 19, 19);

for subi = 1:length(subNum)
    
    % load result
    load([pwd '/Exp/Exp' num2str(subNum(subi)) '_ResultsOverFreq.mat'])

    % calculate decoding accuracy
    for i = 1:size(ResultsOverFreq, 1)
        for j = 1:size(ResultsOverFreq, 2)
            if ~(j < i)
                Exp_Att(subi,i,j)=mean(ResultsOverFreq(i,j).result.acc);
            end
        end
    end
end   

% Concat Data
Att = cat(1, AAK_Att, AADC_Att);
Att = cat(1, Att, Exp_Att);

%% A. Frequency Search Heatmap for Online AAD Model

% get screen size
screensize = get( 0, 'Screensize' );

% axes
xs = 0.5:0.5:9.5; ys = 1:0.5:10;

% figure
figure(1)
clf
% set(gcf, 'color', 'w')
set(gcf, 'color', 'none', 'Position', [1 1 screensize(3)/2 screensize(4)])

% heatmap
dat=squeeze(mean(Att,1))'*100;
imagesc(xs(1:15), ys(1:15), dat(1:15,1:15))
axis xy; axis square
c_h=colorbar;
c_h.Label.String = ' Decoding Accuracy (%)';
colormap jet
xticks(xs); yticks(ys);
set(gca, 'clim', [50 100], 'FontSize', 20, 'color', 'none')
xlabel(' Low Frequency Cut-off (Hz) ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' High Frequency Cut-off (Hz) ', 'FontSize', 25, 'FontWeight', 'bold')

% print('heatmap', '-depsc')

%% B. Bar Plots for Comparisons

% figure
figure(2)
clf
% set(gcf, 'color', 'w')
set(gcf, 'color', 'none', 'Position', [1 1 screensize(3)/3 screensize(4)*2/3])


% Previous
LowPrev  = find(xs==2);
HighPrev = find(ys==8);

% New Band
LowIdx_Att = find(xs==0.5); 
HighIdx_Att = find(ys==8); 


% data settings
Acc = [squeeze(Att(:, LowPrev, HighPrev)) squeeze(Att(:, LowIdx_Att, HighIdx_Att))]*100;
ave    = mean(Acc);
% stdev  = std(Acc, [], 1); % sd
stdev  = std(Acc, [], 1)/sqrt(size(Acc,1)); % sem
lower  = -stdev; upper = stdev;

% barplot
b=bar( ave ); hold on;
b.FaceColor = 'flat';
b.CData(2,:) = [0.8500 0.3250 0.0980];

% errorbar
plot(1*ones(1,2), ave(1)+[lower(1),upper(1)], 'color', 'k', 'LineWidth', 5)
plot(2*ones(1,2), ave(2)+[lower(2),upper(2)], 'color', 'k', 'LineWidth', 5)

% % Individual Plot
% jitter = randn(size(Acc,1),2)/10;
% for i = 1:size(Acc,1)
%     plot([1,2]+jitter(i,:), Acc(i,:), 'o--', 'color', [0.7 0.7 0.7], 'LineWidth', 2, 'MarkerSize', 10); 
% end

% Elaboration
xticklabels({'2-8 Hz', '0.5-7.5 Hz'})
set(gca, 'ylim', [50 100], 'FontSize', 18, 'color', 'none')
xlabel(' Frequency Band ', 'FontSize', 25, 'FontWeight', 'bold')
ylabel(' Decoding Performance ', 'FontSize', 25, 'FontWeight', 'bold')

% print('barplot', '-depsc')

%% Fixing high cut-off at 8Hz
% axes
xs = 0.5:0.5:9.5; ys = 1:0.5:10;
LowIdx = find(xs == 7.5);
HighIdx = find(ys == 8);

% Mean Acc across subjects
Only8hz_Acc = Att(:,1:LowIdx,HighIdx)*100;
Only8hz_AveAcc = mean(Att(:,1:LowIdx,HighIdx),1)*100;
sems = std(Only8hz_Acc)/sqrt(length(Only8hz_Acc));

% Flotting
figure(5)
clf
X = categorical(string(xs(1:LowIdx)));
bar(X, Only8hz_AveAcc); hold on
errorbar(X, Only8hz_AveAcc', -sems, sems, 'k','linestyle','none', 'lineWidth', 1);
ylim([50, 85])










