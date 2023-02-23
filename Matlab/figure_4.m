
%%

for i = 1:9
    
    RealT(i) = DATA(i).RealTime_all;
    EMA(i) = DATA(i).EMA_all;
    
end
%% box plot
a = [EMA; RealT]';

boxplot(a, 'ColorGroup', [0,1])
ylim([50 100])
set(gca,'color', 'none', 'Ytick',[50 60 70 80 90 100])
set(gcf, 'color', 'none')


%% individual subject at 3 type

sub_type = [RealT_all; RealT_fix; RealT_swi]';
[a,c] = sort(RealT_all, 'descend');

for i =1:9
    idx = c(i);
    realT_a(i) = RealT_all(idx);
    realT_f(i) = RealT_fix(idx);
    realT_s(i) = RealT_swi(idx);
end

sub_type = [realT_a; realT_f; realT_s]';

% bar plot
b = bar( sub_type, 'FaceColor', 'flat' ); hold on

% b.CData(2,:) = [0.8500 0.3250 0.0980];
% b.CData(3,:) = [0.9290 0.6940 0.1250];
set(gcf, 'color', 'none')
set(gca,'color','none', 'Ytick',[50 60 70 80 90 100]);
ylim([50 100])
yline(52, '--', 'lineWidth', 1);

%% figure-c / bar & individual

data = [RealT_all_mean; EMA_all_mean];w

% bar plot
b = bar([RealT_all_mean; EMA_all_mean], 'FaceColor', 'flat' ); hold on
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [0.4660 0.6740 0.1880]; %[0.9290 0.6940 0.1250];

% STD
plot(1*ones(1,2), data(1)+[-RealT_all_std, RealT_all_std], 'color', 'k', 'LineWidth', 3)
plot(2*ones(1,2), data(2)+[-EMA_all_std, EMA_all_std], 'color', 'k', 'LineWidth', 3)
% individual
plot([RealT_all; EMA_all], '--o', 'color', [0.5 0.5 0.5] ,'MarkerFaceColor',[1,1,1],'LineWidth', 1);


set(gcf, 'color', 'none')
set(gca,'color','none', 'Ytick',[50 60 70 80 90 100]);
ylim([50 100])
yline(52, '--', 'lineWidth', 1);






