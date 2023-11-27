function [] = plot_violin(data, colors, transparency, tit, conditions)
% ERPs sanity check script -
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% INPUTS
% data = data to plot, format : datapoints x number of conditions (n)
% colors = colors for plots, format : {'color 1, 'color 2',..., 'color n'}
% transparency =   ;
% tit = title of the plot or each subplot, format:  {'title 1, 'title 2',..., 'title n'}
% conditions = conditions to compare (n conditions), format : {'condition 1, 'condition 2',..., 'condition n'}

% get data distribution based on histogram
[y1, x1] = hist(squeeze(data(:,1)),20) ;
[y2, x2] = hist(squeeze(data(:,2)),20) ;
y1 = smooth(y1,5)';
y2 = smooth(y2,5)';
y1 = y1./max(y1) ;
y2 = y2./max(y2) ;

% compute inter-quartile range
dataR1 = tiedrank(data(:,1))./size(data,1) ;
dataR2 = tiedrank(data(:,2))./size(data,1) ;
IQR1 = data([dsearchn(dataR1, .25) dsearchn(dataR1, .75)]) ;
IQR2 = data([dsearchn(dataR2, .25) dsearchn(dataR2, .75)]) ;

% main plot is a patch, start with condition 1
p1 = patch([y1 -y1(end:-1:1)] , [x1 x1(end:-1:1)], colors{1}, 'facealpha', .3) ;
hold on ;
% plot other descriptive stats
p3 = plot([0 0], IQR1, 'k', 'linew', 3) ;
p4 = plot(0, mean(data(:,1)), 'ks', 'markerfacecolor', 'r', 'markersize', 10) ;
p5 = plot(0, median(data(:,1)), 'ko', 'markerfacecolor', 'g', 'markersize', 10) ;

% do the same for condition 2
p2 = patch([y2 -y2(end:-1:1)]+3 , [x2 x2(end:-1:1)], colors{2}, 'facealpha', .3) ;
% plot other descriptive stats
plot([3 3], IQR2, 'k', 'linew', 3) ;
plot(3, mean(data(:,2)), 'ks', 'markerfacecolor', 'r', 'markersize', 10) ;
plot(3, median(data(:,2)), 'ko', 'markerfacecolor', 'g', 'markersize', 10) ;

% add legend and title
conditions{3} = 'IQR' ;
conditions{4} = 'Mean' ;
conditions{5} = 'Median' ;
legend([p1 p2 p3 p4 p5],conditions) ;
title(tit) ;
grid on ;


%% Original script

% all(1,:,:) = data.all_lat ;
% all(2,:,:) = data.all_amp ;
% all(3,:,:) = data.all_auc ;
%
% tit = {'peak latency','peak amplitude','auc'} ;
% leg = {'COND1', 'COND2'} ;
% figure ;
%
% for dd=1:3
%     [y1, x1] = hist(squeeze(all(dd,:,1)),20) ;
%     [y2, x2] = hist(squeeze(all(dd,:,2)),20) ;
%     y1 = smooth(y1,5)';
%     y2 = smooth(y2,5)';
%     y1 = y1./max(y1) ;
%     y2 = y2./max(y2) ;
% %     dataR1 = tiedrank(all(dd,1,:))./size(data_avg,2) ;
% %     dataR2 = tiedrank(all(dd,2,:))./size(data_avg,2) ;
% %
%     subplot(1,3,dd) ;
%     patch([y1 -y1(end:-1:1)] , [x1 x1(end:-1:1)], 'm', 'facealpha', .3) ;
%     hold on ;
%
%     patch([y2 -y2(end:-1:1)]+3 , [x2 x2(end:-1:1)], 'r', 'facealpha', .3) ;
%     legend(leg) ;
%     title(tit(dd)) ;
%     grid on ;
%
% end

end
