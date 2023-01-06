%% ===== PLOT HALO PATCH =====
% Adpated from Braintorm function PlotHaloPatch
% 2021/07/23 :  Anne-Sophie Dubarry
function hPatch = plotHaloPatchMAD(hAxes, vTime, d, C)

% Compute median absolute deviation (MAD)
% Note: mad = median(abs(X - median(X)))
mad = median(abs(d-repmat(median(d'),size(d,2),1)')');

% Set mean, low and high halo borders
avg=mean(d,2) ;
Lhi=mean(d,2) + mad' ; 
Llow=mean(d,2) - mad' ; 
    
% Eliminate Negative & Zero Data
% idx = Lhi>0 & Llow>0 & avg>0; 

% Draw halo patch 
hPatch = patch([vTime fliplr(vTime)], [Llow'  fliplr(Lhi')], [0.6 0.6 0.6],...
        'FaceAlpha', 0.1, ...
        'FaceColor', C/255, ...
        'EdgeColor', 'none', ...
        'Parent',    hAxes);

% hPatch = patch([vTime fliplr(vTime)], [Llow'  fliplr(Lhi')], [0.6 0.6 0.6],...
%         'FaceAlpha', 0.3, ...
%         'EdgeColor', 'none', ...
%         'Parent',    hAxes);

    % Skip the name of the previous plot from the legend
hPatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

% % Plot timeserie
% hold on ; plot(hAxes, vTime,avg, 'color',C/255); hold off ; 

end
