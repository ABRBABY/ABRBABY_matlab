function [] = plot_pitchtrack(OPTIONS, flag_sub_to_create, vTime, data, vTime_stim, data_stim) 

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects_to_process = subjects(flag_sub_to_create) ; 

FONTSZ = 12 ; 

figure('Units','normalized','Position',[0,0.4,0.2,1]) ; 
hplot = gca; 
nPlots = length(OPTIONS.groups); 

% For each group
for iGrp=1:length(OPTIONS.groups)
    flag_grp = contains(subjects_to_process,OPTIONS.groups{iGrp});
    hplots(iGrp) = subplot(nPlots,1,iGrp) ; 
    plot(vTime,mean(data(flag_grp,:),1),'s', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6) ; 
    hold on ; plot(vTime_stim, mean(data_stim(flag_grp,:),1), 'k', 'LineWidth', 2); 
    legend('Response', 'F0') ;
    conditions{iGrp} = sprintf('Group %s (n=%d)',cell2mat(strrep(OPTIONS.groups{iGrp},'_','')),sum(flag_grp));
    title(conditions{iGrp});
    grid on ;
end

ylim(hplots,[80 120]); 
% 
% set(hplot,'XTick',ticklabels, 'XTickLabels',groups_names, 'FontSize',FONTSZ) ; 
% xtickangle(hplot,40) ; 

end