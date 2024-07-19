function [] = plot_variable_nb_cond(OPTIONS, flag_sub_to_create, snr, neural_lag) 

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects_to_process = subjects(flag_sub_to_create) ; 

FONTSZ = 12 ; 

figure('Units','normalized','Position',[0,0.4,0.6,0.6]) ; 
hplot = gca; 

for iGrp=1:length(OPTIONS.groups)

    flag_grp = contains(subjects,OPTIONS.groups{iGrp});
    % [hplot,p1,p2,p3,p4] = plot_patch_violin(hplot,snr(flag_grp),OPTIONS.colors{iGrp},count_violin, OPTIONS.groups{iGrp});
    plot(hplot,snr(flag_grp),neural_lag(flag_grp),'*','MarkerFaceColor',OPTIONS.colors{iGrp}, 'MarkerSize',12) ; hold on ; 
    
    % HERE YOU CAN DISPLAY WHATEVER YOU WANT BY GROUP using flag_grp to
    % restrict to the current group in the loop

    % [X,N] = hist(hplot,snr(flag_grp),10); hold on ; 
    % h = bar(N,X);
    % h.FaceColor = OPTIONS.colors{iGrp}; 
    % h.FaceAlpha = 0.2;

    conditions{iGrp} = sprintf('Group %s (n=%d)',cell2mat(strrep(OPTIONS.groups{iGrp},'_','')),sum(flag_grp));

end

% legend([p1 p2 p3 p4],conditions) ;
legend(conditions) ;
title(OPTIONS.title);
grid on 

end

% plot each violin (within the same plot) 
function [hplot,p1,p2,p3,p4] = plot_patch_violin(hplot,data, color, nb_violin, group_name)

    [y1, x1] = hist(data,20) ;
    y1 = smooth(y1,5)';
    y1 = y1./max(y1) ;
    dataR1 = tiedrank(data(:,1))./size(data,1) ;
    IQR1 = data([dsearchn(dataR1, .25) dsearchn(dataR1, .75)]) ;
  
    % main plot is a patch, start with condition 1
    p1 = patch([y1 -y1(end:-1:1)]+nb_violin , [x1 x1(end:-1:1)], color, 'facealpha', .3, 'Parent',hplot) ;
    hold on ;
    % plot other descriptive stats
    p2 = plot([0 0]+nb_violin, IQR1, 'k', 'linew', 3, 'Parent',hplot) ;
    p3 = plot(0+nb_violin, mean(data,1), 'ks', 'markerfacecolor', 'r', 'markersize', 10, 'Parent',hplot) ;
    p4 = plot(0+nb_violin, median(data,1), 'ko', 'markerfacecolor', 'g', 'markersize', 10, 'Parent',hplot) ;
    
    grid on ;

end
