function [] = display_group_comparison(subjects_grp1,subjects_grp2, OPTIONS)

%Get grand averages for each condition and subject
%(dimensions = subjects x channels x timepoints)
[grpA.DEV1_avg, grpA.DEV2_avg, grpA.STD1_avg, grpA.STD2_avg, timepoints, labels, xlim] = extract_averages_DEV_STD(subjects_grp1,OPTIONS);
[grpB.DEV1_avg, grpB.DEV2_avg, grpB.STD1_avg, grpB.STD2_avg, ~, ~, ~] = extract_averages_DEV_STD(subjects_grp2,OPTIONS);

grpA_color = [0.2 0.2 1]; %blue
%grpB_color = [0.2 0.7765 0.2]; %green
grpB_color = [0.8902 0 0]; %red

% Get number of prticipants in both group to display in legend 
nb_young = size(grpA.DEV1_avg,1);
nb_old = size(grpB.DEV1_avg,1);

% Get indices of electrodes in subset
[sharedvals,idxA,idxB] = intersect(labels, OPTIONS.elec_subset(:),'stable') ; 
[~,id1,id2] = intersect(OPTIONS.elec_subset(:), sharedvals,'stable') ;
OPTIONS.elec_indices = reshape(idxA(id2),size(OPTIONS.elec_subset)) ;

OPTIONS.vTime = timepoints; 
OPTIONS.xlim = xlim ;  

%% First display : STD
%--------------------------------------------------------------------------

% La permutation (ci-dessous) sert à visualiser la variance au travers du nb de
% participants (on change par rapport à l'appel avec la matrice des essais : les essais se trouve dans la dimension 3)

signals =           {   permute(cat(1,grpA.STD1_avg,grpA.STD2_avg),[2,3,1]),...
                        permute(cat(1,grpB.STD1_avg,grpB.STD2_avg),[2,3,1])} ; 
OPTIONS.color =     {grpA_color,                                                                        grpB_color} ; 
OPTIONS.legend =    {sprintf('Std /DA/ 6-10mo \n n= %d',nb_young),                                      sprintf('Std /DA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title =     sprintf('Response to standard /DA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('group_display_STD_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('group_display_STD_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig, fullfile(OPTIONS.fig_folder,strcat('group_display_STD_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : DEV1
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV1_avg,[2,3,1]),                permute(grpB.DEV1_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('Dev /BA/ 6-10mo \n %d',nb_young),     sprintf('Dev /BA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('Response to deviant /BA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('group_display_DEV1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('group_display_DEV1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig, fullfile(OPTIONS.fig_folder,strcat('group_display_DEV1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : DEV2
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV2_avg,[2,3,1]),                permute(grpB.DEV2_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('Dev /GA/ 6-10mo \n %d',nb_young),     sprintf('Dev /GA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('Response to deviant /GA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('group_display_DEV2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('group_display_DEV2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('group_display_DEV2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : MMN DEV1
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV1_avg-grpA.STD1_avg,[2,3,1]),  permute(grpB.DEV1_avg-grpB.STD1_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('MMN /BA/ 6-10mo \n %d',nb_young),     sprintf('MMN /BA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('MMN in response to deviant /BA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('group_display_MMN1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('group_display_MMN1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('group_display_MMN1_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : MMN DEV2
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV2_avg-grpA.STD2_avg,[2,3,1]),                permute(grpB.DEV2_avg-grpB.STD2_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('MMN /GA/ 6-10mo \n %d',nb_young),     sprintf('MMN /GA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('MMN in response to deviant /GA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('group_display_MMN2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('group_display_MMN2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('group_display_MMN2_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupA - DEV1 (/BA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpA,1,OPTIONS,'6-10 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('6_10mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('6_10mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('6_10mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupA - DEV2 (/GA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpA,2,OPTIONS,'6-10 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('6_10mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('6_10mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('6_10mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupA - mean(DEV1+DEV2) (/BA/ and /GA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpA,3,OPTIONS,'6_10 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('6_10mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('6_10mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('6_10mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupB - DEV1 (/BA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpB,1,OPTIONS,'18-24 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('18_24mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('18_24mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('18_24mo_display_STD_DEV_MMN_BA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupB - DEV2 (/GA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpB,2,OPTIONS,'18-24 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('18_24mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('18_24mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('18_24mo_display_STD_DEV_MMN_GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%% Display : GroupB - mean(DEV1+DEV2) (/BA/ and /GA/)
%--------------------------------------------------------------------------
[fig] = display_one_group(grpB,3,OPTIONS,'18-24 mo'); 

% Save figure in svg + png
print('-dsvg',fullfile(OPTIONS.svg_folder,strcat('18_24mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.png_folder,strcat('18_24mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.png')));
saveas(fig,fullfile(OPTIONS.fig_folder,strcat('18_24mo_display_STD_DEV_MMN_BA-GA_',OPTIONS.balance_STD,'_',OPTIONS.params,'.fig')));

%--------------------------------------------------------------
% FUNCTION that displays STD, DEV and MMN at group level
%--------------------------------------------------------------

function [fig] = display_one_group(grp,cond_num,OPTIONS,grp_label)

cond.label = {'/BA/', '/GA/', '/BA/-/GA/'} ;
cond.colors = {[255,215,0]/255,[255,130,0]/255,[255,0,255]/255} ;

mean_STD_grp = (grp.STD1_avg + grp.STD2_avg) / 2 ; 
mean_DEV_grp = (grp.DEV1_avg + grp.DEV2_avg) / 2 ; 
mean_MMN_grp = mean_DEV_grp - mean_STD_grp ; 

STD_options = {grp.STD1_avg, grp.STD2_avg,mean_STD_grp} ;
DEV_options = {grp.DEV1_avg, grp.DEV2_avg,mean_DEV_grp} ;
MMN_options = {grp.DEV1_avg-grp.STD1_avg, grp.DEV2_avg-grp.STD2_avg,mean_MMN_grp} ;

color_dev = cond.colors{cond_num} ;
STD = STD_options{cond_num} ;
DEV = DEV_options{cond_num} ;
MMN = MMN_options{cond_num} ;

signals =           {   permute(STD,[2,3,1]),...
                        permute(DEV,[2,3,1]),...
                        permute(MMN,[2,3,1])};
                    
OPTIONS.color =     {   [147,112,219]/255,...   % Purple
                        color_dev,...           % Yellow, Orange or Pink
                        [0,0,0]} ;              % Black
                    
OPTIONS.legend =    {   sprintf('STD /DA/ %s \n %d',grp_label,size(STD,1)),...
                        sprintf('DEV %s %s \n n= %d',cond.label{cond_num},grp_label,size(DEV,1)),...
                        sprintf('MMN to %s %s \n n= %d',cond.label{cond_num},grp_label,size(MMN,1))} ; 
                    
OPTIONS.title = 	sprintf('MMN %s for group %s | %s number of STD',cond.label{cond_num},grp_label,OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
[fig] = plot_electrodes_subset(signals,OPTIONS) ; 

function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj, timepoints, labels, xlim] = extract_averages_DEV_STD(subjects,OPTIONS)

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject

    %dimensions XXX.data input (for 1 subject) = channels x timepoints x trials
    %dimensions datasets (.set) output (average) =  subjects x channels x timepoints
        
    [STD1_subj(loopnum,:,:), timepoints, labels, xlim]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*STD*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 
    
    [STD2_subj(loopnum,:,:),~,~,~]   = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*STD*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 

    [DEV1_subj(loopnum,:,:),~,~,~]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*DEV1*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 
    
    [DEV2_subj(loopnum,:,:),~,~,~]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*DEV2*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 
end



function [data, timepoints, labels, xlim] = get_data(datafile) 

    FileToLoad = dir(datafile) ; 
    if size(FileToLoad,1) > 1
        rman = find(contains({FileToLoad.name}, '_gd_avg')) ;
        FileToLoad = FileToLoad(rman,1) ;
    end
    if size(FileToLoad,1) > 1
        rman = find(contains({FileToLoad.name}, 'rman')) ;
        FileToLoad = FileToLoad(rman,1) ;
    end
    EEG = pop_loadset(FileToLoad.name,FileToLoad.folder) ; 
    data = mean(EEG.data,3) ;
    
    timepoints = EEG.times;
    labels = {EEG.chanlocs.labels};
    xlim = [EEG.xmin, EEG.xmax]*1000 ;  


end

end
