function [] = display_group_comparison(subjects_grp1,subjects_grp2, OPTIONS)

%Get grand averages for each condition and subject
%(dimensions = subjects x channels x timepoints)
[grpA.DEV1_avg, grpA.DEV2_avg, grpA.STD1_avg, grpA.STD2_avg, timepoints, labels, xlim] = extract_averages_DEV_STD(subjects_grp1,OPTIONS);
[grpB.DEV1_avg, grpB.DEV2_avg, grpB.STD1_avg, grpB.STD2_avg, timepoints, labels, xlim] = extract_averages_DEV_STD(subjects_grp2,OPTIONS);

grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.2 0.7765 0.2]; %green

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

% Estelle : la permutation (ci-dessous) sert à visualiser la variance au travers du nb de
% participants (on change par rapport à l'appel avec la marice des essais : les essais se trouve dans la dimension 3)

signals =           {   permute(cat(1,grpA.STD1_avg,grpA.STD2_avg),[2,3,1]),...
                        permute(cat(1,grpB.STD1_avg,grpB.STD2_avg),[2,3,1])} ; 
OPTIONS.color =     {grpA_color,                                                                        grpB_color} ; 
OPTIONS.legend =    {sprintf('Std /DA/ 6-10mo \n n= %d',nb_young),                                      sprintf('Std /DA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title =     sprintf('Response to standard /DA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 

%% Display : DEV1
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV1_avg,[2,3,1]),                permute(grpB.DEV1_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('Dev /BA/ 6-10mo \n %d',nb_young),     sprintf('Dev /BA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('Response to standard /DA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 

%% Display : DEV2
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV2_avg,[2,3,1]),                permute(grpB.DEV2_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('Dev /GA/ 6-10mo \n %d',nb_young),     sprintf('Dev /GA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('Response to standard /GA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 

%% Display : MMN DEV1
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV1_avg-grpA.STD1_avg,[2,3,1]),  permute(grpB.DEV1_avg-grpB.STD1_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('MMN /BA/ 6-10mo \n %d',nb_young),     sprintf('MMN /BA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('MMN in response to deviant /BA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 


%% Display : MMN DEV2
%--------------------------------------------------------------------------
signals =           {permute(grpA.DEV2_avg-grpA.STD2_avg,[2,3,1]),                permute(grpB.DEV2_avg-grpB.STD2_avg,[2,3,1])} ; % Permute allow to plot variability across participants 
OPTIONS.color =     {grpA_color,                                    grpB_color} ; 
OPTIONS.legend =    {sprintf('MMN /GA/ 6-10mo \n %d',nb_young),     sprintf('MMN /GA/ 18-24m \n n= %d',nb_old)} ; 
OPTIONS.title = 	sprintf('MMN in response to deviant /GA/ group comparison %s number of STD',OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 

%% Display : GroupA
%--------------------------------------------------------------------------
display_one_group(grpA,OPTIONS,'6-10 mo'); 

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 

%% Display : GroupB
%--------------------------------------------------------------------------
display_one_group(grpB,OPTIONS,'18-24 mo');

% Ici : Estelle -> todo rajouter 1 ligne pour sauver le svg 


function [] = display_one_group(grp,OPTIONS,grp_label)

%% Display : GroupA
%--------------------------------------------------------------------------
mean_STD_grp = (grp.STD1_avg + grp.STD2_avg) / 2 ; 
mean_DEV_grp = (grp.DEV1_avg + grp.DEV2_avg) / 2 ; 
mean_MMN_grp = mean_DEV_grp - mean_STD_grp ; 

signals =           {   permute(mean_STD_grp,[2,3,1]),...
                        permute(mean_DEV_grp,[2,3,1]),...
                        permute(mean_MMN_grp,[2,3,1])};
                    
OPTIONS.color =     {   [147,112,219]/255,...   % Purple
                        [255,215,0]/255,...     % Yellow
                        [0,0,0]} ;              % Black
                    
OPTIONS.legend =    {   sprintf('STD BA/GA/ %s \n %d',grp_label,size(mean_STD_grp,1)),...
                        sprintf('DEV BA/GA/ %s \n n= %d',grp_label,size(mean_DEV_grp,1)),...
                        sprintf('MMN BA/GA/ %s \n n= %d',grp_label,size(mean_MMN_grp,1))} ; 
                    
OPTIONS.title = 	sprintf('MMN for group %so | %s number of STD',grp_label,OPTIONS.balance_STD) ; 

% Create figure 
figure('Units','normalized','Position',[0,0,1,1]) ; 

% Call visualisation function (grids with electrode subset) 
plot_electrodes_subset(signals,OPTIONS) ; 
