% ERPs Pipeline analysisfor ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
% custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/EEG_data_revised_by_participant_rejA'; 
custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data'; 
% custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data');
% indir = fullfile(custom_path,'FFR_65rej_epoch');

%Get list of subjects in indir
list_subjects = get_subjects(indir,[]);

plot_dir = fullfile(custom_path, 'plot_dir');

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path,~] = get_custom_path();

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : filter, reref, epoch, set chan positions
OPTIONS_rfe.indir = indir ;
OPTIONS_rfe.hp = 1;                                         % high-pass (Hz) (APICE)
OPTIONS_rfe.lp = 30;                                        % low-pass (Hz) (APICE) 
OPTIONS_rfe.mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; 
OPTIONS_rfe.trig = {'Erg1'};                                % Ref and trigger channels 
OPTIONS_rfe.baseline = [-99, 0] ; 
OPTIONS_rfe.win_of_interest = [-0.1, 0.5] ; 
% OPTIONS_rfe.baseline = [-199, 0] ; 
% OPTIONS_rfe.win_of_interest = [-0.2, 0.5] ; 
OPTIONS_rfe.conditions = {'STD','DEV1','DEV2'} ; 
OPTIONS_rfe.eeg_elec = 1:16 ; 
OPTIONS_rfe.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_rfe.varhistory = 'EEG.history_rfe' ; 
OPTIONS_rfe.analysis = 'ERP';

suffix_rfe = strcat('_',OPTIONS_rfe.analysis,'_stepA') ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rfe, count_rfe]= test_existance_of_params_in_db(OPTIONS_rfe, suffix_rfe, '') ; 

%Subjects to process : when whant to choose
subj_to_process = get_subjects(indir,[]);

flag_sub_to_create_rfe = (contains(list_subjects,subj_to_process))';

% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch(ALLEEG, OPTIONS_rfe,flag_sub_to_create_rfe, count_rfe, suffix_rfe) ;

%% ------------------- Preprocess : Select trials per condition and reject BAD trials
OPTIONS_rej.indir = indir ;                             % directory path
OPTIONS_rej.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_rej.rej_high = 150 ;                            % 150 infants; 120 adults     
OPTIONS_rej.bloc = repelem(1:30,30) ;                   % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
OPTIONS_rej.varhistory = 'EEG.history_rej' ;            % indicates index of rfe set of parameters to use
OPTIONS_rej.analysis = 'ERP';

suffix_rej = '_stepB' ;

set_of_param = '1' ; 

RFE_num = strcat(suffix_rfe,set_of_param) ;              % set of RFE parameters to use for this step
RFE_test_existance = strcat('_stepA',set_of_param);

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rej, count_rej]= test_existance_of_params_in_db(OPTIONS_rej, suffix_rej, RFE_test_existance) ; 

%Subjects to process : when whant to choose one subject otherwise comment
%the following line 
% subj_to_process = {'DVL_046_T24'}  ;
flag_sub_to_create_rej = (contains(list_subjects,subj_to_process))';

% Reject bad trials and save new .set file
[preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'balanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 
[preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'unbalanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 


%% ------------------- Manual rejection of bad trials
% OPTIONS_rman.indir = indir ;                             % directory path
OPTIONS_rman.indir = 'E:\preprocessed_data_EEG\RFE1_REJ1'  ;
OPTIONS_rman.suffix_rman = '_rman' ;
OPTIONS_rman.RFE_num = '_stepA1' ;              % set of RFE parameters to use for this step
OPTIONS_rman.REJ_num = '_stepB1' ;

%Get all subjects
subj_to_rman = get_all_subjects(indir) ; 
% subj_to_rman = {'DVL_046_T24'} ;

% Reject bad trials and save new .set file
[preproc_filenames_balanced] = reject_bad_trials_manual(ALLEEG, OPTIONS_rman, 'balanced', subj_to_rman) ; 
[preproc_filenames_balanced] = reject_bad_trials_manual(ALLEEG, OPTIONS_rman, 'unbalanced', subj_to_rman) ; 

%% ------------------- Create grand averages per condition per subject

OPTIONS_average.indir = 'E:\EEG_ANALYSES\to_gd_avg';                                  % directory path of files to process
% spl = split(OPTIONS_average.indir, '\') ;
% OPTIONS_average.param = num2str(cell2mat(spl(end))) ;
OPTIONS_average.param = 'RFE1_REJ1' ;
OPTIONS_average.opt_balance = 'unbalanced' ;
OPTIONS_average.conditions = {'STD1', 'DEV1', 'DEV2'} ;
OPTIONS_average.srate = 256 ; 
OPTIONS_average.keyword = 'gd_avg' ; 
compute_and_save_grand_averages(ALLEEG, OPTIONS_average) ;

%% ------------------- Display results at individual level
OPTIONS_disp.params = 'RFE1_REJ1';                            % option of preprocess to consider
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_disp.elec_subset = {'F3','Fz','F4','Fp1','Fp2','T7','T8','O1';'C3','Cz','C4','Oz','O2','P3','Pz','P4'};   % electrodes to display
OPTIONS_disp.indir = indir ;                                  % directory path of files to process
OPTIONS_disp.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_disp.balance_STD = 'unbalanced';                        % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp.ylim = [-5,5] ;                                % limits of y axis
OPTIONS_disp.savefigs = 1 ; 

if OPTIONS_disp.savefigs ==1 ; create_plot_dirs_if_does_not_exist(plot_dir); end 

% Display one participant results  
subjects_to_process = {'DVL_046_T24'} ;
% subjects_to_process = get_subjects(indir,[]) ;

display_individual_subjects(subjects_to_process, OPTIONS_disp) ; 


%% ------------------- Display results at GROUP level
OPTIONS_disp_contrast.params = 'RFE1_REJ1';                            % option of preprocess to consider
OPTIONS_disp_contrast.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_disp_contrast.indir = indir ;                                  % directory path of files to process
OPTIONS_disp_contrast.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp_contrast.balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp_contrast.ylim = [-5, 5] ;                                 % limits of y axis
OPTIONS_disp_contrast.png_folder = fullfile(plot_dir,'png_folder') ;                          % path to save png files of plots
OPTIONS_disp_contrast.svg_folder = strrep(OPTIONS_disp_contrast.png_folder,'png','svg') ;
OPTIONS_disp_contrast.fig_folder = strrep(OPTIONS_disp_contrast.png_folder,'png','fig') ;      % path to save fig files of plots
OPTIONS_disp_contrast.writecsv = 0 ;
OPTIONS_disp_contrast.savefigs = 0 ; 

OPTIONS_disp_contrast.indir = 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ;
indir = 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ;

if OPTIONS_disp_contrast.savefigs ==1 ; create_plot_dirs_if_does_not_exist(plot_dir); end 

% subjects_to_process_grp1 = {'DVL_013_T10','DVL_005_T18'} ;
% subjects_to_process_grp2 = {'DVL_013_T10','DVL_005_T18'} ;

OPTIONS.suffix = {'_T6','_T8','_T10'} ;
% OPTIONS.suffix = {'_T6','_T8','_T10','_T18','_T24'} ;
subjects_to_process_grp1 = get_subjects(indir,OPTIONS) ;

% OPTIONS.suffix = {'_T18','_T24'} ;
subjects_to_process_grp2 = get_subjects(indir,OPTIONS) ;
% subjects_to_process_grp2 = subjects_to_process_grp1 ;

% Remove subjects based on number of trial rejected 
thresh = 0.33; %(20 DEV kept in each condtion)
% subjects_to_process_grp1 = filter_subjects_based_rejection(subjects_to_process_grp1, thresh, OPTIONS_disp_contrast) ;
% subjects_to_process_grp2 = filter_subjects_based_rejection(subjects_to_process_grp2, thresh, OPTIONS_disp_contrast) ;

% Display results
display_group_comparison(subjects_to_process_grp1, subjects_to_process_grp2, OPTIONS_disp_contrast)

%% ------------------- MMN search
OPTIONS_mmn.params = 'RFE1_REJ1';                            % option of preprocess to consider
OPTIONS_mmn.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_mmn.indir = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/EEG_data_revised_by_participant_rejA' ;                                  % directory path of files to process
OPTIONS_mmn.indir = 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ;
OPTIONS_mmn.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_mmn.balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
OPTIONS_mmn.ylim = [-15,15] ;                                % limits of y axis
OPTIONS_mmn.savefigs = 0 ; 
OPTIONS_mmn.conditions = {'DEV1','DEV2','STD1'};
OPTIONS_mmn.disp = 0 ;                                       % 1 if want to display local peak figure, 0 otherwise
OPTIONS_mmn.auc_delta = 5 ;                                  % time window to compute auc around peak
OPTIONS_mmn.file = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Analyses\mmn_participants_ok80.csv' ;
% Create folder for plots if doesn't exist
if OPTIONS_mmn.savefigs ==1 ; create_plot_dirs_if_does_not_exist(plot_dir); end 

% Lookk for local peak at group then individual level 
% subjects_to_process = {'DVL_012_T24'} ;
% subjects_to_process = get_subjects(OPTIONS_mmn.indir,[]) ;
subjects_to_process = get_subjects(OPTIONS_mmn.indir,OPTIONS_mmn) ;

% Comment ASD : for structuring the code it would be better to have two
% functions : 
% search et detect gr dmoyenne 
% search et detect accross subject 

OPTIONS_mmn.win_gd_mmn = [150, 250] ; 
OPTIONS_mmn.win_mmn = [-120, 120] ; 
OPTIONS_mmn.keyword = 'gd_avg' ; 

[all_lat, all_amp, all_auc] = search_for_mmn(subjects_to_process, OPTIONS_mmn) ; 

% Display violin plot
figure; subplot(1,3,1);
plot_violin(all_lat, {'r','m'}, '', 'Peak latency', {'COND1', 'COND2'}) ;
hold on ; subplot(1,3,2);
plot_violin(all_amp, {'r','m'}, '', 'Peak amplitude', {'COND1', 'COND2'}) ;
hold on ; subplot(1,3,3);
plot_violin(all_auc, {'r','m'}, '', 'AUC (Area Under the Curve)', {'COND1', 'COND2'}) ;

% figure ; plot(1,all_amp(:,1),'r*') ; hold on ; plot(1,all_amp(:,2),'b*') ; 

%% BRAINSTORM PROC
% Call Brainstorm and populate database 
% indir = 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ;
indir = 'E:\EEG_ANALYSES\to_add_in_bst' ;
ABRBABY_populate_BST_DB_averages('RFE1_REJ1','unbalanced',indir)

