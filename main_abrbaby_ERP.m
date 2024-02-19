% ERPs Pipeline analysisfor ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project
% To run this script, data should be organized as follows: 

% - One level-1 folder (custom_path) containing data folder (indir) and
% plot folder (plot_dir)
% 
% - One level-2 folder (indir) containing level-3 folders: for n
% recordings, there should be n level-3 folders (called "participants
% folders")
%
% - In which participant folder, there must be:
%           > a "participantID.bdf" file containing the EEG data 
%           > a "participantID_trials_description.txt" file containing
%           information about events, which is a table with 3 variables:
%           event labels, rejection choice (1 for kept, 0 for rejected),
%           and latency. No variable name needed.
%
% - Participant ID must have the following format: "STUDYCODE_PARTICIPANTNUMBER_AGE"

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
% custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/EEG_data_revised_by_participant_rejA'; 
custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data'; 
% custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data');

%Get list of subjects in indir
list_subjects = get_subjects(indir,[]);

plot_dir = fullfile(custom_path, 'plot_dir');

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path,~] = get_custom_path();

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%%
%% Here for first execution run automatic_trigger_detection for fixing the trigger issues (input : .bdf + ergstim.txt if ERGO was used)
%% This will produce the initial __trial_description.txt file

%% ------------------- Preprocess : filter, reref, epoch, set chan positions
OPTIONS_stepA.indir = indir ;
OPTIONS_stepA.hp = 1;                                         % high-pass (Hz) (APICE)
OPTIONS_stepA.lp = 30;                                        % low-pass (Hz) (APICE) 
OPTIONS_stepA.mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; 
OPTIONS_stepA.trig = {'Erg1'};                                % Ref and trigger channels 
OPTIONS_stepA.baseline = [-99, 0] ; 
OPTIONS_stepA.win_of_interest = [-0.1, 0.5] ; 
% OPTIONS_stepA.baseline = [-199, 0] ; 
% OPTIONS_stepA.win_of_interest = [-0.2, 0.5] ; 
OPTIONS_stepA.conditions = {'STD','DEV1','DEV2'} ; 
OPTIONS_stepA.eeg_elec = 1:16 ; 
OPTIONS_stepA.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_stepA.varhistory = 'EEG.history_stepA' ; 
OPTIONS_stepA.analysis = 'ERP';
OPTIONS.file = fullfile(indir,'force_rerun_participants.csv') ;

% suffix_stepA = strcat('_',OPTIONS_stepA.analysis,'_stepA') ;
suffix_stepA = '_stepA' ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_stepA, count_stepA]= test_existance_of_params_in_db(OPTIONS_stepA, suffix_stepA, '') ; 

%Subjects to process : when whant to choose
if exist(OPTIONS.file,'file')
   subj_to_process  = get_subjects(indir,OPTIONS);
    flag_sub_to_create_stepA = (contains(list_subjects,subj_to_process))';
end

% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch(ALLEEG, OPTIONS_stepA,flag_sub_to_create_stepA, count_stepA, suffix_stepA) ;

%% ------------------- Preprocess : Select trials per condition and reject BAD trials
OPTIONS_stepB.indir = indir ;                             % directory path
OPTIONS_stepB.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_stepB.rej_high = 150 ;                            % 150 infants; 120 adults     
OPTIONS_stepB.bloc = repelem(1:30,30) ;                   % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
OPTIONS_stepB.varhistory = 'EEG.history_stepB' ;            % indicates index of rfe set of parameters to use
OPTIONS_stepB.analysis = 'ERP';
OPTIONS.file = fullfile(indir,'force_rerun_participants.csv') ;

suffix_stepB = '_stepB' ;

stepA_num = '_stepA1' ;              % set of RFE parameters to use for this step

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_stepB, count_stepB]= test_existance_of_params_in_db(OPTIONS_stepB, suffix_stepB, stepA_num) ; 

%Subjects to process : when whant to choose
if exist(OPTIONS.file,'file')
   subj_to_process  = get_subjects(indir,OPTIONS);
   flag_sub_to_create_stepB = (contains(list_subjects,subj_to_process))';
end

% Reject bad trials and save new .set file
% [preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_stepB, 'balanced', flag_sub_to_create_stepB, count_stepB, suffix_stepB,RFE_num) ; 
reject_bad_trials(ALLEEG, OPTIONS_stepB, 'unbalanced', flag_sub_to_create_stepB, count_stepB, suffix_stepB,stepA_num) ; 

%% ------------------- Update trial_description with manual rejection of bad trials
OPTION_rman.manualdir = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/manually_marked';
% OPTION_rman.manualdir = 'E:\EEG_ANALYSES\manually_marked';
OPTION_rman.indir = indir ; 

% Here update the trial_desrption files in the corresponding participants
% folder (found in initial indir) 
% Pick any .set (without parameters constraints) and add a column to the
% trial_description w header 'manual_rej'
add_manual_bad_trial_detection(ALLEEG, OPTION_rman) ;

%% ------------------- Create grand averages per condition per subject
OPTIONS_average.indir = indir;                                  % directory path of files to process
% spl = split(OPTIONS_average.indir, '\') ;
% OPTIONS_average.param = num2str(cell2mat(spl(end))) ;
OPTIONS_average.param = 'stepA1_stepB1' ;
OPTIONS_average.opt_balance = 'unbalanced' ;
OPTIONS_average.conditions = {'STD1', 'DEV1', 'DEV2'} ;
OPTIONS_average.srate = 256 ; 
OPTIONS_average.keyword = 'gd_avg' ; 
OPTIONS_average.analysis = 'ERP';
compute_and_save_grand_averages(ALLEEG, OPTIONS_average) ;

%% ------------------- Display results at individual level
OPTIONS_disp.params = 'stepA1_stepB1' ;                         % option of preprocess to consider
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_disp.elec_subset = {'F3','Fz','F4','Fp1','Fp2','T7','T8','O1';'C3','Cz','C4','Oz','O2','P3','Pz','P4'};   % electrodes to display
OPTIONS_disp.indir = indir ;                                  % directory path of files to process
% OPTIONS_disp.indir = 'E:\EEG_ANALYSES\EEGdata_CF_revised_excluded' ;                                  % directory path of files to process
OPTIONS_disp.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_disp.balance_STD = 'unbalanced';                        % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp.ylim = [-5,5] ;                                % limits of y axis
OPTIONS_disp.savefigs = 1 ; 
OPTIONS_disp.keyword = 'gd_avg' ; 
OPTIONS_disp.analysis = 'ERP';

if OPTIONS_disp.savefigs ==1 ; create_plot_dirs_if_does_not_exist(plot_dir); end 

% Display one participant results  
%subjects_to_process = {'DVL_055_T18'} ;
%subjects_to_process = {'DVL_046_T24'} ;
subjects_to_process = get_subjects(indir,[]) ;

display_individual_subjects(subjects_to_process, OPTIONS_disp) ; 


%% ------------------- Display results at GROUP level
OPTIONS_disp_contrast.params = 'stepA1_stepB1' ;                            % option of preprocess to consider
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
OPTIONS_mmn.params = 'stepA1_stepB1';                            % option of preprocess to consider
OPTIONS_mmn.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_mmn.indir = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/EEG_data_revised_by_participant_rejA' ;                                  % directory path of files to process
OPTIONS_mmn.indir = indir % 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ;
OPTIONS_mmn.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_mmn.balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
OPTIONS_mmn.ylim = [-15,15] ;                                % limits of y axis
OPTIONS_mmn.savefigs = 0 ; 
OPTIONS_mmn.conditions = {'DEV1','DEV2','STD1'};
OPTIONS_mmn.disp = 0 ;                                       % 1 if want to display local peak figure, 0 otherwise
OPTIONS_mmn.auc_delta = 5 ;                                  % time window to compute auc around peak
OPTIONS_mmn.file = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/participants_to_compute.csv' ;
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

