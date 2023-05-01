% ERPs Pipeline analysisfor ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
%custom_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data') ;
plot_dir = fullfile(custom_path, 'png_folder');
 
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
suffix_rfe = '_reref_filtered_epoched_RFE' ;
interpol.subj = {'DVL_021_T18','DVL_045_T8'} ;
interpol.chan = {{'F3';'O2'},{'Cz'}} ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rfe, count_rfe]= test_existance_of_params_in_db(OPTIONS_rfe, suffix_rfe) ; 

% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS_rfe,flag_sub_to_create_rfe, count_rfe, suffix_rfe, interpol) ;

%% Do it again with new parameters of baseline (rfe2)
OPTIONS_rfe.baseline = [-199, 0] ; 
OPTIONS_rfe.win_of_interest = [-0.2, 0.5] ; 
count_rfe = 2 ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rfe, count_rfe]= test_existance_of_params_in_db(OPTIONS_rfe, suffix_rfe, '') ; 

% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS_rfe,flag_sub_to_create_rfe, count_rfe, suffix_rfe) ;

%% ------------------- Preprocess : Select trials per condition and reject BAD trials 
OPTIONS_rej.indir = indir ;                             % directory path
OPTIONS_rej.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_rej.rej_high = 150 ;                            % 150 infants; 120 adults     
OPTIONS_rej.bloc = repelem(1:30,30) ;                   % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
OPTIONS_rej.varhistory = 'EEG.history_rej' ;            % indicates index of rfe set of parameters to use
suffix_rej = '_REJ' ;
RFE_num = '_reref_filtered_epoched_RFE2' ;              % set of RFE parameters to use for this step
RFE_test_existance = RFE_num(24:28);

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rej, count_rej]= test_existance_of_params_in_db(OPTIONS_rej, suffix_rej, RFE_test_existance) ; 

% Reject bad trials and save new .set file
[preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'balanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 
% [preproc_filenames_unbalanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'unbalanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 

%% ------------------- Display results at individual level
OPTIONS_disp.params = 'RFE2_REJ1';                            % option of preprocess to consider
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% OPTIONS_disp.elec_subset = {'F3','Fz','F4','Fp1','Fp2','T7','T8','O1';'C3','Cz','C4','Oz','O2','P3','Pz','P4'};   % electrodes to display
OPTIONS_disp.indir = indir ;                                  % directory path of files to process
OPTIONS_disp.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_disp.balance_STD = 'balanced';                        % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp.ylim = [-15,15] ;                                % limits of y axis

% Display one participant results 
% subjects_to_process = {'DVL_003_T10', 'DVL_003_T6', 'DVL_007_T8', 'DVL_008_T10', 'DVL_018_T8', 'DVL_029_T10', 'DVL_032_T10', 'DVL_021_T18'} ;
%subjects_to_process = {'DVL_004_T10','DVL_004_T8','DVL_006_T10','DVL_007_T10','DVL_011_T10','DVL_012_T10','DVL_013_T10','DVL_013_T8','DVL_018_T10','DVL_018_T6','DVL_024_T6','DVL_030_T10','DVL_037_T6','DVL_037_T8'} ;
subjects_to_process = {'DVL_046_T18'} ;
% subjects_to_process = get_all_subjects(indir) ;

display_individual_subjects(subjects_to_process, OPTIONS_disp) ; 


%% ------------------- Display results at GROUP level
OPTIONS_disp_contrast.params = 'RFE2_REJ1';                            % option of preprocess to consider
OPTIONS_disp_contrast.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
OPTIONS_disp_contrast.indir = indir ;                                  % directory path of files to process
OPTIONS_disp_contrast.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp_contrast.balance_STD = 'balanced';                      % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp_contrast.ylim = [-10, 10] ;                               % limits of y axis
OPTIONS_disp_contrast.png_folder = plot_dir ;                          % path to save png files of plots
OPTIONS_disp_contrast.svg_folder = strrep(plot_dir,'png','svg') ;      % path to save svg files of plots

% subjects_to_process_grp1 = {'DVL_013_T10','DVL_005_T18'} ;
% subjects_to_process_grp2 = {'DVL_013_T10','DVL_005_T18'} ;

suffix1 = {'_T3','_T6','_T8','_T10'} ;
suffix2  = {'_T18','_T24'} ;

subjects_to_process_grp1 = get_subjects_by_suffix(indir,suffix1) ;
subjects_to_process_grp2 = get_subjects_by_suffix(indir,suffix2) ;

% Remove subjects based on number of trial rejected 
thresh = 0.20;
subjects_to_process_grp1 = filter_subjects_based_rejection(subjects_to_process_grp1, thresh, OPTIONS_disp_contrast) ;
subjects_to_process_grp2 = filter_subjects_based_rejection(subjects_to_process_grp2, thresh, OPTIONS_disp_contrast) ;

% Display results
display_group_comparison(subjects_to_process_grp1, subjects_to_process_grp2, OPTIONS_disp_contrast)
