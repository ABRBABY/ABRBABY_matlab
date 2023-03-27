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
%OPTIONS_rfe.baseline = [-199, 0] ; 
%OPTIONS_rfe.win_of_interest = [-0.2, 0.5] ; 
OPTIONS_rfe.conditions = {'STD','DEV1','DEV2'} ; 
OPTIONS_rfe.eeg_elec = 1:16 ; 
OPTIONS_rfe.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_rfe.varhistory = 'EEG.history_rfe' ; 
suffix_rfe = '_reref_filtered_epoched_RFE' ; 

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rfe, count_rfe]= test_existance_of_params_in_db(OPTIONS_rfe, suffix_rfe) ; 

% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS_rfe,flag_sub_to_create_rfe, count_rfe, suffix_rfe) ;

%% ------------------- Preprocess : Select trials per condition and reject BAD trials 
OPTIONS_rej.indir = indir ;                             % directory path
OPTIONS_rej.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_rej.rej_high = 150 ;                            % 150 infants; 120 adults     
OPTIONS_rej.bloc = repelem(1:30,30) ;                   % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
suffix_rej = '_REJ' ;
RFE_num = '_reref_filtered_epoched_RFE1' ;              % set of RFE parameters to use for this step
OPTIONS_rej.varhistory = 'EEG.history_rej' ;            % indicates index of rfe set of parameters to use

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rej, count_rej]= test_existance_of_params_in_db(OPTIONS_rej, suffix_rej) ; 

% When current REJ set of parameter as been run on another RFE set or
% parameters, needs to check for the existance of RFE_REJ combination
if sum(flag_sub_to_create_rej) == 0
    [flag_sub_to_create_rej]= test_existance_of_combination(OPTIONS_rej,flag_sub_to_create_rej, count_rej,RFE_num,suffix_rej) ; 
end

% Reject bad trials and save new .set file
[preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'balanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 
[preproc_filenames_unbalanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'unbalanced', flag_sub_to_create_rej, count_rej, suffix_rej,RFE_num) ; 

%% ------------------- Display results at individual level
OPTIONS_disp.params = 'RFE1_REJ1';                            % option of preprocess to consider
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
OPTIONS_disp.indir = indir ;                                  % directory path of files to process
OPTIONS_disp.diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
OPTIONS_disp.plot_dir = plot_dir ;                            % path to save png files of plots
OPTIONS_disp.balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
OPTIONS_disp.ylim = [-20,20] ;                                % limits of y axis

% Display one participant results 
subjects_to_process = {'DVL_005_T18'} ;
% subjects_to_process = get_all_subjects(indir) ;

display_individual_subjects(subjects_to_process, OPTIONS_disp) ; 


%% ------------------- Display results at GROUP level
OPTIONS_disp_contrast.params = 'RFE1_REJ1';                            % option of preprocess to consider
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
