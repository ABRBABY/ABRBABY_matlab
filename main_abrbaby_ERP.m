% ERPs Pipeline analysisfor ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code
% name = 'EH';
% name = 'ASD';

% DATA directory 
custom_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
% custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data';

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
OPTIONS_rej.indir = indir ;                             % options st for 'select trials'
OPTIONS_rej.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_rej.rej_high = 150 ;                            % 150 infants; 120 adults
OPTIONS_rej.RFE = '_reref_filtered_epoched_RFE1' ;      % indicates index of rfe set of parameters to use
OPTIONS_rej.bloc = repelem(1:30,30) ;                              % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
suffix_rej = '_REJ' ;
OPTIONS_rej.varhistory = 'EEG.history_rej' ; 

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rej, count_rej]= test_existance_of_params_in_db(OPTIONS_rej, suffix_rej) ; 
    
% Reject bad trials and save new .set file
[preproc_filenames_balanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'balanced', flag_sub_to_create_rej, count_rej, suffix_rej) ; 
[preproc_filenames_unbalanced] = reject_bad_trials(ALLEEG, OPTIONS_rej, 'unbalanced', flag_sub_to_create_rej, count_rej, suffix_rej) ; 

%%%%%%%%%%% RUNS UNTIL HERE on 2023/02/16

%% ------------------- Display results

% Display one participant results 
subjects_to_process = {'DVL_013_T10','DVL_005_T18'} ;
% subjects_to_process = get_all_subjects(indir) ;

OPTIONS_disp.params = 'RFE1_REJ3'; 
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
OPTIONS_disp.indir = indir ; 
OPTIONS_disp.diff_display = 1 ; 
OPTIONS_disp.plot_dir = plot_dir ; 
OPTIONS_disp.balance_STD = 'unbalanced'; 
OPTIONS_disp.ylim = [-20,20] ; 

display_individual_subjects(subjects_to_process, OPTIONS_disp) ; 

% Display group result 
