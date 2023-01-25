

%% ------------------- Set environment 
% Variables to enter manually before running the code
% name = 'EH';
name = 'ASD';

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path,indir,plot_dir,~] = get_custom_path(name);

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : filter, reref, epoch, set chan positions
OPTIONS.indir = indir ;
OPTIONS.hp = 1; % high-pass (Hz) (APICE)
OPTIONS.lp = 30; % low-pass (Hz) (APICE) 
OPTIONS.mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; OPTIONS.trig = {'Erg1'}; % Ref and trigger channels 
OPTIONS.baseline = [-99, 0] ; OPTIONS.win_of_interest = [-0.1, 0.5] ; 
OPTIONS.conditions = {'STD','DEV1','DEV2'} ; 
OPTIONS.eeg_elec = 1:16 ; 
OPTIONS.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS.overwrite = 1 ; % this option allow to overwrite (=1) or not (=0) 
[preproc_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS);

%% ------------------- Preprocess : Reject BAD trials 
rej_low = -150; %150 infants; 120 adults
rej_high = 150; %150 infants; 120 adults
bloc = repelem(1:30,30) ; % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]

% Reject bad trials and save new .set file
select_and_save_trials_per_condition(ALLEEG, preproc_filenames, eeg_elec, win_of_interest, rej_low, rej_high, 'balanced') ; 
select_and_save_trials_per_condition(ALLEEG, preproc_filenames, eeg_elec, win_of_interest, rej_low, rej_high, 'unbalanced') ; 

% Write csv file directly into the subject dir
[~] = reject_trials_produce_report(preproc_filenames, eeg_elec, bloc, win_of_interest, rej_low, rej_high,'') ; 


%% ------------------- Display : 
elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
%plot_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder' ; 
display_timeseries_by_condition(preproc_filenames, elec_subset, 'balanced',plot_dir) ; 
