% % Potential things to improve : 
%
% Create an history_abr variable to append to the EEG structure to document
% the operations and implement the overwrite option

%% ------------------- Set environment 
% Variables to enter manually before running the code

name = 'EH';
%name = 'ASD';

[eeglab_path, biosig_installer_path, erplab_path,indir,plot_dir, BT_toolbox] = get_custom_path(name);

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : reref, epoch, set chan positions, reject BAD trials 
OPTIONS.indir = indir;
OPTIONS.mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; OPTIONS.trig = {'Erg1'}; OPTIONS.abr= {'Left','Right'};  % Ref and trigger channels 
OPTIONS.baseline = [-39, 0] ; OPTIONS.win_of_interest = [-0.04, 0.2] ; 
OPTIONS.eeg_elec = 1:16 ; 
OPTIONS.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS.rej_low = -45;
OPTIONS.rej_high = 45;
overwrite = 0 ; % this option allow to overwrite (=1) or not (=0) 
bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]

[preproc_filenames] = reref_epoch_ffr(ALLEEG, OPTIONS, bloc, overwrite);

%% ------------------- Preprocess : Filter 
OPTIONSFFR.indir = indir;
OPTIONSFFR.hp =80; % high-pass (Hz) (APICE)
OPTIONSFFR.lp = 3000; % low-pass (Hz) (APICE) 
OPTIONSFFR.RERBT = 1;
overwrite = 0 ; % this option allow to overwrite (=1) or not (=0) 
tube_length = 0.27  ; 
propag_sound = 340 ; 

[preproc_filt_filenames] = filter_and_prepare_input_brainstem(ALLEEG, EEG, OPTIONSFFR, overwrite, tube_length, propag_sound,BT_toolbox) ; 

%% ------------------- Display : 
elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
%function to update!!
%display_timeseries_by_condition(preproc_filenames, elec_subset, 'balanced',plot_dir) ; 
