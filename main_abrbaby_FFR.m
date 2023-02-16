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
OPTIONS_rerbt.indir = indir;
OPTIONS_rerbt.mastos = {'Lmon','Rmon','MASTOG','MASTOD'};   %Labels of the mastoids electrodes
OPTIONS_rerbt.trig = {'Erg1'};                                    %Label of trigger channel
OPTIONS_rerbt.abr= {'Left','Right'};                              %Label of abr channels for formula 
OPTIONS_rerbt.baseline = [-39, 0] ;                         %Baseline
OPTIONS_rerbt.win_of_interest = [-0.04, 0.2] ;                    %Epoching window
OPTIONS_rerbt.eeg_elec = 1:16 ;                             %Cortical electrodes (to get cortical FFRs)
OPTIONS_rerbt.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_rerbt.rej_low = -45;                                
OPTIONS_rerbt.rej_high = 45;
bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]
suffix_rerbt = '_reref_epoched_FFR_RERBT';

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rerbt, count_rerbt]= test_existance_of_params_in_db(OPTIONS_rerbt, suffix_rerbt) ; 

%Reref data, compute FFR formula and epoch
[preproc_filenames] = reref_epoch_ffr(ALLEEG, OPTIONS_rerbt,flag_sub_to_create_rerbt, count_rerbt,suffix_rerbt,bloc);

% Reject bad trials and write a report
%EEG = reject_trials_produce_report(out_filenames(jj), find(ismember({EEG.chanlocs.labels},'ABR')), bloc, win_of_interest, rej_low, rej_high,'FFR') ; 
%in ERPS:[~] = reject_trials_produce_report(preproc_filenames, eeg_elec, bloc, win_of_interest, rej_low, rej_high,'') ; 

%% ------------------- Preprocess : Filter 
OPTIONS_f.indir = indir;
OPTIONS_f.hp =80;                       % high-pass (Hz) (APICE)
OPTIONS_f.lp = 3000;                    % low-pass (Hz) (APICE)
OPTIONS_f.bt_toolbox = '\\Filer\home\Invites\hervé\Mes documents\GitHub\ABRBABY_matlab\ToolBox_BrainStem\BT_2013' ;  % bt_toolbox path
OPTIONS_f.RERBT = 1;                    %Set of rerbt parameters to use for filtering
tube_length = 0.27  ; 
propag_sound = 340 ; 
suffix_f = '_filtered_FFR_F' ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_f, count_f]= test_existance_of_params_in_db(OPTIONS_f, suffix_f) ; 

%Filter epoched data and prepare input for brainstem toolbox
[preproc_filt_filenames] = filter_and_prepare_input_brainstem(ALLEEG, OPTIONS_f,tube_length, propag_sound,flag_sub_to_create_f, count_f,suffix_f);

%% ------------------- Display : 
elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
%function to update!!
%display_timeseries_by_condition(preproc_filenames, elec_subset, 'balanced',plot_dir) ; 
