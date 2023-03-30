% FFR Pipeline analysis for ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
%custom_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
 custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data') ;
plot_dir = fullfile(custom_path, 'png_folder');

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path, BT_toolbox] = get_custom_path();

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : reref, epoch, set chan positions, reject BAD trials 
OPTIONS_rerbt.indir = indir;
OPTIONS_rerbt.mastos = {'Lmon','Rmon','MASTOG','MASTOD'};   %Labels of the mastoids electrodes
OPTIONS_rerbt.trig = {'Erg1'};                              %Label of trigger channel
OPTIONS_rerbt.abr= {'Left','Right'};                        %Label of abr channels for formula 
OPTIONS_rerbt.baseline = [-39, 0] ;                         %Baseline
OPTIONS_rerbt.win_of_interest = [-0.04, 0.2] ;              %Epoching window
OPTIONS_rerbt.eeg_elec = 1:16 ;                             %Cortical electrodes (to get cortical FFRs)
OPTIONS_rerbt.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_rerbt.rej_low = -45; %initial value = -45                               
OPTIONS_rerbt.rej_high = 45; %initial value = 45
OPTIONS_rerbt.bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]
OPTIONS_rerbt.varhistory = 'EEG.history_rerbt' ;
suffix_rerbt = '_reref_epoched_FFR_RERBT';

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rerbt, count_rerbt]= test_existance_of_params_in_db(OPTIONS_rerbt, suffix_rerbt) ; 

%Reref data, compute FFR formula, epoch, reject bad trials and produce
%report
if sum(flag_sub_to_create_rerbt)~=0
    [preproc_filenames] = reref_epoch_rej_ffr(ALLEEG, OPTIONS_rerbt,flag_sub_to_create_rerbt, count_rerbt,suffix_rerbt) ;
end
%%/!\ some improvement to make -> add possibility to compute FFR on cortical electrodes %%

%% ------------------- Preprocess : Filter and Prepare input for BTtoolbox
OPTIONS_fbt.indir = indir ;
OPTIONS_fbt.hp = 30 ;                          % high-pass (Hz) initial value = 80
OPTIONS_fbt.lp = 3000 ;                        % low-pass (Hz) initial value = 3000
OPTIONS_fbt.bt_toolbox = BT_toolbox ; 
OPTIONS_fbt.varhistory = 'EEG.history_fbt' ;
tube_length = 0.27  ; 
propag_sound = 340 ; 
suffix_fbt = '_FBT' ;
RERBT_num = 1 ;                                %Set of rerbt parameters to use for filtering

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_fbt, count_fbt]= test_existance_of_params_in_db(OPTIONS_fbt, suffix_fbt) ; 

%Filter epoched data and prepare input for brainstem toolbox
if sum(flag_sub_to_create_fbt)~=0
    [preproc_filt_filenames] = filter_and_prepare_input_brainstem(ALLEEG, OPTIONS_fbt,tube_length, propag_sound,flag_sub_to_create_fbt, count_fbt,suffix_fbt, RERBT_num);
end

%% ------------------- Display : 

% Display one participant results 
subjects_to_process = {'DVL_003_T10'} ;
%subjects_to_process = get_all_subjects(indir) ;

OPTIONS_disp.params = 'RERBT1_FBT1'; 
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
OPTIONS_disp.indir = indir ; 
OPTIONS_disp.plot_dir = plot_dir ; 
OPTIONS_disp.ylim = [-0.5, 0.5] ; 
OPTIONS_disp.fs = 16384 ; 

display_individual_subjects_FFR(subjects_to_process, OPTIONS_disp) ; 

% Display group result 
%OPTIONS_group.groups_labels = {};

%display_group_comparison(subjects_to_process, OPTIONS_group)
%%%TODO !!
%% -------------------Reject bad participants and compute analyses

% Reject bad participants based on number of trials rejected
OPTIONS_rej.indir = indir ;
OPTIONS_rej.threshold = 500 ;                                  % minimum number of artifact-free trials to keep a participant
OPTIONS_rej.suffix_csv = '_infos_trials_low_-45_high_45' ;     % suffix for CVS file containing trial rejection info
OPTIONS_rej.param = '_RERBT1';                                 % suffix for set of parameters to process
OPTIONS_rej.visu = 1 ;                                         % 1 to display rejection rates, otherwise 0

all_subjects = get_all_subjects(indir) ;
[subjects_to_analyse] = reject_participants_FFR(all_subjects, OPTIONS_rej) ;

% %%
% % Reject bad participants based neural lag
% OPTIONS_rej.indir = indir ;
% OPTIONS_rej.threshold = 3 ;                                    % minimum value of neural lag
% %OPTIONS_rej.suffix_csv = '_infos_trials_low_-45_high_45' ;     % suffix for CVS file containing trial rejection info
% OPTIONS_rej.param = '_RERBT1';                                 % suffix for set of parameters to process
% OPTIONS_rej.visu = 0 ;                                         % 1 to display rejection rates, otherwise 0
% 
% [subjects_to_analyse] = reject_participants_FFR(subjects_to_analyse, OPTIONS_rej) ;

OPTIONS_analysis.indir = indir ;
OPTIONS_analysis.param = '_RERBT1_FBT1';
grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};
OPTIONS_analysis.groups = {grpA, grpB} ;

%OPTIONS_analysis.
%OPTIONS_analysis.


%Run FFR analysis only on kept subjects
FFR_analysis(subjects_to_analyse,OPTIONS_analysis);
