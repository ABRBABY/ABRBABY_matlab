% FFR Pipeline analysis for ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
custom_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
%  custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data') ;
plot_dir = fullfile(custom_path, 'png_folder');

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path, BT_toolbox] = get_custom_path();

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : reref, set chan positions, filter
OPTIONS_stepA.indir = indir;
OPTIONS_stepA.mastos = {'Lmon','Rmon','MASTOG','MASTOD'};   %Labels of the mastoids electrodes
OPTIONS_stepA.trig = {'Erg1'};                              %Label of trigger channel
OPTIONS_stepA.abr= {'Left','Right'};                        %Label of abr channels for formula 
OPTIONS_stepA.baseline = [-39, 0] ;                         %Baseline
OPTIONS_stepA.win_of_interest = [-0.04, 0.2] ;              %Epoching window
OPTIONS_stepA.eeg_elec = 1:16 ;                             %Cortical electrodes (to get cortical FFRs)
OPTIONS_stepA.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_stepA.hp = 80 ;                          % high-pass (Hz) initial value = 80
OPTIONS_stepA.lp = 3000 ;                        % low-pass (Hz) initial value = 3000
OPTIONS_stepA.bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]
OPTIONS_stepA.varhistory = 'EEG.history_stepA' ;
suffix_stepA = '_FFR_stepA';

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_stepA, count_stepA]= test_existance_of_params_in_db(OPTIONS_stepA, suffix_stepA,'') ; 

%Reref data, compute FFR formula, epoch, reject bad trials and produce
%report
if sum(flag_sub_to_create_stepA)~=0
    [preproc_filenames] = reref_filter_epoch_ffr(ALLEEG, OPTIONS_stepA,flag_sub_to_create_stepA, count_stepA,suffix_stepA) ;
end
%%/!\ some improvement to make -> add possibility to compute FFR on cortical electrodes %%

%% ------------------- Preprocess : Reject bad trials and Prepare input for BTtoolbox
OPTIONS_stepB.indir = indir ;
OPTIONS_stepB.rej_low = -45 ;                         %initial value = -45                               
OPTIONS_stepB.rej_high = 45 ;                         %initial value = 45
OPTIONS_stepB.bt_toolbox = BT_toolbox ; 
OPTIONS_stepB.varhistory = 'EEG.history_stepB' ;
OPTIONS_stepB.win_of_interest = [-0.04, 0.2] ;       %Epoching window
OPTIONS_stepB.bloc = repelem(1:30,170) ;             % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]
tube_length = 0.27  ; 
propag_sound = 340 ; 
suffix_stepB = '_stepB' ;
stepA_num = 1 ;                                      %Set of stepA parameters to use for filtering

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_stepB, count_stepB]= test_existance_of_params_in_db(OPTIONS_stepB, suffix_stepB, strcat('_stepA',num2str(stepA_num))) ; 

%Filter epoched data and prepare input for brainstem toolbox
if sum(flag_sub_to_create_stepB)~=0
    [preproc_filt_filenames] = rej_and_prepare_input_brainstem(ALLEEG, OPTIONS_stepB,tube_length, propag_sound,flag_sub_to_create_stepB, count_stepB,suffix_stepB, stepA_num);
end

%% ------------------- Display : 

% Display one participant results 
subjects_to_process = {'DVL_046_T18'} ;
%subjects_to_process = get_all_subjects(indir) ;

OPTIONS_disp.params = 'stepA1_stepB1'; 
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
OPTIONS_disp.indir = indir ; 
OPTIONS_disp.plot_dir = plot_dir ; 
OPTIONS_disp.ylim = [-0.3, 0.3] ;      % [-0.5, 0.5]
OPTIONS_disp.fs = 16384 ; 

display_individual_subjects_FFR(subjects_to_process, OPTIONS_disp) ;
%to do : display number of trials kept

%% -------------------Reject bad participants and compute group analyses

% Reject bad participants based on number of trials rejected
OPTIONS_rej.indir = indir ;
OPTIONS_rej.threshold = 1500 ;                                  % minimum number of artifact-free trials to keep a participant
OPTIONS_rej.suffix_csv = '_infos_trials_low_-45_high_45' ;     % suffix for CVS file containing trial rejection info
OPTIONS_rej.param = '_stepB1';                                 % suffix for set of parameters to process
OPTIONS_rej.visu = 1 ;                                         % 1 to display rejection rates, otherwise 0

all_subjects = get_all_subjects(indir) ;
[subjects_to_analyse] = reject_participants_FFR(all_subjects, OPTIONS_rej) ;

%% -------------------- Automatically extrcat FFR indices : neural lag



% %%TO DO !!!
% % Reject bad participants based neural lag
% OPTIONS_rej.indir = indir ;
% OPTIONS_rej.threshold = 3 ;                                    % minimum value of neural lag
% %OPTIONS_rej.suffix_csv = '_infos_trials_low_-45_high_45' ;     % suffix for CVS file containing trial rejection info
% OPTIONS_rej.param = '_stepA1';                                 % suffix for set of parameters to process
% OPTIONS_rej.visu = 0 ;                                         % 1 to display rejection rates, otherwise 0
% 
% [subjects_to_analyse] = reject_participants_FFR(subjects_to_analyse, OPTIONS_rej) ;

OPTIONS_analysis.indir = indir ;
OPTIONS_analysis.param = '_stepA1_stepB1';
grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};
OPTIONS_analysis.groups = {grpA, grpB} ;
OPTIONS_analysis.srate = 16384 ;
OPTIONS_analysis.win_of_interest = [-0.04, 0.2] ;


%OPTIONS_analysis.
%OPTIONS_analysis.

% Reject participants based on visualization
participants_to_reject = {'DVL_003_T6','DVL_007_T10','DVL_010_T24','DVL_032_T10','DVL_034_T18','DVL_011_T18','DVL_029_T10'} ;
% subjects_to_analyse = contains(subjects_to_analyse,participants_to_reject) ;
temp = find(ismember(subjects_to_analyse,participants_to_reject)) ;
subjects_to_analyse(temp) = '' ;

%Run FFR analysis only on kept subjects
FFR_analysis(subjects_to_analyse,OPTIONS_analysis);
