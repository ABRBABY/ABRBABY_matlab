% FFR Pipeline analysis for ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
% custom_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data') ;

plot_dir = fullfile(custom_path, 'plot_dir');

% list_subjects = get_subjects(indir,fullfile(indir, 'participants_included.xlsx')) ;
list_subjects = get_subjects(indir,fullfile(indir, '')) ;

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

% % Choose subject to process
% subj_to_process = {'DVL_021_T18','DVL_040_T24','DVL_043_T8','DVL_048_T10','DVL_051_T8','DVL_052_T6'} ;
% flag_sub_to_create_stepA = (contains(list_subjects,subj_to_process))';

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
tube_length = 0.27 ; 
propag_sound = 340 ; 
suffix_stepB = '_stepB' ;
stepA_num = 1 ;                                      %Set of stepA parameters to use for filtering

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_stepB, count_stepB]= test_existance_of_params_in_db(OPTIONS_stepB, suffix_stepB, strcat('_stepA',num2str(stepA_num))) ; 

% % Choose subject to process
% subj_to_process = {'DVL_003_T18','DVL_012_T24','DVL_018_T18','DVL_021_T18','DVL_040_T24','DVL_043_T8','DVL_048_T10','DVL_051_T8','DVL_052_T6'} ;
% flag_sub_to_create_stepB = (contains(list_subjects,subj_to_process))';

%Filter epoched data and prepare input for brainstem toolbox
if sum(flag_sub_to_create_stepB)~=0
    [preproc_filt_filenames] = rej_and_prepare_input_brainstem(ALLEEG, OPTIONS_stepB,tube_length, propag_sound,flag_sub_to_create_stepB, count_stepB,suffix_stepB, stepA_num);
end

%% ------------------- Display : 
OPTIONS_disp.params = 'stepA1_stepB1';
OPTIONS_disp.polarity = 'avg' ;                             % polarity of the FFR: ('avg', 'pos' or 'neg')
OPTIONS_disp.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
OPTIONS_disp.indir = indir ; 
OPTIONS_disp.plot_dir = plot_dir ; 
OPTIONS_disp.ylim = [-0.5, 0.5] ;      % [-0.5, 0.5]
OPTIONS_disp.fs = 16384 ; 
OPTIONS_disp.file = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Analyses\ffr_participants_todecide.csv';

% Display one participant results 
subjects_to_process = {'DVL_056_T10'} ;

% Or display all subjects
% subjects_to_process = list_subjects ;

% Or choose subjects with csv file
% subjects_to_process = get_subjects(indir, OPTIONS_disp) ;

display_individual_subjects_FFR(subjects_to_process, OPTIONS_disp) ;
%to do : display number of trials kept


%% -------------------Compute neural lag for all subject and write a table
OPTIONS_neural.params = 'stepA1_stepB1'; 
OPTIONS_neural.ffr_polarity = 'avg' ;                %polarity of the ffr ('avg', 'pos' or 'neg')
OPTIONS_neural.indir= indir; 
OPTIONS_neural.stim = 'da_170_kraus_16384_LP3000_HP80.avg' ;
OPTIONS_neural.start = 0 ;
OPTIONS_neural.stop = 169 ;
OPTIONS_neural.lagstart = 3 ;
OPTIONS_neural.lagstop = 9 ;
OPTIONS_neural.polarity = 'POSITIVE' ;                %sign of max correlation value ('POSITIVE', 'NEGATIVE', or 'ABSOLUTE')
OPTIONS_neural.chan =1 ;
OPTIONS_neural.chancomp =1 ;
OPTIONS_neural.BT_toolbox = BT_toolbox ;
OPTIONS_neural.grpA = {'_T3','_T6','_T8','_T10'};
OPTIONS_neural.grpB = {'_T18','_T24'};

subjects_to_process = get_subjects(indir, '') ;
% subjects_to_process = list_subjects ;
% subjects_to_process = get_subj('\\Filer\home\Invites\herve\Mes documents\These\EEG\Analyses\participants_rejA.xlsx') ;
% subjects_to_process = {'DVL_TestGiu'};

compute_neural_lag_report(subjects_to_process, OPTIONS_neural) ; 

%% -------------------Export infos on rejection
OPTIONS_rejinfo.indir = indir ;
OPTIONS_rejinfo.param = '_stepA1_stepB1';
subjects_rejinfo = get_subjects(indir, '') ;

export_rejection_infos_FFR(subjects_rejinfo,OPTIONS_rejinfo) ;

%% -------------------Reject bad participants and compute group analyses
% Set options to reject bad participants based on number of trials rejected and neural lag
OPTIONS_rej.indir = indir ;
OPTIONS_rej.threshold = 3000 ;                                 % minimum number of artifact-free trials to keep a participant
OPTIONS_rej.suffix_csv = '_infos_trials_low_-45_high_45' ;     % suffix for CVS file containing trial rejection info
OPTIONS_rej.param = '_stepB1';                                 % suffix for set of parameters to process
OPTIONS_rej.visu = 1 ;                                         % 1 to display rejection rates, otherwise 0
OPTIONS_rej.neural_lag = 3 ;                                   % neural lag threshold : under this value, subjects are rejected
OPTIONS_rej.ffr_polarity = 'avg' ;                             % which FFR polarity to use
OPTIONS_rej.polarity = 'POSITIVE' ;                            % which correlation value for neural lag to use
% OPTIONS_rej.file = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Analyses\ffr_participants_todecide.csv';
OPTIONS_rej.file = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data\ffr_participants_ok.csv';

% Choose subjects to analyse
% subjects_to_analyse = get_subjects(indir, '') ;                      % get all subjects
subjects_to_analyse = get_subjects(indir, OPTIONS_rej) ;           % get subjects in OPTIONS_rej.file

% % Reject bad participants based on number of trials rejected and neural lag
% all_subjects = get_subjects(indir, '') ;
% [subjects_to_analyse] = reject_participants_FFR(all_subjects, OPTIONS_rej) ;
% ----- OR -------
% % Reject bad participants based on visualization
% participants_to_reject = {'DVL_008_T10','DVL_010_T24', 'DVL_021_T18','DVL_032_T10','DVL_034_T18'} ;
% subjects_to_analyse(contains(subjects_to_analyse,participants_to_reject)) = [] ;

% Set options to compute group analyses
OPTIONS_analysis.indir = indir ;
OPTIONS_analysis.param = '_stepA1_stepB1';
grpA.suffix = {'_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};
OPTIONS_analysis.groups = {grpA, grpB} ;
OPTIONS_analysis.srate = 16384 ;
OPTIONS_analysis.win_of_interest = [-0.04, 0.2] ;
OPTIONS_analysis.neural_lag = OPTIONS_rej.neural_lag ;
OPTIONS_analysis.ffr_polarity = OPTIONS_rej.ffr_polarity ; 
OPTIONS_analysis.polarity = OPTIONS_rej.polarity  ;
OPTIONS_analysis.plot_dir = plot_dir ; 
OPTIONS_analysis.plot_FFT = 0 ; 
OPTIONS_analysis.stim_avg = 'C:\Users\herve\Documents\GitHub\ABRBABY_matlab\ToolBox_BrainStem\BT_2013\da_170_kraus_16384_LP3000_HP80.avg' ;
OPTIONS_analysis.woi_F0 = [90:110];   %add option to compute F0 in woi (in Hz) or no -> []
OPTIONS_analysis.nlag_filename = 'all_neural_lags_avg_ffr_POSITIVE_corr_3_9.csv' ;

% Run FFR analysis only on kept subjects
% FFR_analysis(subjects_to_analyse,OPTIONS_analysis);
FFR_analysis_freq(subjects_to_analyse,OPTIONS_analysis);

