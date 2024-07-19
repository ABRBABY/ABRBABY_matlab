% FFR Pipeline analysis for ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Reads initial database directory 
indir = fileread('indir.txt');
[data_pathname,rawdata_fname,~] = fileparts(indir);

plot_dir = fullfile(fileparts(indir),strcat('plot_dir_',rawdata_fname));

% Get the list of subjects in the databse
list_subjects = get_subjects(indir,fullfile(indir, '')) ;

% This function sets custom path 
[eeglab_path, biosig_installer_path, erplab_path, BT_toolbox] = get_custom_path();

OPTIONS.file = fullfile(indir,'force_rerun_participants.csv') ;

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path, BT_toolbox) ;

%% ------------------- Preprocess : reref, set chan positions, filter
OPTIONS_stepA.indir = indir;
OPTIONS_stepA.mastos = {'Lmon','Rmon','MASTOG','MASTOD'};   %Labels of the mastoids electrodes
OPTIONS_stepA.trig = {'Erg1'};                              %Label of trigger channel
OPTIONS_stepA.abr= {'Left','Right'};                        %Label of abr channels for formula 
OPTIONS_stepA.conditions = {'HF'} ; 
OPTIONS_stepA.baseline = [-39, 0] ;                         %Baseline
OPTIONS_stepA.win_of_interest = [-0.04, 0.2] ;              %Epoching window
OPTIONS_stepA.eeg_elec = 19 ;                             %Cortical electrodes (to get cortical FFRs)
OPTIONS_stepA.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
OPTIONS_stepA.hp = 80 ;                          % high-pass (Hz) initial value = 80
OPTIONS_stepA.lp = 3000 ;                        % low-pass (Hz) initial value = 3000
OPTIONS_stepA.bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]
OPTIONS_stepA.varhistory = 'EEG.history_stepA' ;
OPTIONS_stepA.analysis = 'FFR';
OPTIONS.file = fullfile(indir,'force_rerun_participants.csv') ;

% Update for Clem : now just use this suffix for output (overwrite if
% exist, otherwise creates new one)
output_suffix = '_stepA2'; % change this name to create new dataset with different parameters WARNING : Estelle has generated A1 (do not overwrite?)

%Subjects to process : when whant to choose
if exist(OPTIONS.file,'file') && ~isempty(fileread(OPTIONS.file))
   subj_to_process  = get_subjects(indir,OPTIONS);
else 
   subj_to_process  = get_subjects(indir,[]);
end

flag_sub_to_create_stepA = (contains(list_subjects,subj_to_process))';

%Reref data, compute FFR formula, epoch, reject bad trials and produce report
reref_filter_epoch(ALLEEG, OPTIONS_stepA, flag_sub_to_create_stepA, str2num(output_suffix(end)),output_suffix(1:end-1)) ;

% Prints out message on progress
fprintf('JUST FINISHED STEP A\n');

%%/!\ some improvement to make -> add possibility to compute FFR on cortical electrodes %%

%% ------------------- Preprocess : Reject bad trials and Prepare input for BTtoolbox
OPTIONS_stepB.indir = indir ;
OPTIONS_stepB.analysis = 'FFR' ;
OPTIONS_stepB.rej_low = -25 ;                         %initial value = -45                               
OPTIONS_stepB.rej_high = 25 ;                         %initial value = 45
OPTIONS_stepB.bt_toolbox = BT_toolbox ; 
OPTIONS_stepB.varhistory = 'EEG.history_stepB' ;
OPTIONS_stepB.win_of_interest = [-0.04, 0.2] ;       %Epoching window
OPTIONS_stepB.eeg_elec = 'ABR';
 
% Update for Clem : now just use this suffix for output (overwrite if
% exist, otherwise creates new one)
output_suffix = '_stepA1_stepB2'; % change this name to create new dataset with different parameters WARNING : Estelle has generated A1_B1 (do not overwrite?)

%Subjects to process : when whant to choose
if exist(OPTIONS.file,'file') && ~isempty(fileread(OPTIONS.file))
   subj_to_process  = get_subjects(indir,OPTIONS);
else 
   subj_to_process  = get_subjects(indir,[]);
end

flag_sub_to_create_stepB = (contains(list_subjects,subj_to_process))';

%Filter epoched data and prepare input for brainstem toolbox
reject_bad_trials(ALLEEG, OPTIONS_stepB, 'unbalanced', flag_sub_to_create_stepB, str2num(output_suffix(end)), '_stepB',strcat('_',strtok(output_suffix,'_'))) ; 

% Prints out message on progress
fprintf('JUST FINISHED STEP B\n');

%% -------------------  Prepare output for BT_Toolbox + optionnal display
OPTIONS_abr.indir = indir ; 
OPTIONS_abr.display = 1 ; 
OPTIONS_abr.savefigs = 1 ; 
OPTIONS_abr.abr_disp_scale = [-0.2, 0.2];         % scale for display ABR ave
OPTIONS_abr.plot_dir = plot_dir ; 
OPTIONS_abr.png_folder = fullfile(plot_dir,'png_folder');                          % path to save png files of plots
OPTIONS_abr.svg_folder =  fullfile(plot_dir,'svg_folder');
OPTIONS_abr.fig_folder = fullfile(plot_dir,'fig_folder');
 
tube_length = 0.27 ;  % meter
propag_sound =  340 ; % vitesse propagation son meter / sec

if OPTIONS_abr.savefigs ==1 ; create_plot_dirs_if_does_not_exist(plot_dir); end 

%Subjects to process : when whant to choose
if exist(OPTIONS.file,'file') && ~isempty(fileread(OPTIONS.file))
   subj_to_process  = get_subjects(indir,OPTIONS);
else 
   subj_to_process  = get_subjects(indir,[]);
end

flag_sub_to_create_abr = (contains(list_subjects,subj_to_process))';

% The following line should only prepare input for brainstem 
prepare_input_brainstem(ALLEEG, OPTIONS_abr,tube_length, propag_sound,flag_sub_to_create_abr, str2num(output_suffix(end)),'_stepB', strcat('_',strtok(output_suffix,'_')));

% Prints out message on progress
fprintf('JUST FINISHED PREPARE INPUT BRAINSTEM\n');


%% -------------------Compute neural lag for all subject and write a table
OPTIONS_neural.params = 'stepA1_stepB2'; 
OPTIONS_neural.ffr_polarity = 'avg' ;                %polarity of the ffr ('avg', 'pos' or 'neg')
OPTIONS_neural.indir= indir;
OPTIONS_neural.plot_dir = plot_dir ;
OPTIONS_neural.stim = 'da_170_kraus_16384_LP3000_HP80.avg' ;
OPTIONS_neural.start = 0 ;      % first time point in the simulti to compute neural lag 
OPTIONS_neural.stop = 169 ;     % last time point in the simulti to compute neural lag 
OPTIONS_neural.lagstart = 3 ;   % first time point to search for neural lag
OPTIONS_neural.lagstop = 10 ;   % last time point to search for neural lag
OPTIONS_neural.polarity = 'ABSOLUTE' ;                %sign of max correlation value ('POSITIVE', 'NEGATIVE', or 'ABSOLUTE')
OPTIONS_neural.display = 0 ;         % 1 to display r pearson correlation distribution
OPTIONS_neural.table = 0 ;         % 1 to save table with all neural lags 
% OPTIONS_neural.BT_toolbox = BT_toolbox ;
% OPTIONS_neural.grp = [{'_T3','_T6','_T8','_T10'},{'_T18','_T24'}];

%Subjects to process : when whant to choose
flag_sub_to_compute_nlag = ~test_existance_of_BT_toolbox(OPTIONS_neural) ; 

if exist(OPTIONS.file,'file') && isempty(fileread(OPTIONS.file))
   subj_to_process  = get_subjects(indir,OPTIONS);
   flag_sub_to_compute_nlag  = (contains(list_subjects,subj_to_process))';
end

% Computes the neural lag
neural_lag = compute_neural_lag(OPTIONS_neural,flag_sub_to_compute_nlag ) ;

% Prints out message on progress
fprintf('JUST FINISHED COMPUTE NEURAL LAG\n');

%% -------------------Compute SNRs 
OPTIONS_SNR.params = 'stepA1_stepB1';                       % parameters to run
% OPTIONS_SNR.elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};  % electrode subset for cortical FFR
OPTIONS_SNR.indir = indir ;                                 % indir
OPTIONS_SNR.plot_dir = plot_dir ;                           % path to plots folder
OPTIONS_SNR.ylim = [-0.5, 0.5] ;                            % [-0.5, 0.5] 
OPTIONS_SNR.ffr_polarity = 'avg' ;                          % polarity of the FFR: ('avg', 'pos' or 'neg')
OPTIONS_SNR.winNoise = cat(2,80:1:95,105:1:120);            % windows for noise and signal
OPTIONS_SNR.winSignal = 95:1:105;                           % timewindow +/- 5 around F0
OPTIONS_SNR.win_of_interest = [-0.04, 0.2] ;                % epoch limits
OPTIONS_SNR.timew_F0 = [55 200] ;                           % timewindow of FFR on which to compute F0 (in ms)
OPTIONS_SNR.display = 0 ;                                   % 1 if want to display SNR plots
OPTIONS_SNR.savefig = 0 ;                                   % 1 if want to save figures

%Subjects to process : when whant to choose
flag_sub_to_create_ffr = ~test_existance_of_BT_toolbox(OPTIONS_SNR) ; 

if exist(OPTIONS.file,'file') && isempty(fileread(OPTIONS.file))
   subj_to_process  = get_subjects(indir,OPTIONS);
   flag_sub_to_create_ffr  = (contains(list_subjects,subj_to_process))';
end

% tran = [10 55];  % Time window of the stimulus consonant transition(ms)  
% cons = [55 170]; % Time windows of the stimulus constant portion(ms)
% baseline = [-40 0];  % Prestimulus time windows (ms)

%Filter epoched data and prepare input for brainstem toolbox
[spectral_snr,aWin,freq_harmonics,max_psd] = compute_spectral_snr(OPTIONS_SNR, flag_sub_to_create_ffr, neural_lag) ; 
fprintf('JUST FINISHED COMPUTE SNR\n');


%% ------------------- Display SNR violin
% With the piece of code you can modulate how many viollin plot to display,
% you can combien differently the groups, all depends on the content of the
% filed .groups (and colors must be the colors accordingly tothe groups) 

% OPTIONS_display_violin.groups = {{'_T6','_T8'},{'_T18','_T24'},{'_T10'}};
% OPTIONS_display_violin.groups = {{'_T6'},{'_T8'},{'_T18'},{'_T24'},{'_T10'}};
OPTIONS_display_violin.groups = {{'_T6','_T8'},{'_T10'},{'_T18'},{'_T24'}};
OPTIONS_display_violin.colors = {[0,0,1],[0,1,0],[0.5,0.5,0],[0.5,0,0.5]}; 
OPTIONS_display_violin.indir = indir ; 
% OPTIONS_display_violin.title = {'Contrats between SNR : f0, win transition '} ; 
% OPTIONS_display_violin.title = {'Contrats between SNR : f0, vowel'} ; 
OPTIONS_display_violin.title = {'Contrats between SNR : f0, baseline'} ; 

% Exclude some participant based on max_psd 
% size(flag_sub_to_disp) --> nb subjects in the whole database
% sum(flag_sub_to_disp) --> nb subjects which will be processed from this point
flag_sub_to_disp = (max_psd>100.3-4)&(max_psd<100.3+4) ; 

% Here spectral_snr dimension is nSubj x time window (3) x spectral bin (f0+n) 
% Dim 1 : nSubject
% Dim 2 : window (1=transition, 2=vowel , 3=baseline)
% Dim 3 : harmonics (1=f0, 2= next, etc.) 
% Ex : spectral_snr(:,2,1) : all subjects vowel f0
% Ex : spectral_snr(:,1,2) : all subjects transition, 2nd harmonic
spectral_snr_to_proc = spectral_snr(:,1,1); 

% BELOW some different displays (comment/uncomment the one you prefer) 
% plot_violin_variable_nb_cond(OPTIONS_display_violin, flag_sub_to_create_ffr, spectral_snr(:,3,1));
plot_variable_nb_cond(OPTIONS_display_violin, flag_sub_to_disp, spectral_snr_to_proc, neural_lag);
% plot_hist_nb_cond(OPTIONS_display_violin, flag_sub_to_create_ffr, spectral_snr(:,1,1));
% plot_subplot_nb_cond(OPTIONS_display_violin, flag_sub_to_create_ffr, spectral_snr(:,1,1), neural_lag);


%% ------------------- Compute Pitch tracking 
OPTIONS_pitch.indir = indir; 
OPTIONS_pitch.blocksz = 40 ; 
OPTIONS_pitch.step= 1 ; 
OPTIONS_pitch.startSTIM = 55 ; 
OPTIONS_pitch.endSTIM = 169 ; 
OPTIONS_pitch.expectedNeuralag= 10 ; 
OPTIONS_pitch.stim = 'da_170_kraus_16384_LP3000_HP80.avg' ;
OPTIONS_pitch.BT_toolbox = BT_toolbox; 

OPTIONS_pitch.minFrequencyR = 100; 
OPTIONS_pitch.maxFrequencyR = 120; 
OPTIONS_pitch.minFrequency_stim = 80; 
OPTIONS_pitch.maxFrequency_stim = 120; 

% [PITCH_ERROR_AC,PITCH_ERROR_FFT,  PITCH_STRENGTH2, PITCH_SRCORR, vTime, vFreqAC, vFreqFFT, vTime_stim, vFreqAC_stim, vFreqFFT_stim] = compute_pitchtracking(OPTIONS_pitch, flag_sub_to_create_ffr); 
[PITCH_ERROR_AC,PITCH_ERROR_FFT,  PITCH_STRENGTH2, PITCH_SRCORR, vTime, vFreqAC, vFreqFFT, vTime_stim, vFreqAC_stim, vFreqFFT_stim] = compute_pitchtracking(OPTIONS_pitch, flag_sub_to_create_ffr); 

% 
% %% ------------------- Display Pitch violin
% OPTIONS_display_violin.groups = {{'_T8'},{'_T24'},{'_T10'}};
% OPTIONS_display_violin.colors = {[1,0,0],[0,0,1],[0,1,0]}; 
% OPTIONS_display_violin.indir = indir ; 
% OPTIONS_display_violin.title = {'Pitch errors'} ; 
% 

% %% One subject display 
% ss=1 ; figure ; subplot(2,1,1) ; plot(vTime(ss,:),vFreqFFT(ss,:),'s', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6) ; hold on ; plot(vTime_stim(ss,:),vFreqFFT_stim, 'k', 'LineWidth', 2); subplot(2,1,2) ; plot(vTime(ss,:),vFreqAC(ss,:),'s', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6) ; hold on ; plot(vTime_stim(ss,:),vFreqAC_stim, 'k', 'LineWidth', 2);
% 
% %% Mean group display 
% figure ; plot(vTime(1,:),mean(vFreqAC,1),'s', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6) ; hold on ; plot(vTime_stim(1,:),mean(vFreqAC_stim,1), 'k', 'LineWidth', 2);

%% TODOS next : subplot by group (same nb of subplot than groups) 
% OPTIONS_display_violin.groups = {{'_T6','_T8'},{'_T18','_T24'},{'_T10'}};
% OPTIONS_display_violin.groups = {{'_T6'},{'_T8'},{'_T18'},{'_T24'},{'_T10'}};
OPTIONS_display_violin.groups = {{'_T8','_T24','_T10','T24'}};
OPTIONS_display_violin.groups = {{'_T6','_T8'},{'_T10'},{'_T18'},{'_T24'}};

% OPTIONS_display_violin.groups = {{'_T8'},{'_T24'},{'_T10'}};
OPTIONS_display_violin.colors = {[1,0,0],[0,0,1],[0,1,0], [0.5,1,0]}; 
OPTIONS_display_violin.indir = indir ; 
OPTIONS_display_violin.title = {'ADAPT THIS TITLE'} ; 
% plot_pitchtrack(OPTIONS_display_violin, flag_sub_to_create_ffr, vTime(1,:), vFreqAC, vTime_stim(1,:),vFreqAC_stim);
% plot_pitchtrack(OPTIONS_display_violin, flag_sub_to_create_ffr, vTime(1,:), vFreqFFT, vTime_stim(1,:),vFreqFFT_stim);
% plot_pitchtrack(OPTIONS_display_violin, flag_sub_to_create_ffr, vTime(1,:), vFreqAC, vTime_stim(1,:),vFreqAC_stim);
plot_pitchtrack(OPTIONS_display_violin, flag_sub_to_disp, vTime(1,:), vFreqAC, vTime_stim(1,:),vFreqAC_stim);


%% ------------------- Display Pitch violin
OPTIONS_display_violin.groups = {{'_T8'},{'_T24'},{'_T10'}};
OPTIONS_display_violin.colors = {[1,0,0],[0,0,1],[0,1,0]}; 
OPTIONS_display_violin.indir = indir ; 
OPTIONS_display_violin.title = {'Pitch errors'} ; 
plot_violin_variable_nb_cond(OPTIONS_display_violin, flag_sub_to_create_ffr, PITCH_ERROR_AC');
plot_violin_variable_nb_cond(OPTIONS_display_violin, flag_sub_to_create_ffr, PITCH_SRCORR');

% % Or choose subjects with csv file
% subjects_to_process = get_subjects(indir, []) ;
% 
% snr = FFR_analysis_get_SNR_freq(subjects_to_process, OPTIONS_SNR) ; 
% 
% % Write a table with SNR info
% SNR_all = table(subjects_to_process, snr', 'VariableNames', {'suject_ID', 'SNR'}) ;
% writetable(SNR_all,fullfile(OPTIONS_SNR.indir,strcat('all_SNRs_F0s_', OPTIONS_SNR.ffr_polarity, '_ffr_', num2str(OPTIONS_SNR.timew_F0(1)), '_', num2str(OPTIONS_SNR.timew_F0(2)), 'tw_', OPTIONS_SNR.params,'.csv')), 'WriteVariableNames', true) ;
