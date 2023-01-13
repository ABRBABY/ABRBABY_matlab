% % Potential things to improve : 

% Create an history_abr variable to append to the EEG structure to document
% the operations and implement the overwrite option

%% ------------------- Set environment 
% Variables to enter manually before running the code
eeglab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/eeglab2021.1' ; 
erplab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/erplab8.30';
biosig_installer_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/biosig4octmat-3.8.0/biosig_installer.m' ; 

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : reref, epoch, set chan positions, reject BAD trials 
indir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/' ;
mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; trig = {'Erg1'}; abr= {'Left','Right'};  % Ref and trigger channels 
baseline = [-39, 0] ; win_of_interest = [-0.04, 0.2] ; 
eeg_elec = 1:16 ; 
chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
overwrite = 1 ; % this option allow to overwrite (=11) or not (=0) 

rej_low = -45; %150 infants; 120 adults
rej_high = 45; %150 infants; 120 adults
% abr_elec = 17 ; 
bloc = repelem(1:30,170) ; % creates a vector of [1 1 1 1 (170 times) 2 2 2 2 (170 times) etc. up to 30]

[preproc_filenames] = reref_epoch_ffr(ALLEEG, indir, mastos, trig, abr, eeg_elec, baseline, win_of_interest, chan_dir,rej_low,rej_high, bloc, overwrite);

%% ------------------- Preprocess : Filter 

hp =80; % high-pass (Hz) (APICE)
lp = 3000; % low-pass (Hz) (APICE) 

overwrite = 0 ; % this option allow to overwrite (=1) or not (=0) 

tube_length = 0.27  ; 
propag_sound = 340 ; 

[preproc_filt_filenames] = filter_and_prepare_input_brainstem(out, hp,lp, overwrite) ; 


%% This function should includes the followoin g operatyions : 

% Filtering
for jj=1:length(out) 
    
    % Filter data
    EEG  = pop_basicfilter( out{jj},  1 , 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 ); % GUI: 11-Apr-2022 12:47:48
    
    %Extract mean activity (erp) and replace data
    abr = mean(EEG.data(1,:,:),3);
    EEG.data = abr;

    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    abr_shifted = circshift(abr,nsample_delay) ;

    %% Export ABR data into .txt file
    fname_out = fullfile(fullfile(indir,subjects{jj}),strrep(bdf_filename,'.bdf','_abr_shifted_data_HF.txt')) ;
    fid = fopen(fname_out,'w');
    fprintf(fid,'%c\n',abr_shifted);
    fclose(fid);

    addpath 'C:\Users\hervé\Documents\GitHub\ABRBABY\ToolBox_BrainStem\BT_2013\programFiles';
    bt_txt2avg(fname_out, EEG.srate, epoch_timew(1)*1000, epoch_timew(2)*1000);

end

%% Export times (from any sibject : just timepoints) 
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);

% 

%% ------------------- Display : 
elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
plot_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder' ; 
display_timeseries_by_condition(preproc_filenames, elec_subset, 'balanced',plot_dir) ; 
