function [] = compute_spectral_snr(OPTIONS,flag_sub_to_create)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly do : 

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 


% Get signal of the stimuli
sti.signal = load(fullfile(fileparts(mfilename('fullpath')),'ToolBox_BrainStem'),'da_170_kraus_16384_LP3000_HP80.wav');

vTime = readtable(fullfile(OPTIONS.indir,'ABR_timepoints.txt'));
% resampler à la même freq que le signal ABR ?
% sti.signal = resample(sti.signal,FS,sti.rate); 
   
%% Filtrer 
band_filter = [80,1500]; 
B = fir1(7,[band_filter(1),band_filter(2)]*2/FS,'bandpass'); % Filter created
sti.signalfiltered = filter(B,1,sti.signal); % Stimulus filtered
    