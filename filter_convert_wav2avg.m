

indir = '/Volumes/x9/pourASD/DEVLANG_data';

% load wav file
[y, fsamp2]= audioread(fullfile(indir,'da_170_kraus.wav')) ; 
% duraton = 170ms. 
tmp= load('/Volumes/x9/pourASD/DEVLANG_data/da_170_kraus.mat');
Fs = 16384; 

xmax = (length(y)*(1/fsamp2)*1000); 
vTimes_ms = 0:1000/fsamp2:xmax;
vTimes_ms =vTimes_ms(1:end-1); % trunc so that it is the same size as y 

EEG = tmp.EEG; 
EEG.pnts = length(y); 
EEG.srate = fsamp2;
EEG.xmax = xmax; 
EEG.times = vTimes_ms;
EEG.data = single(y');

% FILTERS the data with ERPLab
EEG  = pop_basicfilter(EEG,  1 , 'Boundary', 'boundary', 'Cutoff', [80 1500], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );

% Resample at 16k
EEG = pop_resample(EEG, Fs);

% Converts stim wav into filtered avg 
fname_out = fullfile(indir,'da_170_kraus_filtered_80_1500.txt');
fid = fopen(fname_out,'w');
fprintf(fid,'%s\n',string(EEG.data));
fclose(fid);

% Converts the output file into BT_Toolbpx compatible data 
bt_txt2avg(fname_out, Fs, 0, 170);
