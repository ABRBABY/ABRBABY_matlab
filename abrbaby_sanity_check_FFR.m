function [] = abrbaby_sanity_check_FFR(eeglab_path, biosig_installer_path,indir,plots_dir) 
% FFR sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
tmp = pwd ; 
cd(eeglab_path) ; 
% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
run(biosig_installer_path) ; 
cd(tmp) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%% Variables to enter manually before running the code

%subject_of_interest = 'DVL_009_T10';

%Set variables for filtering
hp = 80; %value for high-pass filter (Hz)
lp = 3000; %value for low-pass filter (Hz)

%Rejection treshold for bad epochs
rej_low = -45;
rej_high = 45;

%Epoch window
epoch_timew =  [-0.04, 0.2] ; 
baseline_timew = [-39 0]  ; 

%Color for plot
DIFF_color = [0 0 0]; %black
FFR_color = [0.8902 0 0]; %red

% Loop through subjects
for jj=1:length(subjects) 
  
%for jj=find(ismember(subjects, subject_of_interest)) ; 

fname= dir(fullfile(indir,subjects{jj},'*.bdf'));

bdf_filename = fname.name ; 
txt_filename = strrep(bdf_filename,'.bdf','.txt');
filename = strrep(bdf_filename,'.bdf','');
filename_filter = strrep(bdf_filename,'.bdf','ABR_filtered');

%Set filepath for saving mean abr (.txt)
filepath = fullfile(indir,subjects{jj});

%% From loading file to get ABR + trigg channels

% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Modify preferences for multiple datasets loading (off)
pop_editoptions('option_storedisk', 0);

% Select bdf file in the folder
EEG = pop_biosig(fullfile(indir,subjects{jj},bdf_filename));

% Rename file
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',filename,'gui','off');
EEG = eeg_checkset( EEG );

% Select channels to keep
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'channel',{'Erg1','Left','Right'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 

eeglab redraw;
%%
% Compute ABR channel from Left / Right channels
% forumla : {(Left)+(Right)}/-2 = Ref - {(LA+RA)/2}
EEG.nbchan = EEG.nbchan+1;
if ~isempty(EEG.chanlocs)
	EEG.chanlocs(end+1).labels = 'ABR';
    EEG.chanlocs(end).type = 'EEG';
end

for m = 1:size(EEG.data,2)
    EEG.data(4,m) = ((EEG.data(1,m))+(EEG.data(2,m)))/2; 
end
%%
% Keep only ABR channel and triggers (select channel)
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'channel',{'ABR', 'Erg1'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 

% Extract event from trigger channel (Erg1)
EEG = pop_chanevent(EEG, 1,'oper','X>20000','edge','leading','edgelen',1);

eeglab redraw;

%% Reject bad events 
    
% Remove events outliers (e.g. boundaries) or too close events 
EEG.event(find(diff([EEG.event.latency])<0.1*EEG.srate)+1) = [] ; % minimum intretrial duration = 220 ms
EEG.event(find(diff([EEG.event.latency])>2*EEG.srate))= [] ; % maximum intertrial duration = 1 second
%EEG = pop_editeventvals(EEG,'delete',1); %when first event must be removed (6001 events)

%% Replace event type by infos from .txt
% Read .txt
my_events = readtable(fullfile(indir,subjects{jj}, txt_filename), 'ReadVariableNames', 0);

%Insert info from .txt into EEG.event
my_events = table2array(my_events);

temp = struct('latency', {EEG.event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {EEG.event(:).urevent});

EEG.event = temp;          
eeglab redraw;

%% Extract epochs for HF

EEG = pop_epoch( EEG, {  'HF'  }, epoch_timew, 'newname', 'epochs', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off') 

%Remove baseline
EEG = pop_rmbase( EEG, baseline_timew,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 

%Extract condition and latencies of all HF trials
condition = {EEG.event.type} ;
latency = [EEG.event.latency]/EEG.srate;
trial_index = [EEG.event.urevent];

%Epoch rejection
EEG = eeg_checkset( EEG );

binf = EEG.times(1)/1000 ; 
bsup = EEG.times(end)/1000;

EEG = pop_eegthresh(EEG,1,1,rej_low,rej_high,binf,bsup,0,1);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'setname','epochs_rej','gui','off'); 

eeglab redraw;

% Export information about rejected trials
trial_rej = [EEG.reject.rejthresh];
bloc = zeros(5100,1);
tr = 1;
for bc = 1:30
    bloc(tr:tr+169,1)= repmat(bc,170,1);
    bc = bc+1;
    tr = tr+170;
end

% Create table to store these information
list_trial_infos = table(trial_index',condition',latency',trial_rej', bloc) ;
%  Save this table into a csv file (use function writetable)
writetable(list_trial_infos,fullfile(indir,subjects{jj},strcat(filename,'_low_',num2str(rej_low),'_high_',num2str(rej_high),'infos_trials_FFR.csv'))) ; 
%writetable(list_trial_infos,fullfile(indir,subject_of_interest,strcat(filename,'_low_',num2str(rej_low),'_high_',num2str(rej_high),'infos_trials_FFR.csv'))) ; 


%% Filtering

% Filter data
EEG  = pop_basicfilter( EEG,  1 , 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 ); % GUI: 11-Apr-2022 12:47:48
% EEG  = pop_basicfilter( EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 6, 'savenew',fullfile(indir,subjects{jj},filename_filter),'gui','off');
 
%Extract mean activity (erp) and replace data
abr = mean(EEG.data(1,:,:),3);
EEG.data = abr;

eeglab redraw;

% Add tube delay (27 cm x 340 m/s ) 
nsample_delay = fix(EEG.srate * (0.27 / 340) ) ; 

abr_shifted = circshift(abr,nsample_delay) ;


%% Export ABR data into .txt file
fname_out = fullfile(filepath,strrep(bdf_filename,'.bdf','_abr_shifted_data_HF.txt')) ;
fid = fopen(fname_out,'w');
fprintf(fid,'%c\n',abr_shifted);
fclose(fid);

addpath 'C:\Users\hervé\Documents\GitHub\ABRBABY\ToolBox_BrainStem\BT_2013\programFiles';
bt_txt2avg(fname_out, EEG.srate, epoch_timew(1)*1000, epoch_timew(2)*1000);

end

%% Export times
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
beep;


%% VISUALIZATION

%========Time domain============

% Get timepoints
FFR_times_file = fullfile(indir,'ABR_timepoints.txt') ; 
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects) %for each subject
%for loopnum=find(ismember(subjects,'subject_of_interest')) ; 
    FFR_file = fullfile(indir,subjects{loopnum},strcat(subjects{loopnum},'_abr_shifted_data_HF.txt')) ; 
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

% Plot individual FFR
for jj = 1:size(subjects,1)
%for jj=1 
    figure ; 
    plot(timepoints,all_subj(:,jj),'Color', FFR_color, 'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
    grid on ; 
    legend('Individual FFR', subjects(jj), 'Interpreter', 'None');
    %legend(['Individual FFR ', subject_of_interest], 'Interpreter', 'None');
    xlabel('Times (ms)'); ylabel('uV'); title ([' FFR ', subjects(jj)], 'Interpreter', 'None');
    %xlabel('Times (ms)'); ylabel('uV'); title ([' FFR ', subject_of_interest], 'Interpreter', 'None');

%end

print('-dsvg',fullfile(indir,subjects{jj},strcat(subjects{jj},'_FFR_temporal.svg')));
print('-dpng',fullfile(plots_dir,strcat(subjects{jj},'_FFR_temporal.png')));

%%
%========Frequency domain============

% Adaptation of Skoe function (bt_fftsc)
 %function [Freq1 Freq2 Freq3 fftFFR HzScale]=bt_fftsc(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% bt_fftsc computes frequency-domain amplitudes of three user-defined 
% frequency bins of Brainstem response.  Results are scaled to peak µV.
%
% Usage: [F0 F1 HF] = bt_fftsc('filename.avg',10,40,100,150,300,350,600,800);
%    over the range of 10 to 40 ms, finds average frequency amplitude of
%    100-150 Hz, 300-350 Hz and 600-800 Hz.
% 
% Three variables are returned the workspace:
% (1) Freq1: mean amplitude over F0_Lo-F0_Hi Hz
% (2) Freq2: mean amplitude over F1_Lo-F1_Hi Hz
% (3) Freq3: mean amplitude over HF_Lo to HF_Hi Hz 

% % Dependancies: ms2row, openavg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F0_Lo = 100;
F0_Hi = 110;
F1_Lo = 455;
F1_Hi = 740;
HF_Lo = 1060;
HF_Hi = 2750;

%%%  OPEN UP AVG FILE
%avg = openavg(FILENAME);
FS = 16384;
%FFR = all_subj(:,1);
FFR = all_subj(:,jj);
%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD 

%startPoint = ms2row(avg, start);
%endPoint = ms2row(avg, stop);

%FFR = data(startPoint:endPoint);
numPoints = length(FFR);

%**** STEP 2. FFT OF FFR
%******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
rampMS = 4/1000; % length of ramp (on and off) in seconds
% hanPoints = 26;  %hard coded to be the same as Biologic's settings (December 7, 2005);
hanPoints = rampMS.*FS; % length of ramp in points
hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
hanHalfPoints = round(hanPoints./2);
numberOfOnes = numPoints - hanPoints;
FFRhan = hann(hanPoints);  
FFRhan = [FFRhan(1:hanHalfPoints); ones(numberOfOnes,1); FFRhan(hanHalfPoints+1:hanPoints)];

% baseline, window, then baseline again
FFR = detrend(FFR, 'constant');
FFR = detrend(FFR.*FFRhan, 'constant');

%******** STEP 2b. Perform FFT
fftFFR = abs(fft(FFR, round(FS)));
fftFFR = fftFFR(1:round(round(FS)/2));
fftFFR = fftFFR.*(2./numPoints); % scale to peak µV
HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
HzScale = HzScale(1:length(fftFFR));

%% Plot FFT in specific frequency windows
figure ;

subplot(2,2,1);
plot(HzScale,fftFFR);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,2); 
plot(HzScale,fftFFR);
xlim([90 110]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,3); 
plot(HzScale,fftFFR);
xlim([300 500]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,4); 
plot(HzScale,fftFFR);
xlim([1100 1300]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

print('-dsvg',fullfile(indir,subjects{jj},strcat(subjects{jj},'_FFR_frequential.svg')));
print('-dpng',fullfile(plots_dir,strcat(subjects{jj},'_FFR_frequential.png')));

end

end