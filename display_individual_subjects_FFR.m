function [] = display_individual_subjects_FFR(subjects_to_process, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%OPTIONS is a structure containing:
%params = 'RFE1_REJ1';                            % option of preprocess to consider
%elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
%indir = indir ;                                  % directory path of files to process
%plot_dir = plot_dir ;                            % path to save png files of plots
%ylim = [-20,20] ;                                % limits of y axis
%fs = 16384 ;                                     % sampling rate

%Color for plot
%DIFF_color = [0 0 0]; %black
FFR_color = [0.8902 0 0]; %red

%========Time domain============

% Get timepoints
FFR_times_file = fullfile(OPTIONS.indir,'ABR_timepoints.txt') ; 
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects_to_process) %for each subject
    FFR_file = fullfile(OPTIONS.indir,subjects_to_process{loopnum},strcat(subjects_to_process{loopnum},'_',OPTIONS.params,'_abr_shifted_data_HF.txt')) ;
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

% Plot individual FFR
for jj = 1:length(subjects_to_process)
    % Plot timeseries
    figure ; 
    plot(timepoints,all_subj(:,jj),'Color', FFR_color, 'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
    grid on ; 
    %Add legend
    legend('Individual FFR', subjects_to_process(jj), 'Interpreter', 'None');
    % Adjust scales (y-axis and x-axis) (transform in milliseconds)
    xlim([min(timepoints), max(timepoints)]); ylim(OPTIONS.ylim) ; grid on ;
    %Display labels
    xlabel('Times (ms)'); ylabel('uV'); 
    %Add title
    title ([' FFR ', subjects_to_process(jj)], 'Interpreter', 'None');

print('-dsvg',fullfile(OPTIONS.indir,subjects_to_process{jj},strcat(subjects_to_process{jj},'_FFR_temporal.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat(subjects_to_process{jj},'_FFR_temporal.png')));


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

% F0_Lo = 100;
% F0_Hi = 110;
% F1_Lo = 455;
% F1_Hi = 740;
% HF_Lo = 1060;
% HF_Hi = 2750;

%%%  OPEN UP AVG FILE
%avg = openavg(FILENAME);
FS = OPTIONS.fs;
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
title(["Single-Sided Amplitude Spectrum of X(t)", subjects_to_process(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,2); 
plot(HzScale,fftFFR);
xlim([90 110]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects_to_process(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,3); 
plot(HzScale,fftFFR);
xlim([300 500]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects_to_process(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,4); 
plot(HzScale,fftFFR);
xlim([1100 1300]);
grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", subjects_to_process(jj)]);
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

print('-dsvg',fullfile(OPTIONS.indir,subjects_to_process{jj},strcat(subjects_to_process{jj},'_FFR_frequential.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat(subjects_to_process{jj},'_FFR_frequential.png')));

end

end