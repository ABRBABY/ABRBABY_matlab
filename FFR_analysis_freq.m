function [] = FFR_analysis_freq(subjects, OPTIONS)
% ERPs analysis script -
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

grpA = OPTIONS.groups{1} ;
grpB = OPTIONS.groups{2} ;

%Color for plot
FFR_color = [0 0 0]; %black
grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.8902 0 0]; %red

% Get timepoints
FFR_times_file = fullfile(OPTIONS.indir,'ABR_timepoints.txt') ;
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects) %for each subject
    FFR_file = fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.txt')) ;
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR preprocessing steps on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

%% Compute FFR

% Compute grand average
grd_FFR = mean(all_subj,2);

% Compute grand average by group
grpA.subj = subjects(contains(subjects,grpA.suffix));
grpB.subj = subjects(contains(subjects,grpB.suffix));

grpA.data = zeros(size(timepoints,1),1);
grpB.data = zeros(size(timepoints,1),1);

groups = {grpA.subj, grpB.subj};
data_groups = {grpA.data, grpB.data};
for k = 1:length(groups)
    grp_subj = zeros(size(timepoints,1),1);
    indices = find(ismember(subjects,groups{k})==1);
    for l = 1:length(indices)
        grp_subj(:,l) = all_subj(:,indices(l));
    end
    data_groups{k} = grp_subj;
end

FFR_avg_grpA = mean(data_groups{1,1},2);
FFR_avg_grpB = mean(data_groups{1,2},2);

%% Export mean FFRs into .txt files and convert into .avg files
fname_out_grpA = fullfile(OPTIONS.indir,strcat('mean_FFR_grpA_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_grpA,'w');
fprintf(fid,'%c\n',FFR_avg_grpA);
fclose(fid);

fname_out_grpB = fullfile(OPTIONS.indir,strcat('mean_FFR_grpB_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_grpB,'w');
fprintf(fid,'%c\n',FFR_avg_grpB);
fclose(fid);

fname_out_all = fullfile(OPTIONS.indir,strcat('mean_FFR_all_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_all,'w');
fprintf(fid,'%c\n',grd_FFR);
fclose(fid);

BT_toolbox_path = fullfile(pwd, strcat('ToolBox_Brainstem\BT_2013\programFiles')) ;
addpath(BT_toolbox_path) ;
bt_txt2avg(fname_out_grpA, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_grpB, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_all, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);

%% FFT : group analysis

% Adaptation of Skoe function (bt_fftsc)
%function [Freq1 Freq2 Freq3 fftFFR HzScale]=bt_fftsc(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% bt_fftsc computes frequency-domain amplitudes of three user-defined
% frequency bins of Brainstem response.  Results are scaled to peak �V.
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

FS = 16384;
FFR_list = {grd_FFR,FFR_avg_grpA,FFR_avg_grpB};
FFR_name = {'All_data','6-10mo','18-24mo'};
filename_list = {strcat('FFR_frequential_domain_all_subj_', OPTIONS.params,'.svg'),strcat('FFR_frequential_domain_groupA_', OPTIONS.params,'.svg'),strcat('FFR_frequential_domain_groupB_', OPTIONS.params,'.svg')};
color_list = {FFR_color;grpA_color;grpB_color};

%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD
for list_num = 1:3
    FFR = FFR_list{list_num};
    filename_export = filename_list{list_num};
    color = color_list{list_num};

    numPoints = length(FFR);

    %**** STEP 2. FFT OF FFR
    %******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
    rampMS = 4/1000; % length of ramp (on and off) in seconds
    hanPoints = rampMS.*FS; % length of ramp in points
    hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
    hanHalfPoints = round(hanPoints./2);
    numberOfOnes_A = numPoints - hanPoints;
    FFRhan_A = hann(hanPoints);
    FFRhan_A = [FFRhan_A(1:hanHalfPoints); ones(numberOfOnes_A,1); FFRhan_A(hanHalfPoints+1:hanPoints)];

    % baseline, window, then baseline again
    FFR = detrend(FFR, 'constant');
    FFR = detrend(FFR.*FFRhan_A, 'constant');

    %******** STEP 2b. Perform FFT
%     fftFFR = abs(fft(FFR, round(FS)));
    fftFFR = abs(fft(FFR(round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(1)/1000)*16384):round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(2)/1000)*16384)), round(FS)));
    fftFFR = fftFFR(1:round(round(FS)/2));
    fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
    HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
    HzScale = HzScale(1:length(fftFFR));

    % Plot FFT
    figure('Name',FFR_name{list_num}) ;

    subplot(2,2,1);
    plot(HzScale,fftFFR, 'Color', color);
    xlim([0 3000]);
    grid on;
    title("Single-Sided Amplitude Spectrum of X(t)");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude (µV)");

    subplot(2,2,2);
    plot(HzScale,fftFFR, 'Color', color);
    xlim([90 110]);
    grid on;
    title("Single-Sided Amplitude Spectrum of X(t)");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude (µV)");

    subplot(2,2,3);
    plot(HzScale,fftFFR, 'Color', color);
    xlim([300 500]);
    grid on;
    title("Single-Sided Amplitude Spectrum of X(t)");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude (µV)");

    subplot(2,2,4);
    plot(HzScale,fftFFR, 'Color', color);
    xlim([1100 1300]);
    grid on;
    title("Single-Sided Amplitude Spectrum of X(t)");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude (µV)");

    print('-dsvg',fullfile(OPTIONS.indir, filename_export));

end

%Display and save group comparison in frequency domain

FFR_A = FFR_avg_grpA;
FFR_B = FFR_avg_grpB ;

%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD
numPoints_A = length(FFR_A);
numPoints_B = length(FFR_B);

%**** STEP 2. FFT OF FFR
%******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
rampMS = 4/1000; % length of ramp (on and off) in seconds
hanPoints = rampMS.*FS; % length of ramp in points
hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
hanHalfPoints = round(hanPoints./2);
numberOfOnes_A = numPoints_A - hanPoints;
numberOfOnes_B = numPoints_B - hanPoints;
FFRhan_A = hann(hanPoints);
FFRhan_B = hann(hanPoints);
FFRhan_A = [FFRhan_A(1:hanHalfPoints); ones(numberOfOnes_A,1); FFRhan_A(hanHalfPoints+1:hanPoints)];
FFRhan_B = [FFRhan_B(1:hanHalfPoints); ones(numberOfOnes_B,1); FFRhan_B(hanHalfPoints+1:hanPoints)];

% baseline, window, then baseline again
FFR_A = detrend(FFR_A, 'constant');
FFR_A = detrend(FFR_A.*FFRhan_A, 'constant');
FFR_B = detrend(FFR_B, 'constant');
FFR_B = detrend(FFR_B.*FFRhan_B, 'constant');

%******** STEP 2b. Perform FFT
% fftFFR_A = abs(fft(FFR_A, round(FS)));
fftFFR_A = abs(fft(FFR_A(round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(1)/1000)*16384):round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(2)/1000)*16384)), round(FS)));
fftFFR_A = fftFFR_A(1:round(round(FS)/2));
fftFFR_A = fftFFR_A.*(2./numPoints_A); % scale to peak �V
HzScale_A = [0:1:round(FS/2)]'; % frequency 'axis'
HzScale_A = HzScale_A(1:length(fftFFR_A));

% fftFFR_B = abs(fft(FFR_B, round(FS)));
fftFFR_B = abs(fft(FFR_B(round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(1)/1000)*16384):round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(2)/1000)*16384)), round(FS)));
fftFFR_B = fftFFR_B(1:round(round(FS)/2));
fftFFR_B = fftFFR_B.*(2./numPoints_A); % scale to peak �V
HzScale_B = [0:1:round(FS/2)]'; % frequency 'axis'
HzScale_B = HzScale_B(1:length(fftFFR_B));

% Plot FFT in specific window of interest
figure('Name', 'Group_comparison_FFR_frequential') ;
plot(HzScale_A,fftFFR_A, 'Color', grpA_color); hold on;
plot(HzScale_B,fftFFR_B,'Color', grpB_color);
xlim([0 550]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

% Plot FFT in differents windows of interest
figure('Name', 'Group_comparison_FFR_frequential') ;

subplot(2,2,1);
plot(HzScale_A,fftFFR_A, 'Color', grpA_color); hold on;
plot(HzScale_B,fftFFR_B,'Color', grpB_color);
xlim([0 3000]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,2);
plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
plot(HzScale_B,fftFFR_B,'Color', grpB_color);
%xlim([90 110]);
xlim([95 105]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,3);
plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
plot(HzScale_B,fftFFR_B,'Color', grpB_color);
xlim([300 500]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,4);
plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
plot(HzScale_B,fftFFR_B,'Color', grpB_color);
xlim([1100 1300]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

%Add a single legend for 4 plots
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('6-10 mo', '18-24 mo','Location','bestoutside');
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.9;

%Add a single title for 4 plots
sgtitle('FFR Frequential domain, group comparison', 'Fontsize', 16, 'FontWeight', 'bold');

print('-dsvg',fullfile(OPTIONS.indir, strcat('FFR_frequential_domain_group_comparison_bis_', OPTIONS.params,'.svg')));

%% FFT : by subject analysis

% Adaptation of Skoe function (bt_fftsc)
%function [Freq1 Freq2 Freq3 fftFFR HzScale]=bt_fftsc(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% bt_fftsc computes frequency-domain amplitudes of three user-defined
% frequency bins of Brainstem response.  Results are scaled to peak �V.
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

all_F0 = zeros(size(all_subj,2),1) ;
all_ampl = zeros(size(all_subj,2),1) ;

% Loop through subjects to get FFR, F0 and amplitude of F0

for ss=1:size(all_subj,2)

    FFR = all_subj(:,ss) ;

    % baseline, window, then baseline again
    FFR = detrend(FFR, 'constant');
    FFR = detrend(FFR.*FFRhan_A, 'constant');

    % Perform FFT
    fftFFR = abs(fft(FFR(round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(1)/1000)*16384):round((abs(OPTIONS.win_of_interest(1))+OPTIONS.timew_F0(2)/1000)*16384)), round(FS)));
    %     fftFFR = abs(fft(FFR, round(FS)));
    fftFFR = fftFFR(1:round(round(FS)/2));
    fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
    HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
    HzScale = HzScale(1:length(fftFFR));

    % Plot FFT if set in options
    if OPTIONS.plot_FFT == 1
        figure('Name',subjects{ss}) ;

        plot(HzScale,fftFFR);
        xlim([0 1000]);
        grid on;
        title("Single-Sided Amplitude Spectrum of X(t)");
        xlabel("Frequency (Hz)");
        ylabel("Amplitude (µV)");
    end

    % Save F0 and amplitude
    if isempty(OPTIONS.woi_F0)
        A = max(fftFFR) ;             % get amplitude A of F0 (no window of interest)
    else
        A = max(fftFFR(OPTIONS.woi_F0(1):OPTIONS.woi_F0(2),1)) ;   % get amplitude A of F0 within window of interest
    end

    freq = find(fftFFR == A) ;    % get index of F0 in fftFFR list of values
    F = HzScale(freq);            % get value F of F0

    all_F0(ss,1) = F ;
    all_ampl(ss,1) = A ;

end

F0_and_ampl = table(subjects, all_F0, all_ampl,'VariableNames', {'subject','F0_Hz','amplitude_uV'}) ;

% Save table in a .csv file
if isempty(OPTIONS.woi_F0)
    fname = fullfile(OPTIONS.indir, strcat('all_subjects_F0_and_amplitude_no_woi_',num2str(OPTIONS.timew_F0(1)),'_',num2str(OPTIONS.timew_F0(2)),'tw_', OPTIONS.params,'.csv'));
else
    fname = fullfile(OPTIONS.indir, strcat('all_subjects_F0_and_amplitude_', num2str(OPTIONS.woi_F0(1)),'_', num2str(OPTIONS.woi_F0(end)),'woi_',num2str(OPTIONS.timew_F0(1)),'_',num2str(OPTIONS.timew_F0(2)),'tw_', OPTIONS.params,'.csv'));
end
writetable(F0_and_ampl,fname, 'WriteVariableNames', true) ;

%% Pitch error and pitch tracking

% Get inputed values. -------------------------------------------
block = 40 ;
step = 1 ;
startSTIM = 0 ;
endSTIM =  169;
channel =  1;
minFrequency = 0 ;
maxFrequency = 120 ;
minFrequency_stim = 80 ;  % 80
maxFrequency_stim = 120 ;   % 120

% set some more defaults and compute some more values from inputs
stim_channel = 1; % i.e. if stereo file, will only work on 1st channel
n_lags = readtable(fullfile(OPTIONS.indir,OPTIONS.nlag_filename));
%  % this allows path information to be retained from subject to subject
%  set(guifig, 'visible', 'off');

%% Run pitchtrack functions ------------------------------------------------------

autocorr_all = [] ; lag  = [] ; FFT_resp_all = [] ; prestimFFT_all = [] ; totalblocks_all = [] ;
time_stim_all = [] ; autocorr_stim_all = [] ; lag_stim = [] ; FFT_stim_all = [] ; totalblocks_stim_all = [] ;
LagStart_stim = [] ; LagStop_stim = [] ;
R_all = [] ;
FreqAC_all = [] ;
FreqFFT_stim_all = [] ;
FreqFFT_all = [] ;

for ss=1:size(all_subj,2)

    neural_lag = n_lags.neural_lag(ss);
    startRESP = startSTIM + neural_lag;
    % pitchtrack response
    [~, autocorr_all{ss}, lag{ss}, FFT_resp_all{ss}, ~, prestimFFT_all{ss}, totalblocks_all{ss}]= pitchtrack(fullfile(OPTIONS.indir,subjects{ss}, strcat(subjects{ss}, '_', OPTIONS.params, '_abr_avg_shifted_data_HF.avg')), block, step, startRESP, channel, 0) ;

    % pitchtrack stimulus (make conditional)
    [time_stim_all{ss}, autocorr_stim_all{ss}, lag_stim{ss}, FFT_stim_all{ss}, ~, ~, totalblocks_stim_all{ss}]=pitchtrack(OPTIONS.stim_avg, block, step, startSTIM, stim_channel,0) ;

    time_stim_all{ss} = time_stim_all{ss} + neural_lag ;  %for plotting purpose shift stimulus forward in time.

    [~, stopPT]=closestrc(time_stim_all{ss}, endSTIM+(block./2)) ;

    autocorr_stim = autocorr_stim_all{ss} ;
    autocorr_stim(:, stopPT+1:end)=[] ;
    time_stim = time_stim_all{ss} ;
    time_stim(stopPT+1:end)=[] ;
    FFT_stim = FFT_stim_all{ss} ;
    FFT_stim(:, stopPT+1:end)=[] ;
    totalblocks_stim_all{ss} = stopPT ;

    %Extract F0 from stimulus using Autocorrelation Method
    freqAC_vector = 1000./lag_stim{ss} ;
    [~, LagStart_stim{ss}]=closestrc(freqAC_vector,maxFrequency_stim) ;  %s is dummy variable
    [~, LagStop_stim{ss}] =closestrc(freqAC_vector,minFrequency_stim) ;
    [~, index] = max(autocorr_stim(LagStart_stim{ss}:LagStop_stim{ss}, :)) ;
    FreqAC_stim = freqAC_vector(LagStart_stim{ss}+index-1) ;
    clear index s

    %Extract F0 from response using Autocorrelation Method
    freqAC_vector = 1000./lag{ss} ;
    [~, LagStart{ss}]=closestrc(freqAC_vector,maxFrequency);  %s is dummy variable
    [~, LagStop{ss}] =closestrc(freqAC_vector,minFrequency);
    autocorr = autocorr_all{ss} ;
    [R_all{ss}, index] = max(autocorr(LagStart{ss}:LagStop{ss}, :));
    FreqAC_all{ss}=freqAC_vector(LagStart{ss}+index-1);
    LAG = 1000./FreqAC_all{ss};

    %Extract F0 from stimulus using FFT Method
    [~, index]= max(FFT_stim(minFrequency_stim+1:maxFrequency_stim+1, :));
    FreqFFT_stim_all{ss} = minFrequency_stim+index-1;
    FreqFFT_stim_all{ss} = FreqFFT_stim_all{ss}';

    %Extract F0 from response using FFT Method
    FFT_resp = FFT_resp_all{ss} ;
    [~, index]= max(FFT_resp(minFrequency+1:maxFrequency+1, :));
    FreqFFT_all{ss} = (minFrequency+index-1)'; %flip direction so that it matches FreqAC

end

%% Determine whether each extracted frequency was above the noise floor (NF)
% When the pitch-track is plotted a small gray dot will appear over the time
% ranges where the extracted frequency is below the noise floor The dot will be located above the plot and the location
% is based on what was inputted for maxFrequency (maxFrequency+4).  The total (i.e. total_belowNF) that gets exported is based
% only on the total number of blocks in the stimulus and not the total
% number of blocks in the response.  In other words, it only includes the blocks that were used to calculate pitch error.

PITCH_ERROR_autocorrelation = [] ;
PITCH_ERROR_fft = [] ;
PITCH_STRENGTH = [] ;
PITCH_STRENGTH2 = [] ;
PITCH_SRCORR = [] ;

for subnum=1:size(all_subj,2)

    % Initialize variables
    FFT_resp = FFT_resp_all{subnum} ;
    FreqAC =  FreqAC_all{subnum} ;
    totalblocks_stim = totalblocks_stim_all{subnum} ;
    prestimFFT = prestimFFT_all{subnum} ;
    FreqFFT = FreqFFT_all{subnum} ;
    F0_AMPusingACfreqs = 0;
    R = R_all{subnum} ;
    FreqFFT_stim = FreqFFT_stim_all{subnum} ;

    for x = 1:totalblocks_all{subnum}
        F0_AMPusingACfreqs(x,1) = FFT_resp(round(FreqAC(x))-1,x);  %Find amplitude using frequencies extracted using AC method
    end

    mean_F0_AMP = mean(F0_AMPusingACfreqs(1:totalblocks_stim,1));
    F0_AMP_prestim = prestimFFT(round(FreqAC)-1);
    PITCH_SNR = F0_AMPusingACfreqs./F0_AMP_prestim ;

    %Determine whether the extracted frequency was also the spectral maximum
    for x = 1:totalblocks_all{subnum}
        peaks{x}=localmax(FFT_resp(:,x));  %Find amplitude using frequencies extracted using AC method
        [closest_spectralpeak(x), row(x), column(x)]=closestrc(peaks{x}, FreqAC(x));
        if  FreqFFT(x)== (closest_spectralpeak(x)-1)
            notspectralMax(x) = 0;
        else
            notspectralMax(x) = 1;
        end
    end

    % When the pitch-track is plotted a small gray dot will appear over the time
    % ranges where the extracted frequency is not at the spectral max. The dot will be located above the plot and the location
    % is based on what was inputted for maxFrequency.  The total (i.e. total_notatSpectralMax) that gets exported is based
    % only on the total number of blocks in the stimulus and not the total number of blocks in the response.  In other words, \
    % it only includes the blocks that were used to calculate pitch error.
    total_notatSpectralMax = sum(notspectralMax(1:totalblocks_stim));
    plot_notatspectralMax(notspectralMax==1)=maxFrequency+5;
    plot_notatspectralMax(notspectralMax==0)=NaN;

    plot_belowNF(PITCH_SNR<1)=maxFrequency+4;
    plot_belowNF(PITCH_SNR>=1)=NaN;

    % Calculate Final Measures:
    % Autocorrelation
    PITCH_ERROR_autocorrelation{subnum} = mean(abs(FreqAC(1:totalblocks_stim)-FreqAC_stim));  %Measured in Hz.
    %FFT
    PITCH_ERROR_fft{subnum} = mean(abs(FreqFFT(1:totalblocks_stim)-FreqFFT_stim));  %Measured in Hz.

    %on the off chance that one or more of the Rs is exactly 1, we must set
    %these Rs to 0.999999 to get a valid number for fisher (i.e. not inf)If you are are reading this you are as annoyed as we are.
    PITCH_STRENGTH{subnum} = mean(R(1:totalblocks_stim));
    Rtemp=R;
    Rtemp(Rtemp==1)=0.999999;
    PITCH_STRENGTH2{subnum} = fisherinv(mean(fisher(Rtemp(1:totalblocks_stim))));
    CORR = corrcoef(FreqAC(1:totalblocks_stim), FreqAC_stim);  %correlation between stimulus and response f0 contour
    PITCH_SRCORR{subnum} =CORR(1,2); %the first number is always 1, need to take second
    total_belowNF =  sum(PITCH_SNR(1:totalblocks_stim)<1);
    total_notatSpectralMax = sum(notspectralMax(1:totalblocks_stim));

end

%Extract variables of interest in table and save table
pitch_error_infos = table(subjects,PITCH_ERROR_autocorrelation',PITCH_ERROR_fft',PITCH_STRENGTH2',PITCH_SRCORR', ...
    'VariableNames',{'subject_ID', 'PITCH_ERROR_autocorrelation','PITCH_ERROR_fft','PITCH_STRENGTH_2','PITCH_SRCORR'}) ;
writetable(pitch_error_infos,fullfile(OPTIONS.indir,strcat('pitch_error_group_data_', OPTIONS.params,'.csv'))) ;

%% Response-to-response correlation

%compare grpA and grpB

%% Response consistency

%compare 2 subaverages of the same FFR (ucse all subjects from eacg group)


