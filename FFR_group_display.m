function [] = FFR_group_display(subjects, OPTIONS) 
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
%for loopnum=find(ismember(subjects,'DVL_003_T10')) ; 
    FFR_file = fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.txt')) ; 
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR preprocessing steps on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

 
%% Time domain : FFR visualization

% Compute grand average and plot
grd_FFR = mean(all_subj,2);
figure ; 
plot(timepoints,grd_FFR,'Color',FFR_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR, 6-24 mo');

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

% Visualization by group

% Temporal
fig = figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');

%Print plot into .svg, png and fig files
if OPTIONS.savefig == 1
    print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grp_comparison_', OPTIONS.params,'.svg')));
    print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grp_comparison_', OPTIONS.polarity,'_FFR_temporal_', OPTIONS.params,'.png'))); 
    saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grp_comparison_', OPTIONS.polarity,'_FFR_temporal_', OPTIONS.params,'.fig')));
end

% Visualization of each group
fig = figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 6-10mo');

if OPTIONS.savefig == 1
    print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grpA_', OPTIONS.params,'.svg')));
    print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grpA_', OPTIONS.polarity,'_', OPTIONS.params,'.png')));
    saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grpA_', OPTIONS.polarity,'_', OPTIONS.params,'.fig')));
end

fig = figure ; 
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 18-24mo');

print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grpB_', OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grpB_', OPTIONS.polarity,'_', OPTIONS.params,'.png')));
saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grpB_', OPTIONS.polarity,'_', OPTIONS.params,'.fig')));

%% Neural lag

% Estimation of the transmission delay between stimulus and response. 
% Calculated from the time lag that produces the maximum stimulus-to-response cross-correlation magnitude

% Read files that contains neural lag and age information
neural_lags = readtable(fullfile(OPTIONS.indir, OPTIONS.nlag_filename)) ;
age_in_days = readtable(fullfile(OPTIONS.indir, 'age_in_days.csv'), 'Delimiter',';') ;

% Keep only subjects of interest (not rejected)
neural_lags = neural_lags(contains(neural_lags.suject_ID, subjects),:) ;
age_in_days = age_in_days(contains(age_in_days.subjects, subjects),:) ;

%Get variables of interest: neural lags and ages 
IDlist_grA = neural_lags.suject_ID(contains(neural_lags.group,'A')) ;
IDlist_grB = neural_lags.suject_ID(contains(neural_lags.group,'B')) ;
neural_lags_grA = neural_lags.neural_lag(contains(neural_lags.group,'A')) ;
neural_lags_grB = neural_lags.neural_lag(contains(neural_lags.group,'B')) ;
age_grA = age_in_days.age_in_days(contains(age_in_days.subjects,IDlist_grA)) ;
age_grB = age_in_days.age_in_days(contains(age_in_days.subjects,IDlist_grB)) ;

all_info = table(neural_lags.suject_ID, neural_lags.neural_lag, neural_lags.group, age_in_days.age_in_days) ;

 % Display neural lag distribution as a function of age
figure ; scatter(age_grA, neural_lags_grA) ; hold on ; scatter(age_grB, neural_lags_grB) ; legend({'6-10 mo', '18-24mo'}) ;
figure ; boxplot(all_info.Var2, all_info.Var3, 'Notch','on','Labels',{'6-10 mo', '18-24mo'}) ;

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

if OPTIONS.savefig == 1
    print('-dsvg',fullfile(OPTIONS.indir, strcat('FFR_frequential_domain_group_comparison_bis_', OPTIONS.params,'.svg')));
end
