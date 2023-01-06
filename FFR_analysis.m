function [] = FFR_analysis(indir) 
% ERPs analysis script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in INDIR 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Color for plot
DIFF_color = [0 0 0]; %black
FFR_color = [0.8902 0 0]; %red
grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.2 0.7765 0.2]; %green

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
%for loopnum=find(ismember(subjects,'DVL_003_T10')) ; 
    FFR_file = fullfile(indir,subjects{loopnum},strcat(subjects{loopnum},'_abr_shifted_data_HF.txt')) ; 
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR sanity check on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

% % Plot individual FFRs
% for jj = 1:size(subjects,1)
%     figure ; 
%     plot(timepoints,all_subj(:,jj),'b','Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
%     grid on ; 
%     legend('Individual FFR', subjects(jj), 'Interpreter', 'None');
%     xlabel('Times (ms)'); ylabel('uV'); title ([' FFR ', subjects(jj)], 'Interpreter', 'None');
% end

% % Exclude noisy participants from observation
% a = ismember(subjects,noise_list)==1
% for excl = 1:size(subjects,1)
%     if a(excl)==1
%         subjects(excl) = [];
%         all_subj(:,excl) = []; 
%     end
% end

% Exclude noisy participants from max and min values
excluded_subj = {};
ex = 1;
for noise = 1:size(all_subj,2)
    if max(all_subj(:,noise)) > 0.5 || min(all_subj(:,noise)) < -0.5
        excluded_subj(ex) = subjects(noise,1);   
        disp([subjects(noise,1), 'excluded']);
        ex = ex +1;
    end
   
end

% Delete excluded subjects
exclud = ismember(subjects,excluded_subj);
new_all_subj = zeros(size(timepoints,1),1);
new_subjects = {};
n = 1;
for tt = 1:size(exclud,1)
    if exclud(tt)==0
        new_all_subj(:,n) = all_subj(:,tt);
        new_subjects(n) = subjects(tt);
        n = n+1;
    end
end

  
%%
% Compute grand average and plot
grd_FFR = mean(new_all_subj,2);
figure ; 
plot(timepoints,grd_FFR,'Color',FFR_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR, 6-24 mo');

% Classify into groups
grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};

% Compute grand average by group
grpA.subj = new_subjects(contains(new_subjects,grpA.suffix));
grpB.subj = new_subjects(contains(new_subjects,grpB.suffix));

grpA.data = zeros(size(timepoints,1),1);
grpB.data = zeros(size(timepoints,1),1);

groups = {grpA.subj, grpB.subj};
data_groups = {grpA.data, grpB.data};
for k = 1:length(groups)
    grp_subj = zeros(size(timepoints,1),1);
    indices = find(ismember(new_subjects,groups{k})==1);
    for l = 1:length(indices)
        grp_subj(:,l) = new_all_subj(:,indices(l));
    end
    data_groups{k} = grp_subj;
end

FFR_avg_grpA = mean(data_groups{1,1},2);
FFR_avg_grpB = mean(data_groups{1,2},2);

% Visualization by group

% Temporal
figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');

%%
% Frequential
%TODO


% Fs = 16384;            % Sampling frequency                    
% T = 1/Fs;              % Sampling period       
% L = 250;               % Length of signal
% t = (0:L-1)*T;         % Time vector
% 
% Y = fft(grd_FFR);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")


%% Adaptation of Skoe function (bt_fftsc)
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
%avg = openavg(FILENAME);
%data = grd_FFR;
FS = 16384;
FFR = grd_FFR;
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
fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
HzScale = HzScale(1:length(fftFFR));

%% Plot FFT
figure ;

subplot(2,2,1);
plot(HzScale,fftFFR);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,2); 
plot(HzScale,fftFFR);
xlim([90 110]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,3); 
plot(HzScale,fftFFR);
xlim([300 500]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,4); 
plot(HzScale,fftFFR);
xlim([1100 1300]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

% %% This section does not work for now
% % clear variables no longer needed
% clear FFRramp rampMS hanPoints hanHalfPoints numberOfOnes FFRhan
% 
% %**** STEP 3. compute mean magnitudes over F0, F-1 and HF ranges
% % i. F0: F0_Lo-F0_Hi.
% % find freqs nearest F0_Lo and F0_Hi.
% startF = find(F0_Lo-1/2 < HzScale & HzScale < F0_Lo+1/2); % 1 is stepsize
% stopF = find(F0_Hi-1/2 < HzScale & HzScale < F0_Hi+1/2);
% % compute mean
% Freq1 = mean(fftFFR(startF:stopF,1));
% % plot
% figure ; 
% plot(HzScale(startF:stopF),Freq1);
% grid on;
% title("Freq1");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (�V)");
% 
% % ii. F1: F1_Lo-F1_Hi
% % find freqs nearest F1_Lo and F1_Hi.
% startF = find(F1_Lo-1/2 < HzScale & HzScale < F1_Lo+1/2); % 1 is stepsize
% stopF = find(F1_Hi-1/2 < HzScale & HzScale < F1_Hi+1/2);
% % compute mean
% Freq2 = mean(fftFFR(startF:stopF,1));
% % plot
% figure ; 
% plot(HzScale,Freq2);
% grid on;
% title("Freq2");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (�V)");
% 
% % iii. HF: HF_Lo-HF_Hi
% % find freqs nearest HF_Lo and HF_Hi.
% startF = find(HF_Lo-1/2 < HzScale & HzScale < HF_Lo+1/2); % 1 is stepsize
% stopF = find(HF_Hi-1/2 < HzScale & HzScale < HF_Hi+1/2);
% % compute mean
% Freq3 = mean(fftFFR(startF:stopF,1));
% % plot
% figure ; 
% plot(HzScale,Freq3);
% grid on;
% title("Freq3");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (�V)");

