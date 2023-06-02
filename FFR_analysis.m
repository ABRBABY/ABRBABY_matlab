function [] = FFR_analysis(subjects, OPTIONS) 
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
    FFR_file = fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},OPTIONS.param,'_abr_shifted_data_HF.txt')) ; 
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
figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');

%Print plot into .svg pictures
print('-dsvg',fullfile(OPTIONS.indir,'mean_FFR_grp_comparison.svg'));

% Visualization of each group
figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 6-10mo');

print('-dsvg',fullfile(OPTIONS.indir,'mean_FFR_grpA.svg'));

figure ; 
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 18-24mo');

print('-dsvg',fullfile(OPTIONS.indir,'mean_FFR_grpB.svg'));

%% Export mean FFRs into .txt files and convert into .avg files
fname_out_grpA = fullfile(OPTIONS.indir,'mean_FFR_grpA.txt') ;
fid = fopen(fname_out_grpA,'w');
fprintf(fid,'%c\n',FFR_avg_grpA);
fclose(fid);

fname_out_grpB = fullfile(OPTIONS.indir,'mean_FFR_grpB.txt') ;
fid = fopen(fname_out_grpB,'w');
fprintf(fid,'%c\n',FFR_avg_grpB);
fclose(fid);

fname_out_all = fullfile(OPTIONS.indir,'mean_FFR_all.txt') ;
fid = fopen(fname_out_all,'w');
fprintf(fid,'%c\n',grd_FFR);
fclose(fid);

BT_toolbox_path = fullfile(pwd, strcat('ToolBox_Brainstem\BT_2013\programFiles')) ;
addpath(BT_toolbox_path) ;
bt_txt2avg(fname_out_grpA, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_grpB, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_all, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);

%% Root mean square and SNR on averaged FFRs

% Adapted from bt_rms.m from bt_gui toolbox (Kraus & Skoe, 2010)
% Note: std(X,1) is a shortcut RMS calc. It is equivalent to RMS
% on a baselined (demeaned to 0) waveform.  This is what we want here.

%note : the mean FFR includes the [-40;0] ms prestimulus period. In consequence:
%[0;40] = prestimulus period (->[-40;0]);
%[40;210] = whole FFR(->[0;170]);
%[50-95] = consonant-vowel transition period(->[10;55]);
%[95-210] = steady-state vowel period(->[55;170]);

startPrestim = 0.0610; % (in ms) to obtain 1 in samples
endPrestim = 40; % (in ms)

startCVtransition = 50;
endCVtransition = 95;

startVowel = 95;
endVowel = 210;

% Prestim period
Segment_A = FFR_avg_grpA(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000));
Segment_B = FFR_avg_grpB(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000));
RMSprestimA = std(Segment_A,1);
RMSprestimB = std(Segment_B,1);
%RMSprestimA = std(FFR_avg_grpA(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000)),1);
%RMSprestimB = std(FFR_avg_grpB(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000)),1);

% Whole response period
RMS_whole_grpA = std(FFR_avg_grpA(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
RMS_whole_grpB = std(FFR_avg_grpB(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
%RMS_grpA = std(FFR_avg_grpA,1); 
%RMS_grpB = std(FFR_avg_grpB,1); 

%CV transition period
RMS_CV_grpA = std(FFR_avg_grpA(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000)),1);
RMS_CV_grpB = std(FFR_avg_grpB(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000)),1);

%Vowel period
RMS_vowel_grpA = std(FFR_avg_grpA(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
RMS_vowel_grpB = std(FFR_avg_grpB(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);

% Signal-to-noise on whole response
%SNR is obtained by dividing the RMS of a specific time window of the
%response by the RMS of the prestimulus response
% SNR_grpA = RMS_grpA./RMSprestimA;
% SNR_grpB = RMS_grpB./RMSprestimB;
%SNR_whole_grpA = RMS_whole_grpA./RMSprestimA;
%SNR_whole_grpB = RMS_whole_grpB./RMSprestimB;
SNR_whole_grpA = RMS_whole_grpA./RMSprestimA;
SNR_whole_grpB = RMS_whole_grpB./RMSprestimB;

% Signal-to-noise on CV transition period
SNR_CV_grpA = RMS_CV_grpA./RMSprestimA;
SNR_CV_grpB = RMS_CV_grpB./RMSprestimB;

% Signal-to-noise on vowel period
SNR_vowel_grpA = RMS_vowel_grpA./RMSprestimA;
SNR_vowel_grpB = RMS_vowel_grpB./RMSprestimB;

% display results
%disp(["Group A: RMS total = ",RMS_grpA, "RMS prestim = ", RMSprestimA, "SNR = ", SNR_grpA]);
%disp(["Group B: RMS total = ",RMS_grpB, "RMS prestim = ", RMSprestimB, "SNR = ", SNR_grpB]);

Group = ["6-10 mo";"18-24 mo"];
RMS_Prestim = [RMSprestimA;RMSprestimB];
RMS_Total = [RMS_whole_grpA;RMS_whole_grpB];
RMS_CV = [RMS_CV_grpA;RMS_CV_grpB];
RMS_Vowel = [RMS_vowel_grpA;RMS_vowel_grpB];
SNR_Total = [SNR_whole_grpA;SNR_whole_grpB];
SNR_CV = [SNR_CV_grpA;SNR_CV_grpB];
SNR_Vowel = [SNR_vowel_grpA;SNR_vowel_grpB];

RMS_and_SNR = table(Group,RMS_Prestim,RMS_Total,RMS_CV,RMS_Vowel,SNR_Total,SNR_CV,SNR_Vowel);
RMS_and_SNR

% Display distribution of RMS and SNR
figure ; scatter(age_grA, neural_lags_grA) ; hold on ; scatter(age_grB, neural_lags_grB) ; legend({'6-10 mo', '18-24mo'}) ;
figure ; boxplot(all_info.Var2, all_info.Var3, 'Notch','on','Labels',{'6-10 mo', '18-24mo'}) ;

% Save table in a .csv file
fname = fullfile(OPTIONS.indir, 'RMS_and_SNR_group_comparison.csv');
writetable(RMS_and_SNR,fname, 'WriteVariableNames', true) ;

% %% Root mean square and SNR on indivudual FFRs
% 
% Subj = [subjects(:)];
% 
% RMS_Prestim_all = zeros(1,length(subjects));
% RMS_Total_all = zeros(1,length(subjects));
% RMS_CV_all = zeros(1,length(subjects));
% RMS_Vowel_all = zeros(1,length(subjects));
% SNR_Total_all = zeros(1,length(subjects));
% SNR_CV_all = zeros(1,length(subjects));
% SNR_Vowel_all = zeros(1,length(subjects));
% 
% for s =1:length(subjects)
%     RMS_Prestim_all(s) = std(all_subj(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000),s),1);
%     RMS_Total_all(s) = std(all_subj(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000),s),1);
%     RMS_CV_all(s) = std(all_subj(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000),s),1);
%     RMS_Vowel_all(s) = std(all_subj(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000),s),1);
%     SNR_Total_all(s) = RMS_Total_all(s)/RMS_Prestim_all(s);
%     SNR_CV_all(s) = RMS_CV_all(s)/RMS_Prestim_all(s);
%     SNR_Vowel_all(s) = RMS_Vowel_all(s)/RMS_Prestim_all(s);
% end
% 
% RMS_Prestim_all = RMS_Prestim_all';
% RMS_Total_all = RMS_Total_all';
% RMS_CV_all = RMS_CV_all';
% RMS_Vowel_all = RMS_Vowel_all';
% SNR_Total_all = SNR_Total_all';
% SNR_CV_all = SNR_CV_all';
% SNR_Vowel_all = SNR_Vowel_all';
% 
% RMS_and_SNR_all = table(Subj,RMS_Prestim_all,RMS_Total_all,RMS_CV_all,RMS_Vowel_all,SNR_Total_all,SNR_CV_all,SNR_Vowel_all);
% RMS_and_SNR_all
% 
% %todo : organiser data par groupe d'âge
% 

%% Neural lag

% Estimation of the transmission delay between stimulus and response. 
% Calculated from the time lag that produces the maximum stimulus-to-response cross-correlation magnitude

% Read files that contains neural lag and age information
neural_lags = readtable(fullfile(OPTIONS.indir, 'all_neural_lags.csv')) ;
age_in_days = readtable(fullfile(OPTIONS.indir, 'age_in_days.xlsx')) ;

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

 %% Pitch tracking

 %function [time autocorr LAG FFT freqAxis preFFT blocks]= pitchtrack(avgname, block, step, startAnalysis,channel,exportData)
% 
% % block = 1;
% % step = 1;
% % startAnalysis = 1;
% % %channel = 1;
% % %exportData = 1;
% % file.xmin =1; 
% % file.xmax=240;
% % file.pnts=2932;
% % 
% % % if block<40
% % %     display('Block Size is too small. Using default: 40 ms');
% % %     block = 40;
% % % end
% % 
% % % %Open average file 
% % % [file]= openavg(avgname);
% % %Define time axis 
% % timeaxis = linspace(file.xmin, file.xmax, file.pnts)';
% % 
% %  % extract signal
% %  %SIGNAL = file.signal(:,channel);
% %  SIGNAL = FFR_avg_grpB;
% %  % get sampling rate
% %  fs = 16384;
% % 
% %  
% %  % ---------------------------PRESTIM, ---------------------------
% %  % extract portion of prestimulus time, and detrend. The size of the prestim portion is dependent on the block size.
% %  % To do SNR calculations block and prestim must be same the same number of ms.
% %  % If the block size is 40ms, then only 40ms of the prestim will be
% %  % extracted. 
% %  
% %  %PRESTIM = SIGNAL(1: ms2row(file, block));   %start with the very first point.
% %  PRESTIM = Segment_B;
% %  % ramp
% %  r = hann(size(PRESTIM, 1));  % the entire prestim is ramped.
% %  % ramp and detrend
% %  PRESTIM = detrend(PRESTIM.*hann(size(PRESTIM,1)), 'constant');
% %  % FFT (zero-padded to sampling rate);
% %  preFFT = abs(fft(PRESTIM, fs));
% %  %scale preFFT
% %  preFFT= preFFT*(2./length(PRESTIM));
% %  preFFT=preFFT(1:1001,1);  %truncate above 1000 Hz
% %  
% %    
% % %  ------------------ FFTS of RESPONSE CHUNKS-----------------------
% % j = startAnalysis; % each time through loop j increases by step size;
% % 
% %     chunks = 5000;    % an arbitrary maximum number of blocks that the program will create. 
% %                         
% %     for k = 1:chunks;   %the program knows to stop once file.xmax is exceeded 
% % 
% %         % variables created: 
% %         ramptime = (block/1000);   % ramp the entire chunk
% %         start = j;
% %         stop = j+block;
% % 
% %         if stop>(file.xmax)  % if stop exceeds the maximum ms time then abort and break out from loop
% %             k=k-1;
% %             j=j-step;
% %             break;
% %         else
% %             signal = detrend(SIGNAL(ms2row(file, start):ms2row(file, stop)), 'constant');   % de-mean to zero           
% %         end
% %         
% %         midpoint(k) = mean(ms2row(file, start):ms2row(file,stop));  %calculates the time corresponding to the midpoint of the chunk
% % 
% %         % generate ramp
% %         ramp = hann(size(signal,1));
% %         % ramp and de-mean
% %         signal = detrend(signal.*ramp, 'constant');
% % 
% % 
% %         % autocorrelation (see Boersma 1993)
% %         [c lag]=xcorr(signal, 'coeff');
% %        
% %         
% %         [cwin lagwin]=xcorr(ramp, 'coeff');
% %         LAG = linspace(-block, block, length(c));
% %     
% %         autoc=c./cwin;
% %         autoc(autoc>1)=1;  %this handles the very rare case that the remainder of the previous step is >1.
% %         % only plot the first 15 ms; % lowest frequency is ~66 Hz.
% %         
% %         
% %             
% %             startlag = find(LAG==0);  
% %             endlag = find(LAG==closestrc(LAG, 15));
% % 
% %         
% %         % truncate lag and r value matrices to only include values up to first 15 ms.
% %         autocorr(:,k)=autoc(startlag:endlag)';
% %         LAG = LAG(startlag:endlag)';
% %      
% %         ostartlag = startlag;
% %         
% %         ostoplag = endlag;
% %          
% %          %%% Now do FFTs;
% %         % fft, pads to sampling rate
% %         fftsignal{k} = abs(fft(signal, fs));
% %         % only go up to 1000 Hz;
% %         FFT{k} = fftsignal{k}(1:1001,1);
% %         FFT{k}= FFT{k}*(2./length(signal));
% %         freqAxis = linspace(0, 1000, 1001);
% %         
% %         j = j+step;  % loop through next time chunk
% %         
% %        
% %        
% %     end
% %     time = timeaxis(round(midpoint));
% %     
% %     blocks = k;
% %     
% %     FFT = cell2mat(FFT);
%     
%     
% %     if exportData == 1
% %         [fpath fname ext]=fileparts(avgname);
% %         FFTfile = [fpath, '\', fname, '-FFTmatrix.xls'];
% %         ACfile = [fpath, '\', fname, '-ACmatrix.xls'];
% %         
% %         xlswrite(FFT, fname, {'FFT matrix'}, FFTfile, 'Sheet1');
% %         xlswrite(freqAxis', fname, {'Frequency Axis'}, FFTfile, 'Sheet2');
% %         xlswrite(time, fname, {'Time Axis'}, FFTfile,     'Sheet3');
% %         
% %         xlswrite(autocorr, fname, {'autocorrelation matrix'}, ACfile, 'Sheet1');
% %         xlswrite( LAG, fname, {'Lag Axis'} , ACfile, 'Sheet2');
% %         xlswrite(time, fname, {'Time Axis'} , ACfile, 'Sheet3');
% %        
% %      
% %     end
% %     
% %% Pitch error
% 
% 
% %% Response-to-response correlation    
% 
% %compare grpA and grpB
% 
% %% Response consistency 
% 
% %compare 2 subaverages of the same FFR (use all subjects from eacg group)
% 
% %% Frequency domain
% 
% % Adaptation of Skoe function (bt_fftsc)
%  %function [Freq1 Freq2 Freq3 fftFFR HzScale]=bt_fftsc(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% % bt_fftsc computes frequency-domain amplitudes of three user-defined 
% % frequency bins of Brainstem response.  Results are scaled to peak �V.
% %
% % Usage: [F0 F1 HF] = bt_fftsc('filename.avg',10,40,100,150,300,350,600,800);
% %    over the range of 10 to 40 ms, finds average frequency amplitude of
% %    100-150 Hz, 300-350 Hz and 600-800 Hz.
% % 
% % Three variables are returned the workspace:
% % (1) Freq1: mean amplitude over F0_Lo-F0_Hi Hz
% % (2) Freq2: mean amplitude over F1_Lo-F1_Hi Hz
% % (3) Freq3: mean amplitude over HF_Lo to HF_Hi Hz 
% 
% % % Dependancies: ms2row, openavg
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Originally developed by E.E. Skoe.  
% % Toolbox version by E.E. Skoe & T.G. Nicol
% % eeskoe@northwestern.edu tgn@northwestern.edu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% F0_Lo = 100;
% F0_Hi = 110;
% F1_Lo = 455;
% F1_Hi = 740;
% HF_Lo = 1060;
% HF_Hi = 2750;
% 
% %%%  OPEN UP AVG FILE
% %avg = openavg(FILENAME);
% %data = grd_FFR;
% FS = 16384;
% %FFR = grd_FFR;
% FFR_list = {grd_FFR,FFR_avg_grpA,FFR_avg_grpB};
% FFR_name = {'All_data','6-10mo','18-24mo'};
% filename_list = {'FFR_frequential_domain_all_subj.svg','FFR_frequential_domain_groupA.svg','FFR_frequential_domain_groupB.svg'};
% color_list = {FFR_color;grpA_color;grpB_color};
% %******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD 
% for list_num = 1:3
%     FFR = FFR_list{list_num};
%     filename_export = filename_list{list_num};
%     color = color_list{list_num};
% %startPoint = ms2row(avg, start);
% %endPoint = ms2row(avg, stop);
% 
% %FFR = data(startPoint:endPoint);
% numPoints = length(FFR);
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% %**** STEP 2. FFT OF FFR
% %******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
% rampMS = 4/1000; % length of ramp (on and off) in seconds
% % hanPoints = 26;  %hard coded to be the same as Biologic's settings (December 7, 2005);
% hanPoints = rampMS.*FS; % length of ramp in points
% hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
% hanHalfPoints = round(hanPoints./2);
% numberOfOnes_A = numPoints - hanPoints;
% FFRhan_A = hann(hanPoints);  
% FFRhan_A = [FFRhan_A(1:hanHalfPoints); ones(numberOfOnes_A,1); FFRhan_A(hanHalfPoints+1:hanPoints)];
% 
% % baseline, window, then baseline again
% FFR = detrend(FFR, 'constant');
% FFR = detrend(FFR.*FFRhan_A, 'constant');
% 
% %******** STEP 2b. Perform FFT
% fftFFR = abs(fft(FFR, round(FS)));
% fftFFR = fftFFR(1:round(round(FS)/2));
% fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
% HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
% HzScale = HzScale(1:length(fftFFR));
% 
% % Plot FFT
% figure('Name',FFR_name{list_num}) ;
% 
% subplot(2,2,1);
% plot(HzScale,fftFFR, 'Color', color);
% xlim([0 3000]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,2); 
% plot(HzScale,fftFFR, 'Color', color);
% xlim([90 110]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,3); 
% plot(HzScale,fftFFR, 'Color', color);
% xlim([300 500]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,4); 
% plot(HzScale,fftFFR, 'Color', color);
% xlim([1100 1300]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% print('-dsvg',fullfile(OPTIONS.indir, filename_export));
% 
% end
% 
% %Display and save group comparison in frequency domain
% 
% FFR_A = FFR_avg_grpA;
% FFR_B = FFR_avg_grpB ;
% 
% %******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD
% numPoints_A = length(FFR_A);
% numPoints_B = length(FFR_B);
% 
% %**** STEP 2. FFT OF FFR
% %******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
% rampMS = 4/1000; % length of ramp (on and off) in seconds
% hanPoints = rampMS.*FS; % length of ramp in points
% hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
% hanHalfPoints = round(hanPoints./2);
% numberOfOnes_A = numPoints_A - hanPoints;
% numberOfOnes_B = numPoints_B - hanPoints;
% FFRhan_A = hann(hanPoints);
% FFRhan_B = hann(hanPoints);
% FFRhan_A = [FFRhan_A(1:hanHalfPoints); ones(numberOfOnes_A,1); FFRhan_A(hanHalfPoints+1:hanPoints)];
% FFRhan_B = [FFRhan_B(1:hanHalfPoints); ones(numberOfOnes_B,1); FFRhan_B(hanHalfPoints+1:hanPoints)];
% 
% 
% % baseline, window, then baseline again
% FFR_A = detrend(FFR_A, 'constant');
% FFR_A = detrend(FFR_A.*FFRhan_A, 'constant');
% FFR_B = detrend(FFR_B, 'constant');
% FFR_B = detrend(FFR_B.*FFRhan_B, 'constant');
% 
% %******** STEP 2b. Perform FFT
% fftFFR_A = abs(fft(FFR_A, round(FS)));
% fftFFR_A = fftFFR_A(1:round(round(FS)/2));
% fftFFR_A = fftFFR_A.*(2./numPoints_A); % scale to peak �V
% HzScale_A = [0:1:round(FS/2)]'; % frequency 'axis'
% HzScale_A = HzScale_A(1:length(fftFFR_A));
% 
% fftFFR_B = abs(fft(FFR_B, round(FS)));
% fftFFR_B = fftFFR_B(1:round(round(FS)/2));
% fftFFR_B = fftFFR_B.*(2./numPoints_A); % scale to peak �V
% HzScale_B = [0:1:round(FS/2)]'; % frequency 'axis'
% HzScale_B = HzScale_B(1:length(fftFFR_B));
% 
% % Plot FFT
% figure('Name', 'Group_comparison_FFR_frequential') ;
% 
% subplot(2,2,1);
% plot(HzScale_A,fftFFR_A, 'Color', grpA_color); hold on;
% plot(HzScale_B,fftFFR_B,'Color', grpB_color);
% xlim([0 3000]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,2); 
% plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
% plot(HzScale_B,fftFFR_B,'Color', grpB_color);
% %xlim([90 110]);
% xlim([95 105]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,3); 
% plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
% plot(HzScale_B,fftFFR_B,'Color', grpB_color);
% xlim([300 500]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% subplot(2,2,4); 
% plot(HzScale_A,fftFFR_A,'Color', grpA_color); hold on;
% plot(HzScale_B,fftFFR_B,'Color', grpB_color);
% xlim([1100 1300]);
% grid on;
% title("Single-Sided Amplitude Spectrum of X(t)");
% xlabel("Frequency (Hz)");
% ylabel("Amplitude (µV)");
% 
% %Add a single legend for 4 plots
% fig = gcf;
% fig.Position(3) = fig.Position(3) + 250;
% Lgnd = legend('6-10 mo', '18-24 mo','Location','bestoutside');
% Lgnd.Position(1) = 0.01;
% Lgnd.Position(2) = 0.9;
% 
% %Add a single title for 4 plots
% sgtitle('FFR Frequential domain, group comparison', 'Fontsize', 16, 'FontWeight', 'bold');
% 
% print('-dsvg',fullfile(OPTIONS.indir, 'FFR_frequential_domain_group_comparison_bis.svg'));
% 
% %%
