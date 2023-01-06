function [Freq1 Freq2 Freq3 fftFFR]=bt_fft2(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% bt_fft2 computes frequency-domain amplitudes of three user-defined 
% frequency bins of Brainstem response.  Results are not scaled to peak µV.
%
% Usage: [F0 F1 HF] = bt_fft2('filename.avg',10,40,100,150,300,350,600,800);
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

%%%  OPEN UP AVG FILE
avg = openavg(FILENAME);
data = avg.signal(:,chan);
FS = avg.rate;

%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD 

startPoint = ms2row(avg, start);
endPoint = ms2row(avg, stop);

FFR = data(startPoint:endPoint);
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
HzScale = [0:1:round(FS/2)-1]'; % frequency 'axis'

% clear variables no longer needed
clear FFRramp rampMS hanPoints hanHalfPoints numberOfOnes FFRhan

%**** STEP 3. compute mean magnitudes over F0, F1 and HF ranges
% i. F0: F0_Lo-F0_Hi.
% find freqs nearest F0_Lo and F0_Hi.
startF = find(F0_Lo-1/2 < HzScale & HzScale < F0_Lo+1/2); % 1 is stepsize
stopF = find(F0_Hi-1/2 < HzScale & HzScale < F0_Hi+1/2);
% compute mean
Freq1 = mean(fftFFR(startF:stopF,1));

% ii. F1: F1_Lo-F1_Hi
% find freqs nearest F1_Lo and F1_Hi.
startF = find(F1_Lo-1/2 < HzScale & HzScale < F1_Lo+1/2); % 1 is stepsize
stopF = find(F1_Hi-1/2 < HzScale & HzScale < F1_Hi+1/2);
% compute mean
Freq2 = mean(fftFFR(startF:stopF,1));

% iii. HF: HF_Lo-HF_Hi
% find freqs nearest HF_Lo and HF_Hi.
startF = find(HF_Lo-1/2 < HzScale & HzScale < HF_Lo+1/2); % 1 is stepsize
stopF = find(HF_Hi-1/2 < HzScale & HzScale < HF_Hi+1/2);
% compute mean
Freq3 = mean(fftFFR(startF:stopF,1));
