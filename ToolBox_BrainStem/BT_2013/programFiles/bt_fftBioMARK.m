function [F0 F1 HF fftFFR HzScale]=bt_fftbiomark(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi,xl, chan)
% bt_fftbiomap computes frequency-domain amplitudes of F0, F1 and high-freqency bins of BioMAP response.
%
% Usage: [F0 F1 HF] = bt_fftbiomap('filename.avg',start_latency,stop_latency);
%    Note, BioMAP default latencies are 11.38, 40.58.
% 
% Three variables are returned the workspace:
% (1) F0: mean amplitude over F0_Lo-F0_Hi Hz
% (2) F1: mean amplitude over F1_Lo-F1_Hi Hz
% (3) HF: mean amplitude over HF_Lo to HF_Hi Hz 

% % Dependancies: ms2row, openavg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    disp('Using default FFR period of 11.38-40.58 ms.');
    start = 11.38;
    stop =  40.58;
end

if nargin<10
	xl = 0;
end

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
hanPoints = 26;  %hard coded to be the same as Biologic's settings (December 7, 2005);
hanHalfPoints = round(hanPoints./2);
numberOfOnes = numPoints - hanPoints;
FFRhan = hann(hanPoints);  
FFRhan = [FFRhan(1:hanHalfPoints); ones(numberOfOnes,1); FFRhan(hanHalfPoints+1:hanPoints)];

% baseline, window, then baseline again
FFR = detrend(FFR, 'constant');
FFR = detrend(FFR.*FFRhan, 'constant');

%******** STEP 2b. Perform FFT
Pad = 4096;
fftFFR = abs(fft(FFR, Pad));
fftFFR = fftFFR(1:round(4096/2));
StepSize = FS/4096;
HzScale = [0:StepSize:round(FS/2)-1]'; % frequency 'axis'

% clear variables no longer needed
clear FFRramp rampMS hanPoints hanHalfPoints numberOfOnes FFRhan

%**** STEP 3. compute mean magnitudes over F0, F1 and HF ranges
% i. F0: F0_Lo-F0_Hi.
% find freqs nearest F0_Lo and F0_Hi.
startF = find(F0_Lo-StepSize/2 < HzScale & HzScale < F0_Lo+StepSize/2);
stopF = find(F0_Hi-StepSize/2 < HzScale & HzScale < F0_Hi+StepSize/2);
% compute mean
F0 = mean(fftFFR(startF:stopF,1));

% ii. F1: F1_Lo-F1_Hi
% find freqs nearest F1_Lo and F1_Hi.
startF = find(F1_Lo-StepSize/2 < HzScale & HzScale < F1_Lo+StepSize/2);
stopF = find(F1_Hi-StepSize/2 < HzScale & HzScale < F1_Hi+StepSize/2);
% compute mean
F1 = mean(fftFFR(startF:stopF,1));

% iii. HF: HF_Lo-HF_Hi
% find freqs nearest HF_Lo and HF_Hi.
startF = find(HF_Lo-StepSize/2 < HzScale & HzScale < HF_Lo+StepSize/2);
stopF = find(HF_Hi-StepSize/2 < HzScale & HzScale < HF_Hi+StepSize/2);
% compute mean
HF = mean(fftFFR(startF:stopF,1));
