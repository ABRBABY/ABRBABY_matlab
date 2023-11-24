function [snr] = FFR_analysis_get_SNR_freq(subjects, OPTIONS)
% ERPs analysis script -
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Get timepoints
FFR_times_file = fullfile(OPTIONS.indir,'ABR_timepoints.txt') ;
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
for loopnum = 1:length(subjects) %for each subject
    FFR_file = fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.txt')) ;
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR preprocessing steps on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r'); FFR_subj = fscanf(fileID,'%f');
    
    % Compute FFT
    [vHertz, fft] = compute_fft(timepoints/1000, FFR_subj); 

    % Compute SNR (Ribas-Prats et al. 2021)
    snr(loopnum) = 10*log(mean(fft(OPTIONS.winSignal))/mean(fft(OPTIONS.winNoise)));

end

end

%% ===== GET DEFAULT MEASURE =====
% USAGE:   [vHertz, fft] = compute_fft(timpoints, ffr)
% Inspired by bt_fftsc developed by E.E. Skoe.
function [vHertz, fftFFR] = compute_fft(timepoints, FFR)
 
    FS = 1 / (timepoints(2)-timepoints(1));
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
    fftFFR = abs(fft(FFR, round(FS)));
    fftFFR = fftFFR(1:round(round(FS)/2));
    fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
    vHertz = [0:1:round(FS/2)]'; % frequency 'axis'
    vHertz = vHertz(1:length(fftFFR));
end

