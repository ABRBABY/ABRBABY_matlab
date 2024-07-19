function  [PITCH_ERROR_AC,PITCH_ERROR_FFT, PITCH_STRENGTH2, PITCH_SRCORR, vTime, vFreqAC, vFreqFFT,vTime_stim, vFreqAC_stim, vFreqFFT_stim] = compute_pitchtracking(OPTIONS,flag_sub_to_create )

step = OPTIONS.step ; 
block = OPTIONS.blocksz ; 
startSTIM = OPTIONS.startSTIM; 
neural_lag = OPTIONS.expectedNeuralag ; 
stimulus_name = OPTIONS.stim ;
stimulus_name_path = OPTIONS.BT_toolbox ; 
endSTIM = OPTIONS.endSTIM; 
minFrequency_stim = OPTIONS.minFrequency_stim; 
maxFrequency_stim = OPTIONS.maxFrequency_stim;
minFrequency = OPTIONS.minFrequencyR ; 
maxFrequency = OPTIONS.maxFrequencyR ; 

startRESP = startSTIM + neural_lag;


PITCH_ERROR_AC = [];
PITCH_ERROR_FFT = [];
PITCH_STRENGTH2 = []; 
PITCH_SRCORR = [];
vTime = []; 
vFreqAC = [];
vFreqFFT = []; 
vTime_stim = []; 
vFreqAC_stim = []; 
vFreqFFT_stim = []; 

% This function mainly do : 
% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects_to_process = subjects(flag_sub_to_create) ; 

for ss=1:length(subjects_to_process) %for each subject
    
    response_name_path = fullfile(OPTIONS.indir,subjects_to_process{ss},'BT_toolbox_formatted/'); 
    tmp = dir(fullfile(response_name_path,'*avg_*.avg'));
    
    % pitchtrack response
    [time autocorr lag FFT_resp freqaxis prestimFFT totalblocks]= pitchtrack(fullfile(response_name_path,tmp.name), block, step, startRESP, 1, 0);
    
    % pitchtrack stimulus (make conditional)
    [time_stim autocorr_stim lag_stim FFT_stim null freqaxis_stim totalblocks_stim]=pitchtrack(fullfile(stimulus_name_path ,stimulus_name), block, step, startSTIM, 1,0);
    
    time_stim = time_stim + neural_lag;  %for plotting purpose shift stimulus forward in time.
    
    [v stopPT]=closestrc(time_stim, endSTIM+(block./2));
    
    autocorr_stim(:, stopPT+1:end)=[];
    time_stim(stopPT+1:end)=[];
    FFT_stim(:, stopPT+1:end)=[];
    totalblocks_stim = stopPT;
    
    %Extract F0 from stimulus using Autocorrelation Method 
    freqAC_vector = 1000./lag_stim;
    [s LagStart_stim]=closestrc(freqAC_vector,maxFrequency_stim);  %s is dummy variable
    [s LagStop_stim] =closestrc(freqAC_vector,minFrequency_stim); 
    [R_stim index] = max(autocorr_stim(LagStart_stim:LagStop_stim, :));
    FreqAC_stim=freqAC_vector(LagStart_stim+index-1);
    clear index s
    
    %Extract F0 from response using Autocorrelation Method 
    freqAC_vector = 1000./lag;
    [s LagStart]=closestrc(freqAC_vector,maxFrequency);  %s is dummy variable
    [s LagStop] =closestrc(freqAC_vector,minFrequency); 
    [R index] = max(autocorr(LagStart:LagStop, :));
    FreqAC=freqAC_vector(LagStart+index-1);
    LAG = 1000./FreqAC;
    
    %Extract F0 from stimulus using FFT Method 
    [FFTAMP_stim index]= max(FFT_stim(minFrequency_stim+1:maxFrequency_stim+1, :));
    FreqFFT_stim = minFrequency_stim+index-1;
    FreqFFT_stim = FreqFFT_stim';
    
    %Extract F0 from response using FFT Method 
    [MaxFFTAMP index]= max(FFT_resp(minFrequency+1:maxFrequency+1, :));
    FreqFFT = (minFrequency+index-1)'; %flip direction so that it matches FreqAC
    
    %% Determine whether each extracted frequency was above the noise floor (NF)
    % When the pitch-track is plotted a small gray dot will appear over the time
    % ranges where the extracted frequency is below the noise floor The dot will be located above the plot and the location 
    % is based on what was inputted for maxFrequency (maxFrequency+4).  The total (i.e. total_belowNF) that gets exported is based 
    % only on the total number of blocks in the stimulus and not the total
    % number of blocks in the response.  In other words, it only includes the blocks that were used to calculate pitch error.
    
    for x = 1:totalblocks;
        F0_AMPusingACfreqs(x,1) = FFT_resp(round(FreqAC(x))-1,x);  %Find amplitude using frequencies extracted using AC method
    end
    
    mean_F0_AMP = mean(F0_AMPusingACfreqs(1:totalblocks_stim,1));
    F0_AMP_prestim = prestimFFT(round(FreqAC)-1);
    PITCH_SNR = F0_AMPusingACfreqs./F0_AMP_prestim ;
    
    %Determine whether the extracted frequency was also the spectral maximum
    for x = 1:totalblocks;
        peaks{x}=localmax(FFT_resp(:,x));  %Find amplitude using frequencies extracted using AC method
        [closest_spectralpeak(x) row(x) column(x)]=closestrc(peaks{x}, FreqAC(x));
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
    
    PITCH_ERROR_AC(ss) = mean(abs(FreqAC(1:totalblocks_stim)-FreqAC_stim));  %Measured in Hz.
    PITCH_ERROR_FFT(ss) = mean(abs(FreqFFT(1:totalblocks_stim)-FreqFFT_stim));  %Measured in Hz.

    %on the off chance that one or more of the Rs is exactly 1, we must set
    %these Rs to 0.999999 to get a valid number for fisher (i.e. not inf)If you are are reading this you are as annoyed as we are.
    PITCH_STRENGTH = mean(R(1:totalblocks_stim));
    Rtemp=R;
    Rtemp(Rtemp==1)=0.999999;
    PITCH_STRENGTH2(ss) = fisherinv(mean(fisher(Rtemp(1:totalblocks_stim))));
    CORR = corrcoef(FreqAC(1:totalblocks_stim), FreqAC_stim);  %correlation between stimulus and response f0 contour
    PITCH_SRCORR(ss) =CORR(1,2); %the first number is always 1, need to take second 
    total_belowNF =  sum(PITCH_SNR(1:totalblocks_stim)<1);
    total_notatSpectralMax = sum(notspectralMax(1:totalblocks_stim));

    vTime(ss,:) = time; 
    vTime_stim(ss,:) = time_stim; 
    vFreqAC(ss,:) = FreqAC; 
    vFreqFFT(ss,:) = FreqFFT; 
    vFreqAC_stim(ss,:) = FreqAC_stim; 
    vFreqFFT_stim(ss,:) = FreqFFT_stim; 

end
