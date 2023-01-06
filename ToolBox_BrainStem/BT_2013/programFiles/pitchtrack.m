function [time autocorr LAG FFT freqAxis preFFT blocks]= pitchtrack(avgname, block, step, startAnalysis,channel,exportData)

if block<40
    display('Block Size is too small. Using default: 40 ms');
    block = 40;
end

%Open average file 
[file]= openavg(avgname);
%Define time axis 
timeaxis = linspace(file.xmin, file.xmax, file.pnts)';
    
 % extract signal
 SIGNAL = file.signal(:,channel);
 % get sampling rate
 fs = file.rate;

 
 % ---------------------------PRESTIM, ---------------------------
 % extract portion of prestimulus time, and detrend. The size of the prestim portion is dependent on the block size.
 % To do SNR calculations block and prestim must be same the same number of ms.
 % If the block size is 40ms, then only 40ms of the prestim will be
 % extracted. 
 
 PRESTIM = SIGNAL(1: ms2row(file, block));   %start with the very first point.
 % ramp
 r = hann(size(PRESTIM, 1));  % the entire prestim is ramped.
 % ramp and detrend
 PRESTIM = detrend(PRESTIM.*hann(size(PRESTIM,1)), 'constant');
 % FFT (zero-padded to sampling rate);
 preFFT = abs(fft(PRESTIM, fs));
 %scale preFFT
 preFFT= preFFT*(2./length(PRESTIM));
 preFFT=preFFT(1:1001,1);  %truncate above 1000 Hz
 
   
%  ------------------ FFTS of RESPONSE CHUNKS-----------------------
j = startAnalysis; % each time through loop j increases by step size;

    chunks = 5000;    % an arbitrary maximum number of blocks that the program will create. 
                        
    for k = 1:chunks;   %the program knows to stop once file.xmax is exceeded 

        % variables created: 
        ramptime = (block/1000);   % ramp the entire chunk
        start = j;
        stop = j+block;

        if stop>(file.xmax)  % if stop exceeds the maximum ms time then abort and break out from loop
            k=k-1;
            j=j-step;
            break;
        else
            signal = detrend(SIGNAL(ms2row(file, start):ms2row(file, stop)), 'constant');   % de-mean to zero           
        end
        
        midpoint(k) = mean(ms2row(file, start):ms2row(file,stop));  %calculates the time corresponding to the midpoint of the chunk

        % generate ramp
        ramp = hann(size(signal,1));
        % ramp and de-mean
        signal = detrend(signal.*ramp, 'constant');


        % autocorrelation (see Boersma 1993)
        [c lag]=xcorr(signal, 'coeff');
       
        
        [cwin lagwin]=xcorr(ramp, 'coeff');
        LAG = linspace(-block, block, length(c));
    
        autoc=c./cwin;
        autoc(autoc>1)=1;  %this handles the very rare case that the remainder of the previous step is >1.
        % only plot the first 15 ms; % lowest frequency is ~66 Hz.
        
        
            
            startlag = find(LAG==0);  
            endlag = find(LAG==closestrc(LAG, 15));

        
        % truncate lag and r value matrices to only include values up to first 15 ms.
        autocorr(:,k)=autoc(startlag:endlag)';
        LAG = LAG(startlag:endlag)';
     
        ostartlag = startlag;
        
        ostoplag = endlag;
         
         %%% Now do FFTs;
        % fft, pads to sampling rate
        fftsignal{k} = abs(fft(signal, fs));
        % only go up to 1000 Hz;
        FFT{k} = fftsignal{k}(1:1001,1);
        FFT{k}= FFT{k}*(2./length(signal));
        freqAxis = linspace(0, 1000, 1001);
        
        j = j+step;  % loop through next time chunk
        
       
       
    end
    time = timeaxis(round(midpoint));
    
    blocks = k;
    
    FFT = cell2mat(FFT);
    
    
    if exportData == 1
        [fpath fname ext]=fileparts(avgname);
        FFTfile = [fpath, '\', fname, '-FFTmatrix.xls'];
        ACfile = [fpath, '\', fname, '-ACmatrix.xls'];
        
        xlswrite(FFT, fname, {'FFT matrix'}, FFTfile, 'Sheet1');
        xlswrite(freqAxis', fname, {'Frequency Axis'}, FFTfile, 'Sheet2');
        xlswrite(time, fname, {'Time Axis'}, FFTfile,     'Sheet3');
        
        xlswrite(autocorr, fname, {'autocorrelation matrix'}, ACfile, 'Sheet1');
        xlswrite( LAG, fname, {'Lag Axis'} , ACfile, 'Sheet2');
        xlswrite(time, fname, {'Time Axis'} , ACfile, 'Sheet3');
       
     
    end