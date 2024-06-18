function [SNR_power_Norm_allsubj, aWin,freq_harmonics] = compute_spectral_snr(OPTIONS,flag_sub_to_create, neural_lag)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly do : 
% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects_to_process = subjects(flag_sub_to_create) ; 

FONTSZ = 12 ; 

% Get signal of the stimuli
sti = openavg(fullfile(fileparts(mfilename('fullpath')),'ToolBox_BrainStem','BT_2013','da_170_kraus_16384_LP3000_HP80.avg'));
vTime = table2array(readtable(fullfile(OPTIONS.indir,'ABR_timepoints.txt')));
% resampler à la même freq que le signal ABR ?
% sti.signal = resample(sti.signal,FS,sti.rate); 
   
%% Filtrer 
band_filter = [80,3000]; 
B = fir1(7,[band_filter(1),band_filter(2)]*2/sti.rate,'bandpass'); % Filter created
sti.signalfiltered = filter(B,1,sti.signal); % Stimulus filtered  

tran = [10 55];  % Time window of the stimulus consonant transition(ms)  
cons = [55 170]; % Time windows of the stimulus constant portion(ms)
baseline = [-40 0];  % Prestimulus time windows (ms)

sigSpectralWin = 10; % Frequency windows in which the spectral amplitude of the F0 and HH will be computed. ±5 Hz around each individual peak. (Hz) 
noiseSpectralWin = 15.3; % Frequency windows in which the spectral amplitude of the noise at each side of the sigSpectralWin will be computed. (Hz)

% Defines F0 and harmonics
fund_freq = 100.3;                          % Fundamental frequency (F0)(Hz)
nbHarm = 13;                                % Number of the harmonics (HH) below 1500 that will be analyzed (Hz)
freq_harmonics = fund_freq.*(1:nbHarm+1);   % Row vector with all the harmonics that will be analyzed (Hz)
Ngap = 0 ;                                  % Separation (optional gap) in Hz between signal and noise
max_psd =0;
nplot = [1,3,5];

% Possibility to claissify the SNR values by forman (F0, F1, HH see
% Anderson et al. 2015) 
% WARNING : taking the amplitude values may be tricky because it is not
% normalize and may be impacted by the 1/f effect
% forman_classif_tran = {[100.3];[520 715];[]}; 
% forman_classif_cons = {[100.3];[713];[]}; 
SNR_power_Norm_allsubj = ones(length(subjects_to_process), 3,nbHarm+1);

for ss=1:length(subjects_to_process) %for each subject
    
 
    % Create a folder for files specific to BT_toolbox
    BT_folder = fullfile(OPTIONS.indir, subjects_to_process{ss},'BT_toolbox_formatted');
    fname_avg = fullfile(BT_folder,strcat(subjects_to_process{ss},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.avg')) ;
    
    % Loads FFR rsponse 
    ffr_response = openavg(fname_avg);

    %% NOTE TO OURSELF : WE STOP HERE : the neural lags as comoputed by us (input neural_lag are not the same as the ones obtained with Carles function (below)) 
    % The whole time window
    stim = sti.signalfiltered ; 
    ffr = ffr_response.signal ; 

    % Time window to analyse
    transitionWin = tran + neural_lag(ss) ;    % Time windows of the CV transition (ms)
    steadyWin = cons + neural_lag(ss);       % Time windows of the steady portion (vowel) (ms)
    aWin = {transitionWin; steadyWin; baseline};
    
    % Get number of window
    nbWin = length(aWin);
    
    % Creating the variables of each of the FFR parameters that will be computed.
    rms = ones(nbWin,1); 
    SNR_TD = ones(nbWin,1);
    SNR_power = ones(nbWin,nbHarm+1);
    SNR_power_Norm = ones(nbWin,nbHarm+1);  
    SNR_amp = ones(nbWin,nbHarm+1);

    % Creates SNR summary figure 
    if OPTIONS.display ; h_figsnr = figure('Units','Normalized','Position',[0.3,0,0.5,1],'Name',subjects_to_process{ss}) ; end

    % Loop through the time window to analyse
    for vWin = 1:nbWin 

        % Time to data points
        winSamples = round(((aWin{vWin}+abs(baseline(1)))/1000)*sti.rate)+1;
        time_samples = winSamples(1):winSamples(2);
        winSamples_pres = round(((baseline+abs(baseline(1)))/1000)*sti.rate)+1;
        time_samples_pres = winSamples_pres(1):winSamples_pres(2); % 534 data points
        
        % Root Mean Square (rms)
        rms(vWin) = sqrt(mean(ffr(time_samples).^2));
     
        % Computes pwelch
        X = ffr(time_samples)-mean(ffr(time_samples));
        WINDOW = length(time_samples_pres); 
        NOVERLAP = 542;
        resol = 0.1; % Desired spectral resolution
        NFFT = sti.rate/resol;
        
        % Computes pwelch (power) 
        [pow_spect_density,freqs] = pwelch(X,WINDOW,NOVERLAP,NFFT,sti.rate,'power'); % 497 samples corresponds to 33 ms which is equivalent to 82.5% of 40 ms

        % Get amplitude (square power)
        spectral_total_amp = sqrt(pow_spect_density); 
    
        if OPTIONS.display
                % Display only frequencies <1500 Hz
                vF= freqs<1500;
                
                % % Set current figure to SNR display 
                % set(0,'CurrentFigure',h_figsnr)
                % 
                % Displays power spectal density for this window 
                h_psd(vWin)=subplot(3,2,nplot(vWin),'Parent', h_figsnr) ; 
                plot(freqs(vF), pow_spect_density(vF)) ; 
                title(sprintf('%s, WINDOW --> [%1.1f ; %1.1f] ms',strrep(subjects_to_process{ss},'_','-'),aWin{vWin}(1),aWin{vWin}(2)));
                set(h_psd,'XTick',freq_harmonics,'XTickLabel', string(freq_harmonics),'Fontsize',FONTSZ); grid on ; 
                xlabel('Amplitude uV') ; ylabel('Frequency Hz'); 
        
                % Keep macximum for final display adjustement
                max_psd = max(max(pow_spect_density),max_psd) ; 
            end
    
    
        % For all harmonics (13)
        for vHarm = 1:length(freq_harmonics)
            
            %  SNR (SNR_FD)
            Tw = [freq_harmonics(vHarm)-sigSpectralWin/2,freq_harmonics(vHarm)+sigSpectralWin/2]; % Extremos de la ventana de frecuencias de la señal, localizados a una distancia de 5 Hz arriba y abajo de la frecuencia de interes (f).
            Tw_bool = and(freqs<=Tw(2),freqs>=Tw(1)); % Cada una de las frecuencias que componen la ventana de frecuencias de la señal. 
            Nw_low = [Tw(1) - Ngap - noiseSpectralWin, Tw(1) - Ngap]; % Extremos de la ventana de frecuencias del ruido inferior, localizada por debajo de la ventana de frecuencias de la señal.
            Nw_low_bool = and(freqs<=Nw_low(2),freqs>=Nw_low(1)); % Cada una de las frecuencias que componen la ventana de ruido inferior. 
            Nw_high = [Tw(2) + Ngap, Tw(2) + Ngap + noiseSpectralWin]; % Extremos de la ventana de frecuencias del ruido superior, localizada por encima de la ventana de frecuencias de la señal.
            Nw_high_bool = and(freqs<=Nw_high(2),freqs>=Nw_high(1)); % Cada una de las frecuencias que componen la ventana de ruido superior. 
            Nw_all = or(Nw_low_bool,Nw_high_bool); % Con el "or" incluimos las dos ventanas de ruido
            
            sig_pow_uV2(vWin,vHarm) = mean(pow_spect_density(Tw_bool)); % Extraemos el promedio de power correspondiente a la ventana de frecuencias de la señal.
            noise_pow_uV2 (vWin,vHarm) = mean(pow_spect_density(Nw_all)); % Extraemos el promedio de power correspondiente a las dos ventanas de ruido, a la inferior y a la superior.
            
            SNR_power(vWin,vHarm) = sig_pow_uV2(vWin,vHarm)/noise_pow_uV2(vWin,vHarm); % Dividir la ventana de power correspondiente a la ventana de frecuencias de la señal entre la ventana de power correspondiente a les dos ventanas de ruido.
           
            % Display Noise/SIgnal for 1st harmonic, steadyWin 
            if OPTIONS.display && (vHarm==1) && (vWin==2)
             
                % % Creates figure 
                % h_fig_spect = figure('Units','Normalized','Position',[0.3,0,0.5,1]) ; 
                
                subplot(3,2,2);
                % PLot signal (red) and noise (blue) 
                h_spect = plot(find(Nw_all),pow_spect_density(Nw_all),'b*'); hold on ; 
                plot(find(Tw_bool),pow_spect_density(Tw_bool),'r*') ;
                
                % Arrange axes (display frequencies in X) 
                all_axes= sort(cat(1,find(Nw_all),find(Tw_bool)));
                h_spect.Parent.XTick = fix(linspace(all_axes(1), all_axes(end),30));
                h_spect.Parent.XTickLabel = string(freqs(fix(linspace(all_axes(1), all_axes(end),30))));
                xlabel('Frequency (Hz)'); 
                
                % General documentation (titla grid legend)
                grid on ; legend('noise','signal');
                [~,idx_max] = max(pow_spect_density(Tw_bool|Nw_all)) ; 
                title(sprintf('%s, WINDOW --> [%1.1f ; %1.1f] ms ; Harmonic : %1.1f ; MAX = %1.1f Hz',strrep(subjects_to_process{ss},'_','-'),aWin{vWin}(1),aWin{vWin}(2),freq_harmonics(vHarm),freqs(idx_max+find(Tw_bool|Nw_all,1))));
                
            end
    
            SNR_power_Norm(vWin,vHarm) = 10*log10(SNR_power(vWin,vHarm)); % Passar a dB el SNR_power
           
            mean_amp_uV(vWin,vHarm) = mean((sqrt((pow_spect_density(Tw_bool)))));
            noise_amp_uV(vWin,vHarm) = mean((sqrt((pow_spect_density(Nw_all)))));
            
            SNR_amp(vWin,vHarm) = mean_amp_uV(vWin,vHarm)/noise_amp_uV(vWin,vHarm); % Dividir la ventana de amplitud espectral correspondiente a la ventana de frecuencias de la señal entre la ventana de power correspondiente a las dos ventanas de ruido.

           
            % SNR_amp_Norm(vWin,vHarm) = 20*log10(SNR_amp(vWin,vHarm));
             % Recordar que el SNR obtenido con amplitud y con frecuencias
             % no puede ser el mismo porque al pasar la potencia a amplitud
             % y aplicar el cuadrado, estamos calculando el promedio y
             % despues a ese promedio le aplicamos el cuadrado. Pero no es
             % lo mismo que la suma de elementos elevados al cuadrado
             % (producto notable). Por tanto, decidimos quedarnos con el
             % SNR_power_Norm
        end %end of vWin

    end % end of vHarm
 SNR_power_Norm_allsubj(ss,:,:) = SNR_power_Norm ; 

if OPTIONS.display
    
    set(h_psd,'ylim',[0 max_psd],'Parent', h_figsnr);
    
    % Link the x axis of the two axes together
    linkaxes(h_psd, 'xy')
    
    h_snr = subplot(3,2,[4 6]) ; hp(1) = plot(SNR_power_Norm(1,:),'r+'); hold on ; hp(2) = plot(SNR_power_Norm(2,:),'b*'); hp(3) = plot(SNR_power_Norm(3,:),'mX'); legend('transition','steady','baseline'); grid on ;
    set(hp,'MarkerSize',12,'LineWidth',1.5);
    title(strrep(subjects_to_process{ss},'_','-'));
    set(h_snr,'XTick',1:length(freq_harmonics),'XTickLabel', string(freq_harmonics),'Fontsize',FONTSZ); grid on ; 
    xlabel('Harmonics') ; ylabel('SNR'); ylim([-10 10]);
end

if OPTIONS.savefig == 1 && OPTIONS.display == 1
    % Save figure
    print('-dpng',fullfile(OPTIONS.plot_dir, 'png_folder', strcat('spectral_SNR_',subjects{ss},'_FFR_',OPTIONS.params))); 
end

end


end
