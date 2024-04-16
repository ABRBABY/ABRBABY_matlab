function [lmax] = compute_spectral_snr(OPTIONS,flag_sub_to_create, neural_lag)
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

spectral_snr = []; 

% Get signal of the stimuli
sti = openavg(fullfile(fileparts(mfilename('fullpath')),'ToolBox_BrainStem','BT_2013','da_170_kraus_16384_LP3000_HP80.avg'));
vTime = table2array(readtable(fullfile(OPTIONS.indir,'ABR_timepoints.txt')));
% resampler à la même freq que le signal ABR ?
% sti.signal = resample(sti.signal,FS,sti.rate); 
   
%% Filtrer 
band_filter = [80,1500]; 
B = fir1(7,[band_filter(1),band_filter(2)]*2/sti.rate,'bandpass'); % Filter created
sti.signalfiltered = filter(B,1,sti.signal); % Stimulus filtered  

tran = [10 55];  % Time window of the stimulus consonant transition(ms)  
cons = [55 170]; % Time windows of the stimulus constant portion(ms)
tota = [10 170]; % Time winodws of the stimulus entire duration (ms)
pres = [-40 0];  % Prestimulus time windows (ms)

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

    [r,l,rmax,lmax(ss)] = nbffr_cc(stim,ffr,ffr_response.rate, [0 sti.xmax],[3 13],vTime(1));
    lmax(ss) = lmax(ss)+3;

    % 
    % % CV transition 
    % [r,l,rmax,lmax] = nbffr_cc(sti.signalfiltered,ffr_response.signal,ffr_response.rate,[10, 55],[3, 10],offset);
    % 
    % % Vowel 
    % [r,l,rmax,lmax] = nbffr_cc(sti.signalfiltered,ffr_response.signal,ffr_response.rate,[55,170],[3, 10],offset);
   
end

end

%--------------------------------------------------------------------------------------------------------
% Function from Carles Escera (extracted from
% NBFFR_Analysis_Welch_Amplitude.m)
% stim : stim signal 
% ffr : EEG FFR response 
% FS : sampling rate (they are both the same 16384Hz)
% StimWin : Stimulus time windows in which the cross-correlation with the response will be computed (Carles default = [0 169.80]
% lags_win :  % Time windows in which the neural lag will be computed. BE CAREFULL! The last value should be the value + the first value;
%                        % a limit of 10 ms we should write 13 (10 + 3).
% OUTPUT 
% r : the r value shows the pearson's correlation between the stimulus and the response for each of the delay applied to the response. 
% rmax : rmax is the max. correlation value and index is the data point number that contains rmax.
% lmax : lmax is the delay expressed in ms that corresponds to the max. correlation value (rmax).
%--------------------------------------------------------------------------------------------------------
function [r,l,rmax,lmax] = nbffr_cc(stim,ffr,FS,StimWin,lags_win,offset)

E = stim-mean(stim); % Stimulus demeaned
E = E(round(StimWin(1)*FS/1000)+1:round(StimWin(2)*FS/1000)+1); % % Fragmento del estímulo seleccionado expresado en data points
L = length(E); % % Longitud del fragmento seleccionado del stim expresado en data points
R = ffr-mean(ffr); % Registro demeaned
lags_win = round(FS*(lags_win+abs(offset))/1000)+1; % Data points equivalente a los ms. que corresponderían al primer delay que se retrasa el estimulo(derivado de sumar al valor que determina el inicio del fragmento del estimulo seleccionado, el valor minimo de delay considerado. Ej: inicio fragmento estimulo: 57 ms + valor minimo de delay: 3 ms = 60 ms) y el último permitido (derivado de sumar al valor que determina el inicio del fragmento del estimulo seleccionado, el valor maximo de delay considerado. Ej: inicio fragmento estimulo: 57 ms + valor minimo de delay: 10 ms = 67 ms) sumandoles a ambos valores el offset(para que los data points seleccionados sean sumados al data point correspondiente a 0 ms).
lags_win = lags_win(1):lags_win(2); % Todos los data points existentes comprendidos entre el minimo delay permitido y el máximo.

r = zeros(1,length(lags_win));

for i = 1:length(lags_win)  % Vector que contiene cada uno de los retrasos aplicados a la respuesta
    c1 = sum(E.*E);         % Sumatorio de todos los data points elevados al cuadrado incluidos en la region del estimulo seleccionado.
    c2 = sum(R(lags_win(i):lags_win(i)+L-1).*R(lags_win(i):lags_win(i)+L-1)); % Sumatorio de todos los data points elevados al cuadrado incluidos en el fragmento de respuesta por cada uno de los delays aplicado a la respuesta.
    cc = sum(E.*R(lags_win(i):lags_win(i)+L-1)); % Sumatorio del producto de cada uno de los data points del estimulo por cada uno de los data points de la respuesta por cada uno de los delays aplicado a la respuesta.
    if and(c1~=0,c2~=0)
        r(i) = cc/sqrt(c1*c2); % cada uno de los valores del vector r equivaldrá al valor de correlacion entre estimulo y respuesta obtenido por cada delay aplicado a la respuesta.
    else
        fprintf('toto');
        r(i) = -999;
    end
end
l = 1000*([0:length(lags_win)-1])/FS; % longitud de la ventana temporal expresada en data points.

% figure ; plot(r)

[rmax,index] = max(abs(r)); % El primer output indica el valor de máxima correlacion y el segundo su localización, es decir, el data point que contiene el valor de máxima correlación.
lmax = l(index); % Extraemos el delay expresado en ms que corresponde a la máxima correlacion encontrada entre estimulo y respuesta.
end