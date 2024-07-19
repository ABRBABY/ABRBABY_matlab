function prepare_input_brainstem(ALLEEG, OPTIONS,tube_length, propag_sound,flag_sub_to_create, count,suffix, RFE)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly do : 
% Extract positive and negative FFR separately
% Compute mean activity (FFR) : 
% - mean added FFR ((positive+negative) /2)
% - mean subtracted FFR ((negative-positive) /2)
% - mean positive only
% - mean negative only
% Add tube delay to means
% Export FFR data into .txt file
% Convert .txt file into .avg for BT_toolbox
% Export timepoints from last subject

% Exist if there is no subjects to compute
if sum(flag_sub_to_create)==0
    fprintf('Nothing to compute \n');
    return
end

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

suffix_stepA = strrep(RFE,'_','') ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

%Loop through subjects
for ii=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Input filename
    fname = strcat(subjects{ii},'_FFR_',suffix_stepA,suffix,num2str(count),'.set'); %dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'))) ;
    
     % Error if rfe file does not exist
    if ~exist(fullfile(OPTIONS.indir, subjects{ii},fname),'file') ; error('File %s does not exist for subject %s', fname, subjects{ii}); end
  
    %Load the stepA .set file to work on
    EEG = pop_loadset(fname, fullfile(OPTIONS.indir, subjects{ii})) ; 
    
    % Select only ABR elec
    EEG = pop_select(EEG, 'channel',{'ABR'});
    
    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    % Get and squeeze data 
    abr_trials = squeeze(EEG.data) ; 

    % Shift
    abr_trials = circshift(abr_trials,-nsample_delay,1);

    % Computes the average of ABR signal 
    abr_average = mean(abr_trials,2) ; 
    
    % Creates a vector of alternating flags
    idx_flip1 = ones(1,size(EEG.data,3));
    idx_flip1(2:2:size(EEG.data,3))=0; 
    abr_average1 = mean(abr_trials(:,idx_flip1==1),2);
    df1 = sum(idx_flip1);

    idx_flip2 = ones(1,size(EEG.data,3));
    idx_flip2(1:2:size(EEG.data,3))=0; 
    abr_average2 = mean(abr_trials(:,idx_flip2==1),2);
    df2 = sum(idx_flip2);

    % Compute mean FFR by subtracting 
    abr_subtracted = (abr_average2 - abr_average1) /2 ;
    
    abr_types = {'avg','sub','neg','pos'};
    abr = cat(2,abr_average,abr_subtracted,abr_average2,abr_average1)' ; 

    if OPTIONS.display==1
        % Display FFR sanity check 
        fig = display_temporal_FFR(subjects{ii},EEG.times/1000,abr_average, OPTIONS.abr_disp_scale);

         % Save figures 
        if OPTIONS.savefigs == 1 
            % Save data in vectoriel in subject folder
            print('-dsvg',fullfile(OPTIONS.svg_folder,strcat(subjects{ii},'_FFR_sanity.svg')));
        
            % Save data in png (with same filename as vectoriel) but different directory
            print('-dpng',fullfile(OPTIONS.png_folder,strcat(subjects{ii},'_FFR_sanity.png')));
        
            % Save data in fig (with same filename as vectoriel) but different directory
            saveas(fig, fullfile(OPTIONS.fig_folder,strcat(subjects{ii},'_FFR_sanity.fig')));
        end       

    end

    % Create a folder for files specific to BT_toolbox
    BT_folder = fullfile(OPTIONS.indir, subjects{ii},'BT_toolbox_formatted');
    if ~exist(BT_folder,'dir') ; mkdir(BT_folder);end

    for ff=1:length(abr_types) 
        fname_out = fullfile(BT_folder,strcat(subjects{ii},'_',num2str(suffix_stepA),suffix, num2str(count),'_abr_',abr_types{ff},'_shifted_data_HF.txt')) ;
        fid = fopen(fname_out,'w');
        fprintf(fid,'%c\n',abr(ff,:));
        fclose(fid);
       
        % Converts the output file into BT_Toolbpx compatible data 
        bt_txt2avg(fname_out, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
  
    end
    % % SAVE DATASETS
    % EEG.data = abr_averaged_add_shifted ;
    % EEG.averaged_sub = abr_average ; 
    % EEG.positive_pol = abr_average1;
    % EEG.negative_pol = abr_average2;
    % pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_stepA',num2str(stepA),suffix,num2str(count))),'gui','off');

end

%% Export times (from any subject : just timepoints)
fname_out = fullfile(OPTIONS.indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
end


%--------------------------------------------------------------
% FUNCTION that displays and save FFR (temporal and spectral)
%--------------------------------------------------------------
function [fig] = display_temporal_FFR(subject, vTimes, abr, abr_scale)

%% Temporal display 

% Plot timeseries
fig = figure('units','normalized','position',[0,0,1,1]); 

subplot(3,2, [1 2]); plot(vTimes,abr,'Color', 'r', 'Linewidth',0.5); hold on; grid on; 

%Add legend
legend(strcat(subject, '_avg'), 'Interpreter', 'None');

% Adjust scales (y-axis and x-axis) (transform in milliseconds)
xlim([vTimes(1), vTimes(end)]); ylim(abr_scale) ; grid on ;

%Display labels
xlabel('Times (ms)'); ylabel('uV'); 

%Add title
title(strcat('FFR Avg'));


%% Sprectral display 
% Adaptation of Skoe function (bt_fftsc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sampling rate
Fs = 1 / (vTimes(2) - vTimes(1)) ; 
FFR = abr;

%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD 
numPoints = length(FFR);

%**** STEP 2. FFT OF FFR
%******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
rampMS = 4/1000; % length of ramp (on and off) in seconds
% hanPoints = 26;  %hard coded to be the same as Biologic's settings (December 7, 2005);
hanPoints = rampMS.*Fs; % length of ramp in points
hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
hanHalfPoints = round(hanPoints./2);
numberOfOnes = numPoints - hanPoints;
FFRhan = hann(hanPoints);  
FFRhan = [FFRhan(1:hanHalfPoints); ones(numberOfOnes,1); FFRhan(hanHalfPoints+1:hanPoints)];

% baseline, window, then baseline again
FFR = detrend(FFR, 'constant');
FFR = detrend(FFR.*FFRhan, 'constant');

%******** STEP 2b. Perform FFT
fftFFR = abs(fft(FFR, round(Fs)));
fftFFR = fftFFR(1:round(round(Fs)/2));
fftFFR = fftFFR.*(2./numPoints); % scale to peak µV
HzScale = [0:1:round(Fs/2)]'; % frequency 'axis'
HzScale = HzScale(1:length(fftFFR));

%% Plot FFT in specific frequency windows
subplot(3,2,3); plot(HzScale,fftFFR); grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", strrep(subject,'_','-')]);
xlabel("Frequency (Hz)");ylabel("Amplitude (µV)");

subplot(3,2,4); plot(HzScale,fftFFR); xlim([90 110]); grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", strrep(subject,'_','-')]);
xlabel("Frequency (Hz)");ylabel("Amplitude (µV)");

subplot(3,2,5); plot(HzScale,fftFFR); xlim([300 500]); grid on; 
title(["Single-Sided Amplitude Spectrum of X(t)",strrep(subject,'_','-')]);
xlabel("Frequency (Hz)"); ylabel("Amplitude (µV)");

subplot(3,2,6); plot(HzScale,fftFFR); xlim([1100 1300]); grid on;
title(["Single-Sided Amplitude Spectrum of X(t)", strrep(subject,'_','-')]);
xlabel("Frequency (Hz)"); ylabel("Amplitude (µV)");

end
