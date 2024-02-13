function [out_filenames] = rej_and_prepare_input_brainstem(ALLEEG, OPTIONS,tube_length, propag_sound,flag_sub_to_create, count,suffix, stepA)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%%%%%%%%% TO UPDATE %%%% this comes from old code 

% This function mainly do : 
% Reject bad epoch based on window of rejection in OPTIONS
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

% Get OPTIONS
[indir, rej_low, rej_high, BT_toolbox, bloc, win_of_interest]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Inititalize output parameter
out_filenames = [] ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

%Loop through subjects
for ii=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    %Set stepA file to work on
    file_stepA = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'))) ;
     
    % Error if rfe file does not exist
    if isempty(file_stepA) ; error('File stepA%s does not exist for subject %s', num2str(stepA), subjects{ii}); end
 
    %Get filepath
    filepath = file_stepA.folder ;
    
    % Creates resulting filename
    out_filenames{ii} = fullfile(indir,subjects{ii}, strcat(subjects{ii},suffix,num2str(count),'.set')) ; 
    
    %Load the stepA .set file to work on
    EEG = pop_loadset(strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'), filepath) ; 

    % Select only ABR elec
    EEG = pop_select(EEG, 'channel',{'ABR'});

    % Get positive and negative polarities
    zeros_vector = zeros(1,size(EEG.data, 3)/2);
    ones_vector = ones(1,size(EEG.data, 3)/2);
    negative_altern = [zeros_vector(:) ones_vector(:)]';
    negative_altern = negative_altern(:);
    negative_altern_log = logical(negative_altern);

    positive_altern = [ones_vector(:) zeros_vector(:)]';
    positive_altern = positive_altern(:);
    positive_altern_log = logical(positive_altern);
    
    % Read trials_description.txt file
    fname_trial_desc = fullfile(indir,subjects{ii},strcat(subjects{ii},'_trials_description.txt'));
    
    if ~exist(fname_trial_desc,'file')
        error('\nABRBABY --------- File _trials_description.txt does not exist. You must run the first part of the analysis');
    else
        % Read _trials_description file 
        T1 = readtable(fname_trial_desc); 

        %% Reject bad trials

         % Get indices of the trials which were rejected (in the EEG structure referential)
        [EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,EEG.nbchan,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1); 
    
         %% Identifies (flag) the automatically rejected trials

        % Get indices of variable for rejection at acquisition
        idx_rejacq = find(contains(T1.Properties.VariableNames,'rejection_acq')) ; 
        
        % Get indices of HF
        idx_HF = matches(T1.condition,'HF') ;

        % Get indices of trials from EEG structure in the 6000 referential
        idx_EEG_referential = find(idx_HF.*T1{:,idx_rejacq}) ;% Retrieve indices of rejected trials on the 6000 referential
        idx_rejected_6000_ref = idx_EEG_referential(idx_rejected_all) ;

        % Init a vector 6000 trials
        init_rejected = ones(1,height(T1));
    
        % Update flag values 
        init_rejected(idx_rejected_6000_ref) = 0;
    
        % Update trial_description.txt
        header = char(strcat(subjects{ii},'_infos_trials','_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_stepA',num2str(stepA),suffix,num2str(count))) ;        
        add_flag_column_trials_description(fname_trial_desc, header,init_rejected);


    end
    
%     % Name of the file report 
%     suf = strsplit(suffix,'_');
%     fname = fullfile(file_stepA.folder,strcat(subjects{ii},'_infos_trials','_low_',num2str(rej_low),'_high_',num2str(rej_high),'_',suf(end),num2str(count),'.csv')) ; 
%       
%     % Write csv file directly into the subject dir
%     %produce_report(fname{1}, EEG, 1, bloc, win_of_interest, rej_low, rej_high, suf, count) ; 

    % Compute mean FFR without rejected trials (averaged polarities)
    abr_average = mean(EEG_all.data(1,:,:),3);   %EEG.data = elec x samples x trials     ones_vector = ones(1,size(EEG.data, 3)/2);

    % Reject trials from positive and negative vectors
    negative_altern_log_rej = negative_altern_log ;
    negative_altern_log_rej(idx_rejected_all) = 0;

    positive_altern_log_rej = positive_altern_log ;
    positive_altern_log_rej(idx_rejected_all) = 0;

    % Compute mean positive and negative FFR without rejected trials
    abr_positive = mean(EEG.data(1,:,positive_altern_log_rej),3);
    abr_negative = mean(EEG.data(1,:,negative_altern_log_rej),3);

    % Compute mean FFR by subtracting 
    abr_subtracted = (abr_negative - abr_positive) /2 ;
    
    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 
    abr_averaged_add_shifted = circshift(abr_average,nsample_delay) ;
    abr_averaged_add_shifted(1:nsample_delay) = 0.00001;
    abr_averaged_sub_shifted = circshift(abr_subtracted,nsample_delay) ;
    abr_averaged_sub_shifted(1:nsample_delay) = 0.00001;
    abr_positive_shifted = circshift(abr_positive,nsample_delay) ;
    abr_positive_shifted(1:nsample_delay) = 0.00001;
    abr_negative_shifted = circshift(abr_negative,nsample_delay) ;
    abr_negative_shifted(1:nsample_delay) = 0.00001;

    % Create a custom history variable to keep track of OPTIONS 
    EEG.history_stepB = OPTIONS ;
    
    %% SAVE DATASETS
    EEG.data = abr_averaged_add_shifted ;
    EEG.averaged_sub = abr_averaged_sub_shifted ; 
    EEG.positive_pol = abr_positive_shifted;
    EEG.negative_pol = abr_negative_shifted;
    pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_stepA',num2str(stepA),suffix,num2str(count))),'gui','off');

    %% Export ABR data into .txt file

    fname_out_avg = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_avg_shifted_data_HF.txt')) ;
    fid = fopen(fname_out_avg,'w');
    fprintf(fid,'%c\n',abr_averaged_add_shifted);
    fclose(fid);

    fname_out_sub = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_sub_shifted_data_HF.txt')) ;
    fid = fopen(fname_out_sub,'w');
    fprintf(fid,'%c\n',abr_averaged_sub_shifted);
    fclose(fid);

    fname_out_pos = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_pos_shifted_data_HF.txt')) ;
    fid = fopen(fname_out_pos,'w');
    fprintf(fid,'%c\n',abr_positive_shifted);
    fclose(fid);
    
    fname_out_neg = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_neg_shifted_data_HF.txt')) ;
    fid = fopen(fname_out_neg,'w');
    fprintf(fid,'%c\n',abr_negative_shifted);
    fclose(fid);
    addpath(BT_toolbox);

    %addpath(BT_toolbox,'programFiles');
    %bt_txt2avg(fname_out, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
    bt_txt2avg(fname_out_avg, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
    bt_txt2avg(fname_out_sub, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
    bt_txt2avg(fname_out_pos, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
    bt_txt2avg(fname_out_neg, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);

end

%% Export times (from any subject : just timepoints)
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
end


%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir, rej_low, rej_high, BT_toolbox, bloc, win_of_interest]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
rej_low = OPTIONS.rej_low ;                           
rej_high = OPTIONS.rej_high ; 
%varhistory = OPTIONS.varhistory ;
bloc = OPTIONS.bloc ;
BT_toolbox = OPTIONS.bt_toolbox;
win_of_interest = OPTIONS.win_of_interest ;
end

%--------------------------------------------------------------
% FUNCTION that write a report on rejected trials. Regardless to
% conditions) -> we re-excute pop_eegthresh on all trials
% (we do not save .set but the report)
%--------------------------------------------------------------
function [] = produce_report(fname,EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high, suf, count)

% Get indices of the trials which were rejected (without messing around with the relative indices)
[~,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);

% Extract variables of interest
trial_index = 1:EEG.trials;
trial_num = [EEG.event.urevent];
condition = {EEG.event.type} ;
latency = [EEG.event.latency]/EEG.srate;
rejected_all = ismember(trial_index,idx_rejected_all) ;

% Create table to store these information
list_trial_infos = table(trial_index',condition',latency', trial_num',rejected_all', bloc','VariableNames', {'trial_index', 'condition', 'latency','trial_num','rejected','bloc'}) ;

%  Save this table into a csv file (use function writetable)
writetable(list_trial_infos,fname, 'WriteVariableNames', true) ;

end
