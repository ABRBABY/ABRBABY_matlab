function [out_filenames] = rej_and_prepare_input_brainstem(ALLEEG, OPTIONS,tube_length, propag_sound,flag_sub_to_create, count,suffix, stepA)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%%%%%%%%% TO UPDATE %%%% this comes from old code 

% This function mainly do : 
%Check if stepA(number).set file exists
%Filter data
%Compute mean activity (FFR)
%Add tube delay to mean
%Export FFR data into .txt file
%Convert .txt file into .avg for BT_toolbox
%Export timepoints from last subject

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

%Check if stepA(number).set files exist for all subjects
% for ii=1:length(subjects)
%     setname = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_reref_epoched_FFR_stepA',num2str(stepA),'.set')));
%     if isempty(setname) ; error('_reref_epoched_FFR_stepA%s file does not exist for subject %s', num2str(stepA),subjects{ii}); end
% end

%Loop through subjects
for ii=1:length(subjects) 

     % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    %Set stepA file to work on
    file_stepA = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'))) ;
    
    %Get filepath
    filepath = file_stepA.folder ;
    
    % Creates resulting filename
    out_filenames{ii} = fullfile(indir,subjects{ii}, strcat(subjects{ii},suffix,num2str(count),'.set')) ; 
    
    %Load the stepA .set file to work on
    EEG = pop_loadset(strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'), filepath) ; 
    
    % Reject bad trials and write a report
     
    [EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,1,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1); 
    
    % Name of the file report 
    suf = strsplit(suffix,'_');
    fname = fullfile(file_stepA.folder,strcat(subjects{ii},'_infos_trials','_low_',num2str(rej_low),'_high_',num2str(rej_high),'_',suf(end),num2str(count),'.csv')) ; 
  
    % Write csv file directly into the subject dir
    %produce_report(fname{1}, EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high) ; 
    produce_report(fname{1}, EEG, 1, bloc, win_of_interest, rej_low, rej_high) ; 
    
    %Extract mean activity (erp) and replace data
%     abr = mean(EEG.data(1,:,:),3);   %EEG.data = elec x samples x trials
%     EEG.data = abr;

    %Extract mean activity (ERP) by polarity: averaged, positive only and
    %negative only
   
    % Averaged polarities
    abr_average = mean(EEG.data(1,:,:),3);   %EEG.data = elec x samples x trials     ones_vector = ones(1,size(EEG.data, 3)/2);
    
    % Positive and negative polarities
    zeros_vector = zeros(1,size(EEG.data, 3)/2);
    ones_vector = ones(1,size(EEG.data, 3)/2);
    negative_altern = [zeros_vector(:) ones_vector(:)]';
    negative_altern = negative_altern(:);
    negative_altern_log = logical(negative_altern);

    positive_altern = [ones_vector(:) zeros_vector(:)]';
    positive_altern = positive_altern(:);
    positive_altern_log = logical(positive_altern);

    abr_positive = mean(EEG.data(1,:,positive_altern_log),3);
    abr_negative = mean(EEG.data(1,:,negative_altern_log),3);

    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 
    %abr_shifted = circshift(abr,nsample_delay) ;
    abr_averaged_shifted = circshift(abr_average,nsample_delay) ;
    abr_averaged_shifted(1:nsample_delay) = 0.00001;
    abr_positive_shifted = circshift(abr_positive,nsample_delay) ;
    abr_positive_shifted(1:nsample_delay) = 0.00001;
    abr_negative_shifted = circshift(abr_negative,nsample_delay) ;
    abr_negative_shifted(1:nsample_delay) = 0.00001;

% figure ;
% plot(timepoints,abr_average,'Color',[0 0 0],'Linewidth',0.5); hold on ; plot(timepoints,abr_averaged_shifted,'Color',[1 0 0],'Linewidth',0.5); hold on; set(gca,'YDir','reverse') ;
% grid on ;
% legend('FFR : AVG polarity', 'FFR : Negative polarity, shifted');
% 
% figure ;
% plot(timepoints,abr_positive,'Color',[0 0 0],'Linewidth',0.5); hold on ; plot(timepoints,abr_positive_shifted,'Color',[1 0 0],'Linewidth',0.5); hold on; set(gca,'YDir','reverse') ;
% grid on ;
% legend('FFR : Positive polarity', 'FFR : Positive polarity, shifted');
% 
% figure ;
% plot(timepoints,abr_negative,'Color',[0 0 0],'Linewidth',0.5); hold on ; plot(timepoints,abr_negative_shifted,'Color',[1 0 0],'Linewidth',0.5); hold on; set(gca,'YDir','reverse') ;
% grid on ;
% legend('FFR : Negative polarity', 'FFR : Negative polarity, shifted');
    
    % Create a custom history variable to keep track of OPTIONS 
    EEG.history_stepB = OPTIONS ;
    
    %% SAVE DATASETS
    EEG.data = abr_averaged_shifted ;
    EEG.positive_pol = abr_positive_shifted;
    EEG.negative_pol = abr_negative_shifted;
    pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_stepA',num2str(stepA),suffix,num2str(count))),'gui','off');

    %% Export ABR data into .txt file
%     fname_out = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_shifted_data_HF.txt')) ;
%     fid = fopen(fname_out,'w');
%     fprintf(fid,'%c\n',abr_shifted);
%     fclose(fid);

    fname_out_avg = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_avg_shifted_data_HF.txt')) ;
    fid = fopen(fname_out_avg,'w');
    fprintf(fid,'%c\n',abr_averaged_shifted);
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
function [] = produce_report(fname,EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high)

% Get indices of the trials which were rejected (without messing around with the relative indices)
[~,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);

% Extract variables of interest
trial_index = 1:EEG.trials;
trial_num = [EEG.event.urevent];
condition = {EEG.event.type} ;
latency = [EEG.event.latency]/EEG.srate;
rejected = ismember(trial_index,idx_rejected_all) ;

% Create table to store these information
list_trial_infos = table(trial_index',condition',latency', trial_num',rejected', bloc',...
    'VariableNames', {'trial_index', 'condition', 'latency','trial_num','rejected','bloc'}) ;

%  Save this table into a csv file (use function writetable)
writetable(list_trial_infos,fname, 'WriteVariableNames', true) ;

end
