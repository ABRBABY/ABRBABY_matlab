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
     
      % Extract epochs for HF
%     EEG = pop_epoch( EEG, {  'HF'  }, win_of_interest, 'newname', 'epochs', 'epochinfo', 'yes');
%     
%     %Remove baseline
%     EEG = pop_rmbase( EEG, baseline,[]);

    % Get indices of the trials which were rejected (without messing around with the relative indices)
    %[EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1); 
    [EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,1,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1); 
    
    % Name of the file report 
    suf = strsplit(suffix,'_');
    fname = fullfile(file_stepA.folder,strcat(subjects{ii},'_infos_trials','_low_',num2str(rej_low),'_high_',num2str(rej_high),'_',suf(end),num2str(count),'.csv')) ; 
  
    % Write csv file directly into the subject dir
    %produce_report(fname{1}, EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high) ; 
    produce_report(fname{1}, EEG, 1, bloc, win_of_interest, rej_low, rej_high) ; 
    
    %Extract mean activity (erp) and replace data
    abr = mean(EEG.data(1,:,:),3);
    EEG.data = abr;

    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    abr_shifted = circshift(abr,nsample_delay) ;
    
    % Create a custom history variable to keep track of OPTIONS 
    EEG.history_stepB = OPTIONS ;
    
    %% SAVE DATASET
    pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_stepA',num2str(stepA),suffix,num2str(count))),'gui','off');

    %% Export ABR data into .txt file
    fname_out = fullfile(filepath,strcat(subjects{ii},'_stepA', num2str(stepA),suffix, num2str(count),'_abr_shifted_data_HF.txt')) ;
    fid = fopen(fname_out,'w');
    fprintf(fid,'%c\n',abr_shifted);
    fclose(fid);
    
    addpath(BT_toolbox);
    %addpath(BT_toolbox,'programFiles');
    bt_txt2avg(fname_out, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
end

%% Export times (from any subject : just timepoints)
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
end


%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = read_custom_events(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 0);

% Insert info from .txt into EEG.event
my_events = table2array(my_events);

out_event = struct('latency', {in_event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {in_event(:).urevent});

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
% FUNCTION that check if that sets of param exist 
%--------------------------------------------------------------
function [does_exist, count] = check_exist_set_params(filename, subject, OPTIONS)

% Reads all folders that are in indir 
d = dir(fullfile(OPTIONS.indir,subject, strcat(subject,'stepA',num2str(OPTIONS.stepA),'_filtered_FFR_F*')));
% No file exists 
if isempty(d) ; does_exist=0 ; count =1 ; return ; end

for ff=1:length(d) 
    
    EEG = pop_loadset('filepath',fullfile(OPTIONS.indir,subject, d(ff).name),'loadmode','info') ;  
    
    % Check if the set of param correspond to the current file
    if isequal(EEG.history_f,OPTIONS)
        does_exist = 1 ; 
        tmp = regexp(d(ff).name,'F\d*','Match');
        tmp2 =  regexp(tmp,'\d*','Match');
        count = cell2mat(tmp2{:}); 
        return ; 
    end
    
end

% At this point the set of params does not exist and a new file needs to be
% crated (with a count increment)
does_exist = 0 ; 
count = length(d) +1;

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
