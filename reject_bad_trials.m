function [] = reject_bad_trials(ALLEEG, OPTIONS, opt_balance, flag_sub_to_create, count, suffix, RFE)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project
%INPUTS:
%- ALLEEG = EEGLAB whole structure
%- OPTIONS = rbt options
%- opt_balance = 'balanced' for balancing number of standard with number of
%deviants, or 'unbalanced' to keep all standards that are not rejected
%after trial rejection

%Get options
[indir, rej_low, rej_high, bloc]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

suffix_stepA = strsplit(RFE,'_') ; 

% Loop though subjects
for ii=1:length(subjects)
    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Set rfe file to work on
    file_stepA = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_',OPTIONS.analysis,RFE,'.set'))) ;
    
    % Error if rfe file does not exist
    if isempty(file_stepA) ; error('File %s does not exist for subject %s', RFE, subjects{ii}); end
    
    %Get subDir
    subDir = file_stepA.folder ;
  
     %Load the RFE .set file to work on
    EEG = pop_loadset(file_stepA.name,subDir) ;
    
    %Get eeg_elec and win_of_interest from RFE set of parameters
    eeg_elec = EEG.history_stepA.eeg_elec ;
    win_of_interest = EEG.history_stepA.win_of_interest ;

    % Read trial_description.txt
    fname_trial_desc = fullfile(OPTIONS.indir,subjects{ii},strcat(subjects{ii},'_trials_description.txt'));
    
    if ~exist(fname_trial_desc,'file')
        error('\nABRBABY --------- File _trials_description.txt does not exist. You must run the first part of the analysis');
    else
     
        % Read _trials_description file 
        T1 = readtable(fname_trial_desc); 

         %% Identifies (flag) the first 3 events in blocks 
        begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
        
        % Init a vector 6000 trials
        init_beg_bloc = ones(1,height(T1));
        
        % Get indices of STD, DEV1, DEV2 
        idx_not_HF =  ~matches(T1.condition,'HF') ; 
       
        % Update flag values 
        init_beg_bloc((T1{:,3}~=0)&idx_not_HF) = ~ismember(find(T1{idx_not_HF,3}~=0),begining_of_block);
        add_flag_column_trials_description(fname_trial_desc, 'begining_block',init_beg_bloc);

        %% Identifies (flag) the automatically rejected trials
        idx_rejacq = find(contains(T1.Properties.VariableNames,'rejection_acq')) ; 
        
         % Get indices of the trials which were rejected (without messing around with the relative indices)
        [~,idx_all_rejected] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
        % Init a vector 6000 trials
        init_rejected = ones(1,height(T1));
    
        % Update flag values 
        init_rejected((T1{:,idx_rejacq}~=0)&idx_not_HF) = ~ismember(find(T1{idx_not_HF,idx_rejacq}~=0),idx_all_rejected);
    
        % Update trial_description.txt
        header = char(strcat(subjects{ii},'_infos_trials','_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_',suffix_stepA(end),suffix,num2str(count))) ;        
        add_flag_column_trials_description(fname_trial_desc, header,init_rejected);

        idx_rejauto = find(contains(T1.Properties.VariableNames,header));

        idx_cond = find(contains(T1.Properties.VariableNames,'condition'));

        
        % pop_selectevent(EEG, idx_to_keep) ;


    % Select trials per conditions
    [EEG_DEV1,~] = pop_selectevent(EEG,'type','DEV1');
    [EEG_DEV2,~] = pop_selectevent(EEG,'type','DEV2');
    [EEG_STD,~] = pop_selectevent(EEG,'type','STD');

    [EEG_STD1_thresh,~] = pop_eegthresh(EEG_STD,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV1_thresh,~] = pop_eegthresh(EEG_DEV1,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV2_thresh,~] = pop_eegthresh(EEG_DEV2,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
   
    
    %% Create a custom history variable to keep track of OPTIONS in each
    % .set saved
    EEG_DEV1_thresh.history_stepB = OPTIONS ;
    EEG_DEV2_thresh.history_stepB = OPTIONS ;
    EEG_STD1_thresh.history_stepB = OPTIONS ;
    % EEG_STD2_thresh.history_stepB = OPTIONS ;

    
    DEV1_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV1_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
    DEV2_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV2_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
    STDD_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_STDD_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
    
    % HERE REMOVE BEGINGIN OF BLOCK 

    % Save datasets 
    pop_newset(ALLEEG, EEG_DEV1_thresh, 1, 'setname',DEV1_fname,'savenew', fullfile(subDir,DEV1_fname),'gui','off');
    pop_newset(ALLEEG, EEG_DEV2_thresh, 1, 'setname',DEV2_fname,'savenew', fullfile(subDir, DEV2_fname),'gui','off');
    pop_newset(ALLEEG, EEG_STD1_thresh, 1, 'setname',STDD_fname,'savenew', fullfile(subDir, STDD_fname),'gui','off');

    % Name of the file report 
    fnameReport = fullfile(subDir,strcat(subjects{ii},'_infos_trials','_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_',suffix_stepA(end),suffix,num2str(count),'.csv')) ; 
    fnameTrialDescription = fullfile(subDir,strcat(subjects{ii},'_trials_description.txt'));
    
   
    % Write csv file directly into the subject dir
    % produce_report(fnameReport{1},fnameTrialDescription, EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high, begining_of_block,opt_balance) ; 

    end 
end
   

end

%--------------------------------------------------------------
% FUNCTION that write a report on rejected trials. Regardless to
% conditions) -> we re-excute pop_eegthresh on all trials 
% (we do not save .set but the report)
%--------------------------------------------------------------
function [] = produce_report(fname,fnameTrials, EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high,begining_of_block, opt_balance) 

    NB_TRIALS = 6000 ;

    if strcmp(opt_balance,'balanced')
        warndlg('The report of trial status (rejected/not rejected) is under developement for ''balanced'' option ','Warning')
        return
    end
     % Get indices of the trials which were rejected (without messing around with the relative indices)
    [~,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);

    % Extract variables of interest
    trial_index = 1:EEG.trials;
    trial_num = [EEG.event.urevent];
    condition = {EEG.event.type} ;
    latency = [EEG.event.latency]/EEG.srate;  
    rejected = ismember(trial_index,unique([idx_rejected_all,begining_of_block])) ; 
    
    % Create table to store these information
    list_trial_infos = table(trial_index',condition',latency', trial_num',rejected', bloc',...
        'VariableNames', {'trial_index', 'condition', 'latency','trial_num','rejected','bloc'}) ;

    % Save this table into a csv file (use function writetable)
    writetable(list_trial_infos,fname, 'WriteVariableNames', true) ; 

    % Update trial_despection
    [~,header,~]=fileparts(fname);
 
    % % Create a flag vector for bad trials
    init = zeros(1,height(T1));
    init(T1{idx_not_HF,3}~=0) = ~rejected ;

    % Edit _trial_description.txt file with new column
    add_flag_column_trials_description(fnameTrials, header,init);

end

%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = update_trial_description(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 1);

if size(my_events,2)~=3 
    error('Wrong number of column in file _trial_descriptions.txt');
else
    idx_events = my_events{:,2}==1 ; 
    out_event = struct('latency', num2cell(my_events{idx_events,3}'), ...
                    'type', my_events{idx_events,1}',...
                    'urevent', num2cell(1:height(my_events))) ;                      
end

end

%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir, rej_low, rej_high, bloc]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
rej_low = OPTIONS.rej_low ;
rej_high = OPTIONS.rej_high ;
bloc = OPTIONS.bloc ; 

end

%--------------------------------------------------------------
% FUNCTION that select from EEG_STD ntrial with no repetition with exisitng
% STD 
%--------------------------------------------------------------
function [EEG_STD_ALL] = balance_number_of_STD(EEG,ntrials,std_good,target_indices,begining_of_block,target_indices_std,idx_std_stepB, idx_std)
   
        % Pool of STD without those rejected by threshold detection
        pool_std = setdiff(std_good,target_indices);   
        
        % Pool of STD without beginners in block (3 first trials) 
        pool_std_w_no_beginners = setdiff(pool_std,begining_of_block);
        
        % Trials which were already selected 
        idx_std_already_included = setdiff(target_indices_std, target_indices_std(idx_std_stepB)) ; 
     
        % Trials to add to balance the number of trial to the same number
        % as DEV
        idx_to_add = pool_std_w_no_beginners(randperm(length(pool_std_w_no_beginners),ntrials));
        
        % Select trial : 1) randomly a number = ntrial and 2) those which
        % were already selected and 'good'
        [EEG_STD_ALL,~] = pop_selectevent(EEG,'event',[idx_std_already_included idx_to_add]);
end