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
        
        % Get indices of trials wich are not HF and not rejected at
        % acquisition
        idx_EEG_LF_Noacqrej = find(idx_not_HF.*T1{:,idx_rejacq}) ;

        if length(idx_EEG_LF_Noacqrej)~=EEG.trials
            error('The .set file does not contain the same number of trial than in _trial_description.txt');
        end
        
         % Get indices of the trials which were rejected (in the EEG structure referential)
        [~,idx_rejected] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
        % Retrieve indices of rejected trials on the 6000 referential
        idx_rejected_6000_ref = idx_EEG_LF_Noacqrej(idx_rejected) ;

        % Init a vector 6000 trials
        auto_rejected = ones(1,height(T1));
        auto_rejected(~idx_not_HF) =  0 ; 

        % Update flag values with trial automatically rejected 
        auto_rejected(idx_rejected_6000_ref) = 0;
    
        % Update flag values with trial wich were rejected at acquisition
        auto_rejected = auto_rejected.*T1{:,idx_rejacq}';

        % Update trial_description.txt
        header = char(strcat(subjects{ii},'_infos_trials','_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_',suffix_stepA(end),suffix,num2str(count))) ;        
        add_flag_column_trials_description(fname_trial_desc, header,auto_rejected);

        % Get indices of these rejected events in the EEG structure referential
        all_rej_and_begining_bloc = auto_rejected.*init_beg_bloc ;

        % % Select only kept events
        [EEG_DEV1,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_LF_Noacqrej,find(all_rej_and_begining_bloc))),'type','DEV1');
        [EEG_DEV2,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_LF_Noacqrej,find(all_rej_and_begining_bloc))),'type','DEV2');
        [EEG_STD,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_LF_Noacqrej,find(all_rej_and_begining_bloc))),'type','STD');
    
        %% Create a custom history variable to keep track of OPTIONS in each
        % .set saved
        EEG_DEV1.history_stepB = OPTIONS ;
        EEG_DEV2.history_stepB = OPTIONS ;
        EEG_STD.history_stepB = OPTIONS ;
         
        DEV1_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV1_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
        DEV2_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV2_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
        STDD_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_STDD_',opt_balance,'_',suffix_stepA{end},suffix,num2str(count));
    
        % Save datasets 
        pop_newset(ALLEEG, EEG_DEV1, 1, 'setname',DEV1_fname,'savenew', fullfile(subDir,DEV1_fname),'gui','off');
        pop_newset(ALLEEG, EEG_DEV2, 1, 'setname',DEV2_fname,'savenew', fullfile(subDir, DEV2_fname),'gui','off');
        pop_newset(ALLEEG, EEG_STD, 1, 'setname',STDD_fname,'savenew', fullfile(subDir, STDD_fname),'gui','off');
    
        % Name of the file report 
        fnameReport = fullfile(subDir,strcat(subjects{ii},'_infos_trials','_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_',suffix_stepA(end),suffix,num2str(count),'.csv')) ; 
        fnameTrialDescription = fullfile(subDir,strcat(subjects{ii},'_trials_description.txt'));
        
   end 
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