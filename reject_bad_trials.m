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

suffix_stepA = strrep(RFE,'_','') ; 

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

        % Prepare the header of the column to write in
        % trial_description.txt
        header = char(strcat(subjects{ii},'_',OPTIONS.analysis,'_autorej_low_',num2str(abs(rej_low)),'_high_',num2str(rej_high),'_',suffix_stepA,suffix,num2str(count))) ;        
        
        % Identifies (flag) the first 3 events in blocks 
        begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
        
        % Init a vector 6000 trials
        init_beg_bloc = ones(1,height(T1));
        
        % Get indices of Ssum(TD, DEV1, DEV2 
        flag_sequence =  ~matches(T1.condition,'HF') ; 
         
        % Update flag values 
        init_beg_bloc((T1{:,3}~=0)&~flag_sequence) = ~ismember(find(T1{~flag_sequence,3}~=0),begining_of_block);
        add_flag_column_trials_description(fname_trial_desc, 'begining_block',init_beg_bloc,1);

        % Flip flag if not HF (Condition DEV1, DEV2 or STDD)
        if strcmp(OPTIONS.analysis,'FFR') ; flag_sequence  = ~flag_sequence ; end;

        % Identifies (flag) the automatically rejected trials
        idx_rejacq = find(contains(T1.Properties.VariableNames,'rejection_acq')) ; 
        
        % Get indices of HF trials not rejected at acquisition
        idx_EEG_Noacqrej = find(flag_sequence.*T1{:,idx_rejacq}) ;

        if length(idx_EEG_Noacqrej)~=EEG.trials
            error('The .set file does not contain the same number of trial than in _trial_description.txt');
        end

        % Get indices of the trials which were rejected (in the EEG structure referential)
        [~,idx_rejected] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
        % Retrieve indices of rejected trials on the 6000 referential
        idx_rejected_6000_ref = idx_EEG_Noacqrej(idx_rejected) ;

        % Init a vector 6000 trials
        auto_rejected = ones(1,height(T1));
        
        auto_rejected(~flag_sequence) =  0 ; 

        % Update flag values with trial automatically rejected 
        auto_rejected(idx_rejected_6000_ref) = 0;
    
        % Update flag values with trial wich were rejected at acquisition
        auto_rejected = auto_rejected.*T1{:,idx_rejacq}';
   
        if ~strcmp(OPTIONS.analysis,'FFR')
   
            % Get indices of these rejected events in the EEG structure referential
            all_rej_and_begining_bloc = auto_rejected.*init_beg_bloc ;
    
            % % Select only kept events
            [EEG_DEV1,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_Noacqrej,find(all_rej_and_begining_bloc))),'type','DEV1');
            [EEG_DEV2,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_Noacqrej,find(all_rej_and_begining_bloc))),'type','DEV2');
            [EEG_STD,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_Noacqrej,find(all_rej_and_begining_bloc))),'type','STD');
        
            %% Create a custom history variable to keep track of OPTIONS in each
            % .set saved
            EEG_DEV1.history_stepB = OPTIONS ;
            EEG_DEV2.history_stepB = OPTIONS ;
            EEG_STD.history_stepB = OPTIONS ;
             
            DEV1_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV1_',opt_balance,'_',suffix_stepA,suffix,num2str(count));
            DEV2_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_DEV2_',opt_balance,'_',suffix_stepA,suffix,num2str(count));
            STDD_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_STDD_',opt_balance,'_',suffix_stepA,suffix,num2str(count));
            
            % Save datasets 
            pop_newset(ALLEEG, EEG_DEV1, 1, 'setname',DEV1_fname,'savenew', fullfile(subDir,DEV1_fname),'gui','off');
            pop_newset(ALLEEG, EEG_DEV2, 1, 'setname',DEV2_fname,'savenew', fullfile(subDir, DEV2_fname),'gui','off');
            pop_newset(ALLEEG, EEG_STD, 1, 'setname',STDD_fname,'savenew', fullfile(subDir, STDD_fname),'gui','off');
    
        else
             % Select only kept events
            [EEG_HF,~] = pop_selectevent(EEG,'event', find(ismember(idx_EEG_Noacqrej,find(auto_rejected))),'type','HF');
            
            % Create a custom history variable to keep track of OPTIONS in each .set saved
            EEG_HF.history_stepB = OPTIONS ;
             
            HF_fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_',suffix_stepA,suffix,num2str(count));
            
            % Save datasets 
            pop_newset(ALLEEG, EEG_HF, 1, 'setname',HF_fname,'savenew', fullfile(subDir,HF_fname),'gui','off');
            
        end

        % Update trial_description.txt
        add_flag_column_trials_description(fname_trial_desc, header,auto_rejected,1);
       
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
