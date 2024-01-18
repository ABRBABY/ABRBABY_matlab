function [out_filenames] = reject_bad_trials(ALLEEG, OPTIONS, opt_balance, flag_sub_to_create, count, suffix, RFE)
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

% Inititalize output parameter
out_filenames = [] ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

% Loop though subjects
for ii=1:length(subjects)
    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Set rfe file to work on
    file_rfe = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},RFE,'.set'))) ;

    % Error if rfe file does not exist
    if isempty(file_rfe) ; error('File %s does not exist for subject %s', RFE, subjects{jj}); end
    
    %Get filepath
    filepath = file_rfe.folder ;
  
    %Creates resulting filename
    out_filenames{ii} = fullfile(indir,subjects{ii}, strcat(subjects{ii},suffix,num2str(count),'.set')) ; 

    %Load the RFE .set file to work on
    EEG = pop_loadset(strcat(subjects{ii},RFE,'.set'),filepath) ;
    
    %Get eeg_elec and win_of_interest from RFE set of parameters
    eeg_elec = EEG.history_rfe.eeg_elec ;
    win_of_interest = EEG.history_rfe.win_of_interest ;

    % Select trials per conditions
    [EEG_DEV1,target_indices1] = pop_selectevent(EEG,'type','DEV1');
    [EEG_DEV2,target_indices2] = pop_selectevent(EEG,'type','DEV2');
    [EEG_STD,target_indices_std] = pop_selectevent(EEG,'type','STD');
    
    begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
        
    if strcmp(opt_balance,'balanced')
        idx_std1 = target_indices_std(ismember(target_indices_std,target_indices1-1));
        idx_std2 = target_indices_std(ismember(target_indices_std,target_indices2-1));
    elseif strcmp(opt_balance,'unbalanced')
        idx_std1 = setdiff(target_indices_std,begining_of_block);
        idx_std2 = setdiff(target_indices_std,begining_of_block);
    else
        error('Unknow option : choose ''balanced'' or ''unbalanced''') ;
    end
    
    [EEG_STD1,target_indices_std1] = pop_selectevent(EEG,'event',idx_std1);
    [EEG_STD2,target_indices_std2] = pop_selectevent(EEG,'event',idx_std2);
    
    [EEG_STD1_thresh,idx_std1_rej] = pop_eegthresh(EEG_STD1,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_STD2_thresh,idx_std2_rej] = pop_eegthresh(EEG_STD2,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV1_thresh,idx_dev1_rej] = pop_eegthresh(EEG_DEV1,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV2_thresh,idx_dev2_rej] = pop_eegthresh(EEG_DEV2,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
    %% If we want to select the same number of trial both for DEV and STD
    if strcmp(opt_balance,'balanced')

        [EEG_STD_thresh, idx_removed] = pop_eegthresh(EEG_STD,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
        std_good = setdiff(1:900,target_indices_std(idx_removed)); 
       
        % If nubmber of STD1 < number of DEV1 : randomly select other STD
        if length(EEG_DEV1_thresh.event)>length(EEG_STD1_thresh.event)
            % Find number of trial to add in STD
            ntrials =  length(EEG_DEV1_thresh.event)-length(EEG_STD1_thresh.event) ;
            %Apply function to balance number of STD trial to reach the number of DEV trials
            EEG_STD1_thresh = balance_number_of_STD(EEG,ntrials,std_good,target_indices1,begining_of_block,target_indices_std1,idx_std1_rej,idx_std1) ;
            
        end
        
        % If nubmber of STD2 < number of DEV2 : randomly select other STD
        if length(EEG_DEV2_thresh.event)>length(EEG_STD2_thresh.event)
            % Find number of trial to add in STD
            ntrials =  length(EEG_DEV2_thresh.event)-length(EEG_STD2_thresh.event) ;
            %Apply function to balance number of STD trial to reach the number of DEV trials
            EEG_STD2_thresh = balance_number_of_STD(EEG,ntrials,std_good,target_indices2,begining_of_block,target_indices_std2,idx_std2_rej,idx_std2) ;
            
        end
        
        % If nubmber of STD1 > number of DEV1 : delete random STD preceding missing DEV1
        if length(EEG_DEV1_thresh.event)<length(EEG_STD1_thresh.event)
            % Removes randomly the number of STD exceeding the number of DEV
            idx_to_keep = randperm(length(EEG_STD1_thresh.event), length(EEG_DEV1_thresh.event)) ; 
            %Modify EEG struct
            [EEG_STD1_thresh,~] = pop_selectevent(EEG,'event',[target_indices1(idx_to_keep)]);
        end
        
        % If nubmber of STD2 > number of DEV2 : delete random STD preceding
        % missing DEV2
        if length(EEG_DEV2_thresh.event)<length(EEG_STD2_thresh.event)
            % Removes randomly the number of STD exceeding the number of DEV
            idx_to_keep = randperm(length(EEG_STD2_thresh.event), length(EEG_DEV2_thresh.event)) ;            
            %Modify EEG struct
            [EEG_STD2_thresh,~] = pop_selectevent(EEG,'event',[target_indices2(idx_to_keep)]);   
        end
      
    end
    
    %% Create a custom history variable to keep track of OPTIONS in each
    % .set saved
    EEG_DEV1_thresh.history_rej = OPTIONS ;
    EEG_DEV2_thresh.history_rej = OPTIONS ;
    EEG_STD1_thresh.history_rej = OPTIONS ;
    EEG_STD2_thresh.history_rej = OPTIONS ;

    suffix_rfe = strsplit(RFE,'_') ; 
    
    % Save datasets 
    pop_newset(ALLEEG, EEG_DEV1_thresh, 1, 'setname',strcat(subjects{ii},'_','EEG_DEV1_',OPTIONS.analysis,'_',opt_balance,suffix_rfe(end),suffix,num2str(count)),'savenew', fullfile(filepath, strcat(subjects{ii},'_DEV1_',OPTIONS.analysis,'_',opt_balance,'_',suffix_rfe(end),suffix,num2str(count))),'gui','off');
    pop_newset(ALLEEG, EEG_DEV2_thresh, 1, 'setname',strcat(subjects{ii},'_','EEG_DEV2_',OPTIONS.analysis,'_',opt_balance,suffix_rfe(end),suffix,num2str(count)),'savenew', fullfile(filepath, strcat(subjects{ii},'_DEV2_',OPTIONS.analysis,'_',opt_balance,'_',suffix_rfe(end),suffix,num2str(count))),'gui','off');
    if strcmp(opt_balance,'balanced')
         pop_newset(ALLEEG, EEG_STD1_thresh, 1, 'setname',strcat(subjects{ii},'_','EEG_STD1_',OPTIONS.analysis,'_',opt_balance,suffix_rfe(end),suffix,num2str(count)),'savenew', fullfile(filepath, strcat(subjects{ii},'_STD1_',OPTIONS.analysis,'_',opt_balance,'_',suffix_rfe(end),suffix,num2str(count))),'gui','off');
         pop_newset(ALLEEG, EEG_STD2_thresh, 1, 'setname',strcat(subjects{ii},'_','EEG_STD2_',OPTIONS.analysis,'_',opt_balance,suffix_rfe(end),suffix,num2str(count)),'savenew', fullfile(filepath, strcat(subjects{ii},'_STD2_',OPTIONS.analysis,'_',opt_balance,'_',suffix_rfe(end),suffix,num2str(count))),'gui','off');
    elseif strcmp(opt_balance,'unbalanced')
         pop_newset(ALLEEG, EEG_STD1_thresh, 1, 'setname',strcat(subjects{ii},'_','EEG_STDD_',OPTIONS.analysis,'_',opt_balance,suffix_rfe(end),suffix,num2str(count)),'savenew', fullfile(filepath, strcat(subjects{ii},'_STDD_',OPTIONS.analysis,'_',opt_balance,'_',suffix_rfe(end),suffix,num2str(count))),'gui','off');

    end
    % Name of the file report 
    fname = fullfile(filepath,strcat(subjects{ii},'_infos_trials','_low_',num2str(rej_low),'_high_',num2str(rej_high),'_',suffix_rfe(end),suffix,num2str(count),'.csv')) ; 
    
    % Write csv file directly into the subject dir
    produce_report(fname{1}, EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high, opt_balance) ; 

end

end

%--------------------------------------------------------------
% FUNCTION that write a report on rejected trials. Regardless to
% conditions) -> we re-excute pop_eegthresh on all trials 
% (we do not save .set but the report)
%--------------------------------------------------------------
function [] = produce_report(fname,EEG, eeg_elec, bloc, win_of_interest, rej_low, rej_high, opt_balance) 

    if strcmp(opt_balance,'balanced')
        warndlg('The report of trial status (rejected/not rejected) is under developement for ''balanced'' option ','Warning')
        return
    end
     % Get indices of the trials which were rejected (without messing around with the relative indices)
    [~,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);

    % Reject the STD at the very begining of block
    begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 

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
function [EEG_STD_ALL] = balance_number_of_STD(EEG,ntrials,std_good,target_indices,begining_of_block,target_indices_std,idx_std_rej, idx_std)
   
        % Pool of STD without those rejected by threshold detection
        pool_std = setdiff(std_good,target_indices);   
        
        % Pool of STD without beginners in block (3 first trials) 
        pool_std_w_no_beginners = setdiff(pool_std,begining_of_block);
        
        % Trials which were already selected 
        idx_std_already_included = setdiff(target_indices_std, target_indices_std(idx_std_rej)) ; 
     
        % Trials to add to balance the number of trial to the same number
        % as DEV
        idx_to_add = pool_std_w_no_beginners(randperm(length(pool_std_w_no_beginners),ntrials));
        
        % Select trial : 1) randomly a number = ntrial and 2) those which
        % were already selected and 'good'
        [EEG_STD_ALL,~] = pop_selectevent(EEG,'event',[idx_std_already_included idx_to_add]);
end