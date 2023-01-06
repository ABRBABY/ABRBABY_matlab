function [] = select_and_save_trials_per_condition(ALLEEG, preproc_filenames, eeg_elec, win_of_interest, rej_low, rej_high, opt_balance)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Check balance option
if strcmp(opt_balance,'balanced')
    STD_number = 1 ;
elseif strcmp(opt_balance,'unbalanced')
    STD_number = 2 ;
else 
    error('Unknow option : choose ''balanced'' or ''unbalanced''') ;
end

% Loop though subjects
for ii=1:length(preproc_filenames)

    [filepath,filename,ext] = fileparts(preproc_filenames{ii}) ;

    EEG = pop_loadset(strcat(filename,'.set'),filepath) ;

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
    end
    
    [EEG_STD1,target_indices_std1] = pop_selectevent(EEG,'event',idx_std1);
    [EEG_STD2,target_indices_std2] = pop_selectevent(EEG,'event',idx_std2);
    
    [EEG_STD1_thresh,idx_std1_rej] = pop_eegthresh(EEG_STD1,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_STD2_thresh,idx_std2_rej] = pop_eegthresh(EEG_STD2,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV1_thresh,idx_dev1_rej] = pop_eegthresh(EEG_DEV1,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV2_thresh,idx_dev2_rej] = pop_eegthresh(EEG_DEV2,1,eeg_elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
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
            %Find indices of rejected DEV1 in the 900 trials referential
            list_of_rej_dev1 = target_indices1(idx_dev1_rej);
            %Select random indices amont the rejected DEV1
            random_std1_to_remove = randsample(list_of_rej_dev1-1,size(idx_dev1_rej,2)-size(idx_std1_rej,2));
            %Remove STD corresponding to rejected DEV1
            idx_std_already_included = setdiff(target_indices_std1, target_indices_std1(idx_std1_rej)) ;
            new_list_of_std1 = setdiff(idx_std_already_included,random_std1_to_remove);
            %Modify EEG struct
            [EEG_STD1_thresh,~] = pop_selectevent(EEG,'event',[new_list_of_std1]);
        end
        
        % If nubmber of STD2 > number of DEV2 : delete random STD preceding
        % missing DEV2
        if length(EEG_DEV2_thresh.event)<length(EEG_STD2_thresh.event)
            %Find indices of rejected DEV1 in the 900 trials referential
            list_of_rej_dev2 = target_indices2(idx_dev2_rej);
            %Select random indices amont the rejected DEV1
            random_std2_to_remove = randsample(list_of_rej_dev2-1,size(idx_dev2_rej,2)-size(idx_std2_rej,2));
            %Remove STD corresponding to rejected DEV1
            idx_std_already_included = setdiff(target_indices_std2, target_indices_std2(idx_std2_rej)) ;
            new_list_of_std2 = setdiff(idx_std_already_included,random_std2_to_remove);
            %Modify EEG struct
            [EEG_STD2_thresh,~] = pop_selectevent(EEG,'event',[new_list_of_std2]);
        end
      
    end

    % Save datasets 
    pop_newset(ALLEEG, EEG_DEV1_thresh, 1, 'setname',strcat(filename,'_','EEG_DEV1_thresh_',opt_balance),'savenew', fullfile(filepath, strcat(filename,'_DEV1_thresh_',opt_balance)),'gui','off');
    pop_newset(ALLEEG, EEG_DEV2_thresh, 1, 'setname',strcat(filename,'_','EEG_DEV2_thresh_',opt_balance),'savenew', fullfile(filepath, strcat(filename,'_DEV2_thresh_',opt_balance)),'gui','off');
    pop_newset(ALLEEG, EEG_STD1_thresh, 1, 'setname',strcat(filename,'_','EEG_STD1_thresh_',opt_balance),'savenew', fullfile(filepath, strcat(filename,'_STD1_thresh_',opt_balance)),'gui','off');
    pop_newset(ALLEEG, EEG_STD2_thresh, 1, 'setname',strcat(filename,'_','EEG_STD2_thresh_',opt_balance),'savenew', fullfile(filepath, strcat(filename,'_STD2_thresh_',opt_balance)),'gui','off');

end

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