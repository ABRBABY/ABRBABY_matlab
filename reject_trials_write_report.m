function [] = reject_trials_write_report(preproc_filenames,eeg_elec, win_of_interest, rej_low, rej_high) 
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Loop though subjects
for ii=1:length(preproc_filenames)

    [filepath,filename,ext] = fileparts(preproc_filenames{ii}) ;

    EEG = pop_loadset(strcat(filename,ext),filepath) ;

    % Get indices of the trials which were rejected (without messing around with the relative indices)
    [EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,eeg_elec,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);

    % Extract variables of interest
    trial_index = 1:EEG.trials;
    trial_num = [EEG.event.urevent];
    condition = {EEG.event.type} ;
    latency = [EEG.event.latency]/EEG.srate;  
    rejected = ismember(trial_index,idx_rejected_all) ; 
    bloc = repelem(1:30,30) ; % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
    
    % Create table to store these information
    list_trial_infos = table(trial_index',condition',latency', trial_num',rejected', bloc',...
        'VariableNames', {'trial_index', 'condition', 'latency','trial_num','rejected','bloc'}) ;

    %  Save this table into a csv file (use function writetable)
    writetable(list_trial_infos,...
        fullfile(filepath,strcat(filename,'_low_',num2str(rej_low),'_high_',num2str(rej_high),'infos_trials.csv')),...
        'WriteVariableNames', true) ; 

end
