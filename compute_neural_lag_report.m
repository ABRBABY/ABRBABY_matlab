
function [] = compute_neural_lag_report(subjects_to_process, OPTIONS)

addpath(OPTIONS.BT_toolbox);

for ss = 1:length(subjects_to_process) %for each subject
    
    
    FFR_file = fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_',OPTIONS.params,'_abr_shifted_data_HF.avg')) ;
    
    [LAG_atmaxmin, maxmincor, all_corrs, all_lags] = bt_xcorrelation2(FFR_file, OPTIONS.stim, OPTIONS.start, OPTIONS.stop, OPTIONS.lagstart, OPTIONS.lagstop, OPTIONS.polarity, OPTIONS.chan, OPTIONS.chancomp) ; 
    
    subject_lag(ss) =  LAG_atmaxmin ;

    % Classify subjects into groups
    if contains(subjects_to_process{ss},OPTIONS.grpA)
        group{ss} = 'A';
    elseif contains(subjects_to_process{ss},OPTIONS.grpB)
        group{ss} = 'B';
    end
    
end

% Write a table with all lags 
neural_lags = table(subjects_to_process, subject_lag', group', 'VariableNames', {'suject_ID', 'neural_lag', 'group'}) ;
writetable(neural_lags,fullfile(OPTIONS.indir,'all_neural_lags.csv'), 'WriteVariableNames', true) ;
