function [] = compute_neural_lag_report(subjects_to_process, OPTIONS)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly computes neurl lab based on BT toolbox function
% bt_xcorrelation2_abrbaby and create a csv table 

addpath(OPTIONS.BT_toolbox);

for ss = 1:length(subjects_to_process) %for each subject
    
    
    FFR_file = fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_',OPTIONS.params,'_abr_', OPTIONS.ffr_polarity,'_shifted_data_HF.avg')) ;
    
    %[LAG_atmaxmin, maxmincor, all_corrs, all_lags] = bt_xcorrelation2(FFR_file, OPTIONS.stim, OPTIONS.start, OPTIONS.stop, OPTIONS.lagstart, OPTIONS.lagstop, OPTIONS.polarity, OPTIONS.chan, OPTIONS.chancomp) ; 
    [LAG_atmaxmin, maxmincor, all_corrs, all_lags] = bt_xcorrelation2_abrbaby(FFR_file, OPTIONS.stim, OPTIONS.start, OPTIONS.stop, OPTIONS.lagstart, OPTIONS.lagstop, OPTIONS.polarity, OPTIONS.chan, OPTIONS.chancomp) ; 

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
writetable(neural_lags,fullfile(OPTIONS.indir,strcat('all_neural_lags_',OPTIONS.ffr_polarity, '_ffr_',OPTIONS.polarity,'_corr_', num2str(OPTIONS.lagstart),'_',num2str(OPTIONS.lagstop),'_',OPTIONS.params, '.csv')), 'WriteVariableNames', true) ;

% Display end message
disp('Neural lag table is saved.') ;