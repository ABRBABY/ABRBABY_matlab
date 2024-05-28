function [neural_lag] = compute_neural_lag(OPTIONS, flag_sub_to_create)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly computes neurl lab based on BT toolbox function
% bt_xcorrelation2_abrbaby and create a csv table 

% BT_Toolbox parameters 
CHAN =1;        % Number of channels 
CHANCOMP =1;    % Number of stimuli files?

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Only keeps subjects to process
subjects_to_process = subjects(flag_sub_to_create) ; 

for ss=1:length(subjects_to_process) %for each subject
    
    % Create a folder for files specific to BT_toolbox
    BT_folder = fullfile(OPTIONS.indir, subjects_to_process{ss},'BT_toolbox_formatted');
    fname_avg = fullfile(BT_folder,strcat(subjects_to_process{ss},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.avg')) ;
    
    % Computes the neural lag 
    [neural_lag(ss), idx_neural_lag, maxmincor, all_corrs, all_lags] = bt_xcorrelation2_abrbaby(fname_avg, OPTIONS.stim, OPTIONS.start, OPTIONS.stop, OPTIONS.lagstart, OPTIONS.lagstop, OPTIONS.polarity, CHAN , CHANCOMP) ; 
    
    if OPTIONS.display == 1
        figure; plot(all_lags,all_corrs,'r') ; hold on ; plot(all_lags(idx_neural_lag),all_corrs(idx_neural_lag),'*b') ; title(strrep(subjects_to_process{ss},'_','-')) ; grid on ; xlabel('Neural lag (ms)'); ylabel('Pearsons correlation') ; 
        legend('All correlations (time shifted)','Neural lag','Fontsize',12);
        % Save figure
         print('-dpng',fullfile('E:\EEG\DATA\plot_dir_DEVLANG_data\png_folder\ffr_neural_lags_r_pearsons', strcat(subjects{ss},'_FFR_',OPTIONS.params, '_', num2str(OPTIONS.lagstart),'_', num2str(OPTIONS.lagstop)))); 
    end

end

if OPTIONS.table == 1
    % Write a table with all lags 
    neural_lags = table(subjects_to_process, neural_lag', 'VariableNames', {'suject_ID', 'neural_lag'}) ;
    writetable(neural_lags,fullfile(OPTIONS.indir,strcat('all_neural_lags_',OPTIONS.ffr_polarity, '_ffr_',OPTIONS.polarity,'_corr_', num2str(OPTIONS.lagstart),'_',num2str(OPTIONS.lagstop),'_',OPTIONS.params, '.csv')), 'WriteVariableNames', true) ;
    % Display end message
    warning('Neural lag table is saved.') ;
end