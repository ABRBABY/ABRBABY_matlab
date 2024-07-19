function [subjects_to_analyse] = reject_participants_FFR(subjects, OPTIONS)

% Load .csv with rejection information by participant to merge all rejection information in one single .csv. 
% Then process rejection
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% Export trial rejection information for all participants
HF_rejected = zeros(size(subjects));

%Loop through subjects to extract number of trials rejected
for jj=1:length(subjects)
    %Get individual .csv file with rejected trials information
    csv_name= fullfile(OPTIONS.indir,subjects{jj},strcat(subjects{jj}, OPTIONS.suffix_csv, OPTIONS.param, '.csv'));

    %Open for reading
    temp= readtable(csv_name);
    HF_rejected(jj) = sum((temp.rejected==1)&strcmp(temp.condition,'HF'));
end

% Create table to store these information
rej_info_all = table(subjects,HF_rejected);

% Save this table into a csv file (use function writetable)
writetable(rej_info_all,fullfile(OPTIONS.indir, strcat('trial_rejected_FFR_summary', OPTIONS.param, '.csv')));

%Visualize % of trials rejected for whole subjects
if OPTIONS.visu == 1
    visualize_rejection_rates_FFR(subjects, OPTIONS) ;
end

% Get list of participants rejected based on number of trials rejected
%rej_thresh = 5100 - OPTIONS.threshold ;
rej_trials = subjects(rej_info_all.HF_rejected>=OPTIONS.threshold) ;
% Display list of subjects rejected based on number of trials
if size(rej_trials,1)>0
    for ii = 1:size(rej_trials,1)
        fprintf(sprintf('\nSubjects rejected (number of rej trials): %s\n',cell2mat(rej_trials(ii))));
    end
end

% Get list of participants rejected based on neural lag
neural_table = readtable(fullfile(OPTIONS.indir, strcat('all_neural_lags_',OPTIONS.ffr_polarity, '_ffr_',OPTIONS.polarity,'_corr.csv')),  'Delimiter', ',') ;
rej_neural = subjects(neural_table.neural_lag<OPTIONS.neural_lag) ;
% Display list of subjects rejected based on neural lag value
if size(rej_neural,1)>0
    for ii = 1:size(rej_neural,1)
    fprintf(sprintf('\nSubjects rejected (neural lag): %s\n',cell2mat(rej_neural(ii))));
    end
end

% Get list of participants kept based on the two preceding criteria
% Remove participants from one list if they have been detected for both criteria
rej_neural(contains(rej_neural,rej_trials)) = [] ;
% Remove participants rejected from number of trials
subjects_to_analyse = subjects ;
subjects_to_analyse(contains(subjects_to_analyse,rej_trials)) = [];
% Remove participants rejected from neural lag
subjects_to_analyse(contains(subjects_to_analyse,rej_neural)) = [];

% Create .csv file to store information about rejected participants
rej_ffr_info = table(subjects,contains(subjects,rej_trials), contains(subjects, rej_neural), 'VariableNames',{'subjects','rej_trials','rej_neural'}) ;
writetable (rej_ffr_info, fullfile(OPTIONS.indir,'participants_rejection_ffr.csv')) ;
end
