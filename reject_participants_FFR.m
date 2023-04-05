function [subjects_to_analyse] = reject_participants_FFR(subjects, OPTIONS)

% Load .csv with rejection information by participant to merge all rejection information in one single .csv. 
% Then process rejection
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%Export trial rejection information for all participants
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

% Get list of participants kept based on threshold
rej_thresh = 5100 - OPTIONS.threshold ;
subjects_to_analyse = subjects(find(rej_info_all.HF_rejected<rej_thresh));


end
