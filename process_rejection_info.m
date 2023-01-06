function [] = process_rejection_info(indir, suffix) 
% Load .csv with rejection information by participant to merge all rejection information in one single .csv. 
%This new table is then visualized with the rejected_trials_visu_ERPs function
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Export trial rejection information for all participants
std_rejected = zeros(size(subjects));
dev1_rejected = zeros(size(subjects));
dev2_rejected = zeros(size(subjects));
%Loop through subjects
for jj=1:length(subjects)
    %Get individual .csv file with rejected trials information
    csv_name= fullfile(indir,subjects{jj},strcat(subjects{jj}, suffix));
    %Open for reading
    temp = readtable(csv_name);
    
    std_rejected(jj) =  sum((temp.rejected==1)&strcmp(temp.condition,'STD'));
    dev1_rejected(jj) =  sum((temp.rejected==1)&strcmp(temp.condition,'DEV1'));
    dev2_rejected(jj) =  sum((temp.rejected==1)&strcmp(temp.condition,'DEV2'));
end

% Create table to store these information
rej_info_all = table(subjects,std_rejected, dev1_rejected, dev2_rejected);
% Save this table into a csv file (use function writetable)
writetable(rej_info_all,fullfile(indir, 'trial_rejected_summary.csv'));

%Display barplots and histograms of trial rejection
rejected_trials_visu_ERPs(rej_info_all,sum(strcmp(temp.condition,'STD')),sum(strcmp(temp.condition,'DEV1')));

end
