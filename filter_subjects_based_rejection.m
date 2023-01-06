function [subjects] = filter_subjects_based_rejection(indir, thresh) 
% Filters subjects depending on rejection rate threshold and returns list
% of kept subjects for group analyses.
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

suffix =  '_low_-150_high_150infos_trials.csv';

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

subjects = subjects((std_rejected+dev1_rejected+dev2_rejected)/height(temp)<thresh);

end
