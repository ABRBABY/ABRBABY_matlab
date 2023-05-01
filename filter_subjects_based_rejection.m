function [subjects] = filter_subjects_based_rejection(subjects, thresh, OPTIONS) 
% Filters subjects depending on rejection rate threshold and returns list
% of kept subjects for group analyses.
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project


%Export trial rejection information for all participants
std_rejected = zeros(size(subjects));
dev1_rejected = zeros(size(subjects));
dev2_rejected = zeros(size(subjects));

%Loop through subjects
for ss=1:length(subjects)

     %Get the subject .csv file with rejected trials information 
    tmp = dir(fullfile(OPTIONS.indir, subjects{ss},strcat(subjects{ss},'*',OPTIONS.params,'.csv'))) ; 

    %Open for reading
    temp = readtable(fullfile(tmp.folder, tmp.name));
    
    % Get and store rejection rates per condition
    std_rejected(ss) =  sum((temp.rejected==1)&strcmp(temp.condition,'STD'));
    dev1_rejected(ss) =  sum((temp.rejected==1)&strcmp(temp.condition,'DEV1'));
    dev2_rejected(ss) =  sum((temp.rejected==1)&strcmp(temp.condition,'DEV2'));
end

%Create table with rejection information at group level
rej_info_all = table(subjects,std_rejected,dev1_rejected,dev2_rejected);

% Save this table into a csv file (use function writetable)
writetable(rej_info_all,fullfile(OPTIONS.indir, strcat('trial_rejected_ERPs_summary', OPTIONS.param, '.csv')));

% Get list of kept subjects as a function of the threshold
not_rejected = (std_rejected+dev1_rejected+dev2_rejected)/height(temp)<thresh ; 

% Display list of rejected subjects in command window
if sum(not_rejected ==0) 
    fprintf(sprintf('\nSubjects rejected : %s\n',subjects{not_rejected==0}));
end

% Output = kept subjects
subjects = subjects(not_rejected);

end
