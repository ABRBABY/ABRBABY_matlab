function [flag_sub_to_create] = test_existance_of_combination(OPTIONS,flag_sub_to_create, count,RFE_num,suffix)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Loop through subjects
for jj=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{jj}, '...\n')) ;

    % Get suffix of first set of parameters
    suffix_param1 = strsplit(RFE_num,'_') ;
    if size(dir(fullfile(OPTIONS.indir,subjects{jj},(strcat(subjects{jj},'*','balanced_',suffix_param1{end},suffix,num2str(count),'.set'))))) < 8
        flag_sub_to_create(jj) = 1 ;
    else 
        flag_sub_to_create(jj) = 0 ;
  
    end

end