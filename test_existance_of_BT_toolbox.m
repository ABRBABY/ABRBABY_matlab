function [flag_sub_to_create] = test_existance_of_BT_toolbox(OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Loop through subjects
for jj=1:length(subjects) 

    flag_sub_to_create(jj) = ~exist(fullfile(OPTIONS.indir,subjects{jj},'BT_toolbox_formatted'),'dir'); 
      
end
