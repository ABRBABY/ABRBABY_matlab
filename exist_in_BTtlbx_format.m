function [flag_sub_to_create] = exist_in_BTtlbx_format(OPTIONS,output_suffix)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Loop through subjects
for jj=1:length(subjects) 

    flag_sub_to_create(jj) = ~isempty(dir(fullfile(OPTIONS.indir,subjects{jj},'BT_toolbox_formatted',strcat('*',output_suffix,'*')))); 
      
end
