function [subjects] = get_all_subjects(indir)
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

