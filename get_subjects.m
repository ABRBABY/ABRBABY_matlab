function [subjects] = get_subjects(indir,OPTIONS)
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

if ~isempty(OPTIONS) 

    % If participant selection is based on suffix
    if isfield(OPTIONS.suffix)
    
        subjects = subjects(contains(subjects,OPTIONS.suffix)) ;
    
    % If participant selection is based on a table
    elseif isfield(OPTIONS.file)
    
        % Reads a .csv file that contains list of subjects to process
        d = readtable(OPTIONS.file) ;
        subjects = table2array(d) ;
    
    end

end
