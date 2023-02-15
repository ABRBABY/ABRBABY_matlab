function [flag_sub_to_create, counter] = test_existance_of_params_in_db(OPTIONS, suffix)
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
    fprintf(strcat(subjects{jj}, '...\n'));

    [does_exist, count(jj)] = check_exist_set_params(subjects{jj}, strcat('*',suffix,'*'), OPTIONS) ; 
    flag_sub_to_create(jj) = ~does_exist ; 
      
end

% Here if count for file to process are not the same : error 
if ~all(count(~flag_sub_to_create) == max(count))
    error('Something wrong in the databse');
else 
    counter = max(count) ;
end

end

%--------------------------------------------------------------
% FUNCTION that check if that sets of param exist 
%--------------------------------------------------------------
function [does_exist, count] = check_exist_set_params(subject, suffix, OPTIONS)

% Reads all folders that are in indir 
d = dir(fullfile(OPTIONS.indir,subject,strcat(suffix,'.set'))); 

% No file exists 
if isempty(d) ; does_exist=0 ; count =1 ; return ; end

for ff=1:length(d) 
    
    EEG = pop_loadset('filepath',fullfile(OPTIONS.indir,subject, d(ff).name),'loadmode','info') ;  
  
    if isequal(suffix,'*_reref_filtered_epoched_RFE*')
        history_name = 'RFE\d*' ;
        history = EEG.history_rfe ;
    elseif isequal(suffix,'*_select_trials_ST*')
        history_name = 'ST\d*';
        history = EEG.history_st ; 
    end
   
    % Check if the set of param correspond to the current file
    if isequal(history,OPTIONS)
        does_exist = 1 ; 
        tmp = regexp(d(ff).name,history_name,'Match');
        tmp2 =  regexp(tmp,'\d*','Match');
        count = str2num(cell2mat(tmp2{:})); 
        return ; 
    end
    
end

% At this point the set of params does not exist and a new file needs to be
% created (with a count increment)
does_exist = 0 ; 
count = length(d) +1;

end
