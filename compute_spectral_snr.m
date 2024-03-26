function [] = compute_spectral_snr(OPTIONS,flag_sub_to_create)
% 
% Converts the ABR signal into BT_toolbox readable format + optionnal display 
% 
% Estelle Herve, A.-Sophie Dubarry - 2024 - %80PRIME Project
%
% This function mainly do : 

% Reads all folders that are in indir 
d = dir(OPTIONS.indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

suffix_stepA = strrep(RFE,'_','') ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 
