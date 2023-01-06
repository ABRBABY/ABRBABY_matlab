function [ALLEEG] = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) 
% Prepare the environnement : 
% Set EEGLAB path, set ERPLAB path, install biosig
% A.-Sophie Dubarry - 2022 

% Add ERPLAB to Matlab path
addpath(genpath(erplab_path));

% Save current path
tmp = pwd ; 

% Move to EEGLAB directory
cd(eeglab_path) ; 

% Start eeglab (this add EEGLAB to Matlab path) 
[ALLEEG, ~, ~, ~] = eeglab;

% Install biosig
run(biosig_installer_path) ; 

% Return to initial path
cd(tmp) ;