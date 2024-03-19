%% Variables needed to run the code

function [eeglab_path, biosig_installer_path, erplab_path, BT_toolbox] = get_custom_path()

current_path = pwd ; 

eeglab_path = fullfile(current_path,'eeglab2021.1');
biosig_installer_path = fullfile(current_path,'biosig4octmat-3.8.0','biosig_installer.m');
erplab_path = fullfile(current_path,'erplab8.30');
BT_toolbox = fullfile(current_path,'ToolBox_BrainStem','BT_2013') ;  % bt_toolbox path 
