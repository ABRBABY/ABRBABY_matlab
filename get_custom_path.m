%% Variables needed to run the code

function [eeglab_path, biosig_installer_path, erplab_path, BT_toolbox] = get_custom_path()

current_path = pwd ; 

eeglab_path = fullfile(current_path,'eeglab2021.1');
biosig_installer_path = fullfile(current_path,'biosig4octmat-3.8.0/biosig_installer.m');
erplab_path = fullfile(current_path,'erplab8.30');
BT_toolbox = fullfile(current_path,'ToolBox_BrainStem\BT_2013') ;  % bt_toolbox path 


%name is the initials of the user ('EH' or 'ASD')

% if strcmp(name,'EH')
%     
%     eeglab_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1';
%     biosig_installer_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1\plugins\Biosig3.7.9\biosig_installer.m';
%     erplab_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\erplab8.30';
%     indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';
%     %indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';
%     plot_dir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\png_plots_eeg_data';
%     BT_toolbox = '\\Filer\home\Invites\hervé\Mes documents\GitHub\ABRBABY_matlab\ToolBox_BrainStem\BT_2013' ;  % bt_toolbox path 
%     
% elseif strcmp(name,'ASD')
%    
%     eeglab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/eeglab2021.1' ;
%     biosig_installer_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/biosig4octmat-3.8.0/biosig_installer.m' ;
%     erplab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/erplab8.30';
%     indir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/' ;
%     plot_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder';
%     BT_toolbox = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY_matlab/ToolBox_BrainStem/BT_2013';
% end
% 


    