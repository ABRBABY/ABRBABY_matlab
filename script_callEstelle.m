
% HERE I CAN EXECUTE ANY FUNCTION OF THE PACKAGE WITH MY PATHS

%% Variables to enter manually before running the code
eeglab_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1';
biosig_installer_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1\plugins\Biosig3.7.9\biosig_installer.m';
chan_dir ='\\Filer\\home\\Invites\\hervé\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc';
plots_dir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\png_plots_eeg_data';
indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';
%indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';

%STD_number = 1; %balanced number of trials
%STD_number = 2; %unbalanced number of trials (all STDs except the first 3 of each bloc)

for STD_number = 1:2
    
    %% Sanity checks:
    
    %abrbaby_process_ERP_sanity_exportdata(eeglab_path, biosig_installer_path,indir) ;
    %to delete
    
    %abrbaby_process_ERP_sanity_exportdata_allSTD(eeglab_path, biosig_installer_path,indir,plots_dir,chan_dir,STD_number);
    %modified this function so it considers STD_number
    
    %abrbaby_sanity_check_FFR(eeglab_path, biosig_installer_path,indir,plots_dir);
     
    %process_rejection_info(indir, '_low_-150_high_150infos_trials.csv');
    
    %beep;
    %disp("Done");
    
    %% Analyses:
    
     abrbaby_process_visual_comparison_by_group_age(eeglab_path, biosig_installer_path,indir, STD_setfile, STD_number);
%     
%     FFR_analysis(indir) ;
%     
%     beep;
%     disp("Done");
end


