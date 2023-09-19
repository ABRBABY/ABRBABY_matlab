%function display_raw_plot(indir,subject_of_interest)
%%Script for opening bdf files and display plot of raw data

indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_to_look';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_excluded';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data_issues';
% indir = 'E:\sauvegarde_data\DEVLANG_DATA_excluded' ;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%subject_of_interest = {'DVL_013_T10', 'DVL_013_T8', 'DVL_014_T24', 'DVL_015_T18', 'DVL_015_T24', 'DVL_016_T18', 'DVL_016_T24', 'DVL_017_T18', 'DVL_018_T10', 'DVL_018_T6', 'DVL_018_T8'} ;
% subject_of_interest = {'DVL_035_T18','DVL_045_T10','DVL_027_T24','DVL_043_T6'} ;
subject_of_interest = {'DVL_016_T18'} ;
% subject_of_interest = {'DVL_test_a135','DVL_test_a240','DVL_test_a300','DVL_test_b135','DVL_test_b245','DVL_test_b300'} ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Loop through subjects
    for ind = 1:length(subject_of_interest)
     
    subject = subject_of_interest{ind} ;
    fprintf(strcat(subject, '...\n'));
     
    %% IMPORT
    % Get BDF file
    fname= dir(fullfile(indir,subject,'*.bdf'));
 
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subject, fname.name));
    
    % Save a first dataset in EEGLAB 
    [~,filename,~] = fileparts(fname.name);    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',filename,'gui','off');
    
    EEG = eeg_checkset( EEG );
    pop_eegplot( EEG, 1, 1, 1);

    eeglab redraw ;
    
    end

    
    
%end