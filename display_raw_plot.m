%function display_raw_plot(indir,subject_of_interest)
%%Script for opening bdf files and display plot of raw data

indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_to_look';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_excluded';
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data_issues';
% indir = 'E:\sauvegarde_data\DEVLANG_DATA_excluded' ;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

display_erg_channel = 0 ;    % 1 if want to display ERG channel, 0 otherwise
subject_of_interest = {'DVL_056_T10'} ;

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
    
    if display_erg_channel == 0
    EEG = pop_select( EEG, 'nochannel',{'Erg1'});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',strcat(filename,'_noerg'),'gui','off');
    end 

    EEG = eeg_checkset( EEG );
    pop_eegplot( EEG, 1, 1, 1);

    eeglab redraw ;
    
    end

    
    
%end