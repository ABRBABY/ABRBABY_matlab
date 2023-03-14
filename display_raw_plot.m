function display_raw_plot(indir,subject_of_interest)
%%Script for opening bdf files and display plot of raw data

indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data';
%indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_to_look';
%indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';
%indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_DATA_excluded';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

subject_of_interest = ['DVL_027_T24'];

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..


%Loop through subjects
    %for jj=1:length(subjects)
    for jj=find(ismember(subjects, subject_of_interest))  
        
     fprintf(strcat(subjects{jj}, '...\n'));
     
    %% IMPORT
    % Get BDF file
    fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
 
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));
    
    % Save a first dataset in EEGLAB 
    [~,filename,~] = fileparts(fname.name);    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',filename,'gui','off');
    
    EEG = eeg_checkset( EEG );
    pop_eegplot( EEG, 1, 1, 1);

    end
    
end