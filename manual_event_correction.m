% ERPs sanity check script -
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project
clear all ;
indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data' ;
subject_name = 'DVL_059_T8' ;

%% Open files

% Reads all folders that are in indir
d = dir(indir);
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Check that subject is in folder
if ~ismember(subject_name,subjects)
    error('Subject %s is not in folder.', subject_name)
end

% Inititalize output parameter
out_filenames = [] ;

% Printout the id of the subject in console
fprintf(strcat(subject_name, '...\n'));

% Get trials_description file
fname_txt= dir(fullfile(indir,subject_name, strcat(subject_name, '_trials_description.txt')));

% Error if trials_description file does not exist
if isempty(fname_txt)
    error('trials_description.txt file does not exist for participant %s.', subject_name)
end

% Open EEGLAB
[ALLEEG, CURRENTSET, ALLCOM] = eeglab;

% Select bdf file in the folder and open it
fname_bdf = dir(fullfile(indir,subject_name, strcat(subject_name, '*.bdf')));
EEG = pop_biosig(fullfile(indir, subject_name, fname_bdf.name));

% Save a first dataset in EEGLAB
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',subject_name,'gui','off');

% Open trials_description file
trials_descr = readtable(fullfile(fname_txt.folder, fname_txt.name)) ;

% Save description file with new name
writetable(trials_descr, fullfile(fname_txt.folder, strrep(fname_txt.name,'trials_description', 'OLD_trials_description_')), 'WriteVariableNames',true) ;

% Add filepath in EEG structure for detect_and_create_report function
EEG.filepath = fname_txt.folder ;

% Reject bad events
EEG_rejected = detect_events_and_create_report(EEG, indir, '_ergstim.txt', {'Erg1'}) ;
%% Change manually values of rejection variable

% ######### If last XX event to remove #########

% XX = ... ;
% remove = (6000 - XX) + 1 ;
% trials_descr{remove:6000,3} = 0 ;

% XX = 2604 ;
% remove = (6000 - XX) + 1 ;
% trials_descr{remove:6000,3} = 0;

% ######### If first XX event to remove #########

% XX = ... ;
% trials_descr{1:XX,3} = 0 ;

% XX = 18 ;
% trials_descr{1:XX,3} = 0 ;

% ######### If events XX to YY to remove #########

% XX = ... ;
% YY = ... ;
% trials_descr{XX:YY,3} = 0 ;

% XX = 4407 ;
% YY =  4416;
% trials_descr{XX:YY,3} = 0 ;
% 
% % Check number of kept trials
% sum(trials_descr{:,3}) 

%% Retrieve latencies

% ######### If last XX event to remove #########

% XX = ... ;
% remove = (6000 - XX) + 1 ;
% trials_descr{1:removed,2} = [EEG_rejected.event.latency]' ;

% XX = 2604 ;
% removed = 6000 - XX ;
% trials_descr{1:removed,2} = [EEG_rejected.event.latency]' ;

% ######### If first XX event to remove #########

% XX = ... ;
% remove = XX + 1;
% trials_descr{remove:6000,2} = [EEG_rejected.event.latency]' ;

% XX = 18 ;
% remove = XX + 1 ;
% trials_descr{remove:6000,2} = [EEG_rejected.event.latency]' ;

% ######### If events XX to YY to remove #########
% Depends on rejected indices... case by case
% XX = 4407 ;
% trials_descr{1:XX-1,2} = [EEG_rejected.event(1:XX-1).latency]' ;
% trials_descr{4417:end,2} = [EEG_rejected.event(4418:end).latency]' ;
% 
% % Check sizes
% size(trials_descr{4417:end,2})
% size([EEG_rejected.event(4418:end).latency]')

%% Save new file_description

% Save new description file
writetable(trials_descr, fullfile(fname_txt.folder, fname_txt.name), 'WriteVariableNames',true) ;
