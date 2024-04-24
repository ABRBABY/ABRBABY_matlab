% Automatic trigger detection

% ------------------- Set environment 
% Variables to enter manually before running the code

% DATA directory 
% indir = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data/EEG_data_revised_by_participant_rejA'; 
indir = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data'; 
% indir = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

% selected_subj = {'DVL_013_T24','DVL_011_T10','DVL_044_T8'};
selected_subj = [];

%Get list of subjects in indir
list_subjects = get_subjects(indir,[]);

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path,bt_toolbox] = get_custom_path();

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path, bt_toolbox) ;

trig = {'Erg1'};
TRIG_MODALITY = '_ergstim.txt';
chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 

% Reads all folders that are in indir
d = dir(indir);
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..te

% Inititalize output parameter
out_filenames = [] ;

if ~isempty(selected_subj)
    % Get selected subjects
    subjects = subjects(ismember(subjects,selected_subj)) ; 
end

%Loop through subjects
for jj=1:length(subjects)

    EventDetection = fullfile(indir,subjects{jj},strcat(subjects{jj},'_trials_description.txt'));
    if ~exist(EventDetection,'file')

        % Printout the id of the subject in console
        fprintf(strcat(subjects{jj}, '...\n'));

        %% IMPORT
        % Get BDF file
        %[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
        fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
        [~,filename,~] = fileparts(fname.name);

        % Select bdf file in the folder
        EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

        % Save a first dataset in EEGLAB
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');

        % Add channels information
        EEG = pop_chanedit(EEG, 'lookup',chan_dir) ;

        % Add filepath in EEG structure
        EEG.filepath = fname.folder ;
        [EEG] = detect_events_and_create_report(EEG, indir, TRIG_MODALITY, trig) ;

    end


end


%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir, hp, lp, mastos, trig, eeg_elec, baseline, win_of_interest, conditions, chan_dir]= get_OPTIONS(OPTIONS)

indir = OPTIONS.indir ;
hp = OPTIONS.hp;
lp = OPTIONS.lp;
mastos = OPTIONS.mastos;
trig = OPTIONS.trig;
eeg_elec = OPTIONS.eeg_elec;
baseline = OPTIONS.baseline;
win_of_interest = OPTIONS.win_of_interest;
conditions = OPTIONS.conditions;
chan_dir = OPTIONS.chan_dir;

end