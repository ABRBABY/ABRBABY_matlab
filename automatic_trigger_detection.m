% Automatic trigger detection

custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';
indir = fullfile(custom_path,'DEVLANG_data');

% OPTIONS_rfe.indir = indir ;
% OPTIONS_rfe.hp = 1;                                         % high-pass (Hz) (APICE)
% OPTIONS_rfe.lp = 30;                                        % low-pass (Hz) (APICE)
% OPTIONS_rfe.mastos = {'Lmon','Rmon','MASTOG','MASTOD'};
% OPTIONS_rfe.                                % Ref and trigger channels
% OPTIONS_rfe.baseline = [-99, 0] ;
% OPTIONS_rfe.win_of_interest = [-0.1, 0.5] ;
% OPTIONS_rfe.conditions = {'STD','DEV1','DEV2'} ;
% OPTIONS_rfe.eeg_elec = 1:16 ;
% OPTIONS_rfe.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ;
% OPTIONS_rfe.varhistory = 'EEG.history_rfe' ;
% OPTIONS_rfe.analysis = 'ERP';
% 
% % Get OPTIONS
% [indir, hp, lp, mastos, trig, eeg_elec, baseline, win_of_interest, conditions, chan_dir]= get_OPTIONS(OPTIONS_rfe) ;

trig = {'Erg1'};
TRIG_MODALITY = '_ergstim.txt';

% Reads all folders that are in indir
d = dir(indir);
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Inititalize output parameter
out_filenames = [] ;

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