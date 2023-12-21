function [out_filenames] = reref_filter_epoch(ALLEEG, OPTIONS, flag_sub_to_create, count, suffix)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Get OPTIONS
[indir, hp, lp, mastos, trig, eeg_elec, baseline, win_of_interest, conditions, chan_dir]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Inititalize output parameter
out_filenames = [] ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

%Loop through subjects
for jj=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{jj}, '...\n'));

    %% IMPORT
    % Get BDF file
    %[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
    fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
    [~,filename,~] = fileparts(fname.name);    

    % Creates resulting filename
    out_filenames{jj} = fullfile(indir,subjects{jj}, strcat(filename,'_',OPTIONS.analysis,suffix,num2str(count),'.set'));
            
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

    % Save a first dataset in EEGLAB 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');

    % Compute ABR (to computed before reref)
    [abr_signal, id_left] = compute_ABR(EEG) ; 
    
     % Add channels information
    EEG = pop_chanedit(EEG, 'lookup',chan_dir) ;

    %% Interpolate bad channels if needed
    InterpolChannels = fullfile(indir,subjects{jj},strcat(subjects{jj},'_interp.txt')); 

    if exist(InterpolChannels,'file') 
        % Get the electrode(s) to re-reference
        interp_chans = cellstr(strrep(fileread(InterpolChannels),' ','')) ; 
        % Find electrode(s) to interpolate with index
        interp_elec = find(ismember({EEG.chanlocs.labels},split(interp_chans,','))) ; 
        % Interpolate channel(s)
        EEG = pop_interp(EEG, interp_elec, 'spherical') ;
        % Add interpolation information in EEG structure
        EEG.interpolation = interp_chans ;
    else 
        EEG.interpolation = 'none' ;
    end
 
    %% RE-REF (excluding trig channel)
    % Find TRIG electrodes indices by labels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 

    % Find REF electrodes indices by labels 
    ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 

    % Re-reference data and rename new file
    EEG = pop_reref(EEG, ref_elec, 'exclude',trigg_elec, 'keepref','on');
    
    % Detect events when first time run on participant
    EventDetection = fullfile(indir,subjects{jj},strcat(subjects{jj},'_trials_description.txt')); 
    if ~exist(EventDetection,'file') 
        error('ABRBABY: File "_trials_description.txt" does not exist for participant %s.', subjects{jj})
    end

    % Replace Left channel by a new one called ABR
    EEG.data(id_left,:) = abr_signal ; EEG.chanlocs(id_left).labels = 'ABR' ;

    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','_trials_description.txt'),EEG.event) ;
    EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;

    %% FILTERS the data with ERPLab
    EEG  = pop_basicfilter(EEG,  eeg_elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );

    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); EEG = eeg_checkset( EEG );

    %% Extract ALL conditions epochs
    EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');

    %% Remove baseline
    EEG = pop_rmbase( EEG, baseline,[] );

    %% Create a custom history variable to keep trakc of OPTIONS 
    EEG.history_rfe = OPTIONS ;
    
    %% SAVE DATASET 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_reref_filtered_epoched_',OPTIONS.analysis),'savenew', out_filenames{jj},'gui','off');

end
end


%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = read_custom_events(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 0);

if size(my_events,2)~=3 
    error('Wrong number of column in file _trial_descriptions.txt');
else
    idx_events = my_events{:,2}==1 ; 
    out_event = struct('latency', num2cell(my_events{idx_events,3}'), ...
                    'type', my_events{idx_events,1}',...
                    'urevent', num2cell(1:height(my_events))) ;                      
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

%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [abr_signal,id_left]= compute_ABR(EEG)

    %% Compute FFR
    % forumla : {(Left)+(Right)}/-2 = Ref - {(LA+RA)/2}
    id_left  = find(strcmp({EEG.chanlocs.labels},'Left')) ;
    id_right = find(strcmp({EEG.chanlocs.labels},'Right')) ;    
    abr_signal = (EEG.data(id_left,:) + EEG.data(id_right,:))/-2 ;
end

