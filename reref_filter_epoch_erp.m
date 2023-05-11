function [out_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS, flag_sub_to_create, count, suffix, interpol)
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
    out_filenames{jj} = fullfile(indir,subjects{jj}, strcat(filename,suffix,num2str(count),'.set')) ; 
            
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

    % Find REF electrodes indices by labels 
    ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 

    % Save a first dataset in EEGLAB 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');

   %% Interpolate bad channels if needed
   
   % Add channels information
   EEG = pop_chanedit(EEG, 'lookup',chan_dir) ;

    % Run interpolation
    if contains(subjects{jj},interpol.subj)
        % Find index of interpol that refers to the current subject
        ind = find(contains(interpol.subj,subjects{jj})) ;
        % Find electrode(s) to interpolate with index
        interp_elec = find(ismember({EEG.chanlocs.labels},interpol.chan{ind})) ;
        % Interpolate channel(s)
        EEG = pop_interp(EEG, interp_elec, 'spherical') ;
        % Add interpolation information in EEG structure
        EEG.interpolation = interpol.chan{ind} ;
    else 
        EEG.interpolation = 'none' ;
    end

   
    
    %% RE-REF (excluding trig channel)
    % Find TRIG electrodes indices by labels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 

    % Re-reference data and rename new file
    EEG = pop_reref(EEG, ref_elec, 'exclude',trigg_elec, 'keepref','on');

    %% EVENTS 
    % Extract event from trigger channel (Erg1)
    EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);

    % Identifies outliers events (e.g. boundaries) or too close events 
%     idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate),... % minimum intretrial duration = 220 ms
%                         find(diff([EEG.event.latency])>2*EEG.srate) ];  
%   % Identifies outliers events (e.g. boundaries) or too close events 
    idx_to_remove = [   find(diff([EEG.event.latency])<0.219*EEG.srate),... % minimum intretrial duration = 219 ms
                        find(diff([EEG.event.latency])>1.5*EEG.srate) ];    % maximum intertrial duration = around 1500 m
    
    % Removes outliers events
    EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ; 

    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;

    %% FILTERS the data with ERPLab
    EEG  = pop_basicfilter(EEG,  eeg_elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    % Extract ALL conditions epochs
    EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');

    % Remove baseline
    EEG = pop_rmbase( EEG, baseline,[] );

%     % Add channels information
%     EEG=pop_chanedit(EEG, 'lookup',chan_dir);

    % Create a custom history variable to keep trakc of OPTIONS 
    EEG.history_rfe = OPTIONS ;
    
    %% SAVE DATASET 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_reref_filtered_epoched'),'savenew', out_filenames{jj},'gui','off');

end
end


%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = read_custom_events(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 0);

% Insert info from .txt into EEG.event
my_events = table2array(my_events);

out_event = struct('latency', {in_event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {in_event(:).urevent});

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

