function [out_filenames] = reref_filter_epoch(ALLEEG, indir, hp, lp, mastos, trig, eeg_elec, baseline, win_of_interest, conditions, chan_dir, overwrite)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

    % Reads all folders that are in indir 
    d = dir(indir); 
    isub = [d(:).isdir]; % returns logical vector if is folder
    subjects = {d(isub).name}';
    subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    
    %Loop through subjects
    for jj=1:length(subjects) 
           
        
        fprintf(strcat(subjects{jj}, '...\n'));
         
        %% IMPORT
        % Get BDF file
        %[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
        fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
        [~,filename,~] = fileparts(fname.name);    

        % Creates resulting filename
        out_filenames{jj} = fullfile(indir,subjects{jj}, strcat(filename,'_reref_filtered_epoched.set')) ; 
 
        % Skip if subject rerefe filtered_epochs already exist and we don't
        % want to overwrite
        if exist(out_filenames{jj},'file') && overwrite == 0; continue; end

        % Select bdf file in the folder
        EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));
    
        % Find REF electrodes indices by labels 
        ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 
        
        % Save a first dataset in EEGLAB 
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');
     
        %% RE-REF (excluding trig channel)
        % Find TRIG electrodes indices by labels 
        trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 
    
        % Re-reference data and rename new file
        EEG = pop_reref(EEG, ref_elec, 'exclude',trigg_elec, 'keepref','on');
        
        %% EVENTS 
        % Extract event from trigger channel (Erg1)
        EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);
        
        % Identifies outliers events (e.g. boundaries) or too close events 
        idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate),... % minimum intretrial duration = 220 ms
                            find(diff([EEG.event.latency])>2*EEG.srate) ];  
        % Removes outliers events
        EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ; 
        
        % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
        EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
        EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;
        
        %% FILTERS the data with ERPLab
        EEG  = pop_basicfilter(EEG,  eeg_elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset( EEG );
    
    %     %% SAVE DATASET BEFORE EPOCHING
    %     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_filtered'),'savenew', fullfile(indir,subjects{jj}, strcat(filename,'_filtered')),'gui','off');
    %     
        % Extract ALL conditions epochs
        EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');
    
        % Remove baseline
        EEG = pop_rmbase( EEG, baseline,[] );
        
        % Add channels information
        EEG=pop_chanedit(EEG, 'lookup',chan_dir);
        
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