function [EEG] = detect_events_and_create_report(EEG, indir, trig_modality, trig)

% Inputs : 
% EEG : EEG structure created with EEGLAB
% indir : directory path where participants folders are stored (1 folder
% per recording)
% Trig modality : if a 'participantID_ergstim.txt' exist in participant
% folder, the function considers that the triggers are stored in a 'Erg'
% channels and need to be detected with threshold. Otherwise, the function 
% considers that the triggers were received via DB37 plug.
% trig : name of trigger channel if stored in Erg channel. Generally 'Erg'.

% Check trigger modality (Erg or DB37)
erg_filename = fullfile(EEG.filepath, strcat(EEG.setname, trig_modality)) ;
if exist(erg_filename, 'file')
    % Find TRIG electrodes indices by labels
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig));
    % Erg dat : Extract event from trigger channel (Erg1)
    EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);
    % Resolve bad event detection linked to trigger artefact and detect events to remove
    if size(EEG.event,2)>6000 && sum(contains(readlines(fullfile(indir, 'arte_stim_participants.txt')), EEG.setname))
        [idx_to_remove_trigg] = resolve_event_detection_HF_trigg_artefact(EEG) ;
        % Removes events identified above
        EEG.event(idx_to_remove_trigg) = [] ;  EEG.urevent(idx_to_remove_trigg) = [] ;
    else
        % Identifies outliers events (e.g. boundaries) or too close events
        idx_to_remove = [   find(diff([EEG.event.latency])<0.219*EEG.srate),... % minimum intretrial duration = 219 ms
            find(diff([EEG.event.latency])>1.5*EEG.srate) ];    % maximum intertrial duration = around 1500 m
        % Removes outliers events
        EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ;
    end
else
    % DB37 data : Remove Erg channel
    EEG = pop_select( EEG, 'nochannel', trig) ;

end

% Create  _trial_description.txt file that contains event info : label,
% acquired or not (= to keep or reject), and latency
events_table = readtable(strcat(fullfile(EEG.filepath,EEG.setname),'.txt'),'ReadVariableNames', 0) ;
conditions_description = events_table{:,1};
flag_description = ones(1,height(conditions_description)) ;
latencies_description =  [EEG.event.latency];

if size(latencies_description,2)==height(conditions_description)
    T = table(conditions_description, flag_description', latencies_description','VariableNames', {'condition', 'rejection_acq', 'latency'}) ;
    writetable(T, strcat(fullfile(EEG.filepath,EEG.setname),'_trials_description.txt'), 'WriteVariableNames', 1) ;
else
    latencies_description = zeros(1,height(conditions_description)) ;
    T = table(conditions_description, flag_description', latencies_description') ;
    writetable(T, strcat(fullfile(EEG.filepath,EEG.setname),'_trials_description.txt'), 'WriteVariableNames', 1) ;
    warning('Manual event correction needed for participant %s.', EEG.setname) ;
    open("manual_event_correction.m")
end
end