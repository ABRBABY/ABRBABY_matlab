function [out_filenames] = reref_epoch_ffr(ALLEEG, OPTIONS,flag_sub_to_create, count,suffix,bloc)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% This function mainly do : 
% Select channels
% Compute FFR and create ABR channel
% Detect trigg
% Reject outliers events (flaw with ergo channel)
% Map events (with labels) 
% Epoch

% Get OPTIONS
[indir, mastos, trig, abr, eeg_elec, baseline, win_of_interest, chan_dir,rej_low,rej_high]= get_OPTIONS(OPTIONS) ;

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
    %RERBT = Reref, Epoch, Reject Bad Trial
    out_filenames{jj} = fullfile(indir,subjects{jj}, strcat(filename,suffix,num2str(count),'.set')) ; 
        
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

    % Find REF electrodes indices by labels 
    ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 

    % Save a first dataset in EEGLAB 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');

    %% RE-REF (excluding trig channel)
    % Find TRIG electrodes indices by labels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 
    abr_elec = find(ismember({EEG.chanlocs.labels},abr)); 
    
    % Re-reference data and rename new file
    EEG = pop_reref(EEG, ref_elec, 'exclude',[trigg_elec,abr_elec], 'keepref','on');
    
    %% Compute FFR
    % forumla : {(Left)+(Right)}/-2 = Ref - {(LA+RA)/2}
    id_left  = find(strcmp({EEG.chanlocs.labels},'Left')) ; 
    id_right = find(strcmp({EEG.chanlocs.labels},'Right')) ; 
   
    abr_signal = (EEG.data(id_left,:) + EEG.data(id_right,:))/-2 ; 
    
    % Replace Left channel by a new one called ABR 
    EEG.data(id_left,:) = abr_signal ; EEG.chanlocs(id_left).labels = 'ABR' ; 
    
    %% Select Channels
    EEG = pop_select( EEG, 'channel',{'ABR', 'Erg1', EEG.chanlocs(eeg_elec).labels});
    
    % Re-index the trigger channel after removing some channels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 
   
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

    % Extract epochs for HF
    EEG = pop_epoch( EEG, {  'HF'  }, win_of_interest, 'newname', 'epochs', 'epochinfo', 'yes');
    
    %Remove baseline
    EEG = pop_rmbase( EEG, baseline,[]);

    % Add channels information
    EEG=pop_chanedit(EEG, 'lookup',chan_dir);
    
     % Create a custom history variable to keep trakc of OPTIONS 
    EEG.history_rerbt = OPTIONS ;
 
    %% SAVE DATASET 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_reref_epoched_FFR'),'savenew', out_filenames{jj},'gui','off');

%     % Reject bad trials and write a report
%     EEG = reject_trials_produce_report(out_filenames(jj), find(ismember({EEG.chanlocs.labels},'ABR')), bloc, win_of_interest, rej_low, rej_high,'FFR') ; 

    %% SAVE DATASET 
    %[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG{1}, CURRENTSET, 'setname', strcat(filename,'_reref_epoched_FFR'),'savenew', out_filenames{jj},'gui','off');

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
function [indir, mastos, trig, abr, eeg_elec, baseline, win_of_interest, chan_dir, rej_low, rej_high]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
mastos = OPTIONS.mastos;
trig = OPTIONS.trig;
abr = OPTIONS.abr;
eeg_elec = OPTIONS.eeg_elec;
baseline = OPTIONS.baseline;
win_of_interest = OPTIONS.win_of_interest;
chan_dir = OPTIONS.chan_dir;
rej_low = OPTIONS.rej_low;
rej_high = OPTIONS.rej_high;

end

% %--------------------------------------------------------------
% % FUNCTION that check if that sets of param exist 
% %--------------------------------------------------------------
% function [does_exist, count] = check_exist_set_params(filename, subject, OPTIONS)
% 
% % Reads all folders that are in indir 
% d = dir(fullfile(OPTIONS.indir,subject, strcat(filename,'_reref_epoched_FFR_RERBT*.set')) ); 
% 
% % No file exists 
% if isempty(d) ; does_exist=0 ; count =1 ; return ; end
% 
% for ff=1:length(d) 
%     
%     EEG = pop_loadset('filepath',fullfile(OPTIONS.indir,subject, d(ff).name),'loadmode','info') ;  
%     
%     % Check if the set of param correspond to the current file
%     if isequal(EEG.history_rerbt,OPTIONS)
%         does_exist = 1 ; 
%         tmp = regexp(d(ff).name,'RERBT\d*','Match');
%         tmp2 =  regexp(tmp,'\d*','Match');
%         count = cell2mat(tmp2{:}); 
%         return ; 
%     end
%     
% end
% 
% % At this point the set of params does not exist and a new file needs to be
% % crated (with a count increment)
% does_exist = 0 ; 
% count = length(d) +1;
% 
% end