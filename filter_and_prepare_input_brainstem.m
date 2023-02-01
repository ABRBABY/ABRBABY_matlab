function [out_filenames] = filter_and_prepare_input_brainstem(ALLEEG, EEG, OPTIONS, overwrite, tube_length, propag_sound, BT_toolbox)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%%%%%%%%% TO UPDATE %%%% this comes from old code 

% This function mainly do : 
%Check if RERBT(number).set file exists
%Filter data
%Compute mean activity (FFR)
%Add tube delay to mean
%Export FFR data into .txt file
%Convert .txt file into .avg for BT_toolbox
%Export timepoints from last subject

% Get OPTIONS
[indir, hp, lp, RERBT]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Check if RERBT(number).set files exist for all subjects
for jj=1:length(subjects) 
    fil = dir(fullfile(indir,subjects{jj},strcat(subjects{jj},'_reref_epoched_FFR_RERBT',num2str(RERBT),'.set')));
    if isempty(fil) ; error('_reref_epoched_FFR_RERBT%s file does not exist for subject %s', num2str(RERBT),subjects{jj}); end
end

%Loop through subjects
for jj=1:length(subjects) 

     % Printout the id of the subject in console
    fprintf(strcat(subjects{jj}, '...\n'));
    
    filepath = fullfile(indir,subjects{jj},strcat(subjects{jj},'_reref_epoched_FFR_RERBT',num2str(RERBT),'.set'));
    [filepath,filename,ext] = fileparts(filepath) ;

    %Check if output file with selected parameters already exists
    [does_exist, count] = check_exist_set_params(filename, subjects{jj},OPTIONS) ; 

    if does_exist && overwrite == 0; continue; end
     
    EEG = pop_loadset(strcat(filename, ext),filepath) ;

    % Creates resulting filename
    out_filenames{jj} = fullfile(indir,subjects{jj}, strcat(subjects{jj}, '_RERBT',num2str(RERBT),'_filtered_FFR_F', num2str(count),'.set')) ; 

    % Skip if subject rerefe filtered_epochs already exist and we don't
    % want to overwrite
    if exist(out_filenames{jj},'file') && overwrite == 0; continue; end

    % Select bdf file in the folder
    %EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));
    
    %Filter data
    EEG  = pop_basicfilter(EEG,  1 , 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 ); % GUI: 11-Apr-2022 12:47:48
    %[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_reref_epoched_FFR'),'gui','off');

    %Extract mean activity (erp) and replace data
    abr = mean(EEG.data(1,:,:),3);
    EEG.data = abr;

    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    abr_shifted = circshift(abr,nsample_delay) ;
    
    % Create a custom history variable to keep track of OPTIONS 
    EEG.history_f = OPTIONS ;
    
    %% SAVE DATASET
    pop_newset(ALLEEG, EEG, 1, 'setname', strcat(filename,'_filtered_FFR'),'savenew', out_filenames{jj},'gui','off');
    
    %% Export ABR data into .txt file
    fname_out = fullfile(filepath,strcat(subjects{jj},'_RERBT', num2str(RERBT),'_F', num2str(count),'_abr_shifted_data_HF.txt')) ;
    fid = fopen(fname_out,'w');
    fprintf(fid,'%c\n',abr_shifted);
    fclose(fid);
    
    addpath(BT_toolbox);
    bt_txt2avg(fname_out, EEG.srate, EEG.history_rerbt.win_of_interest(1)*1000, EEG.history_rerbt.win_of_interest(2)*1000);
end

%% Export times (from any subject : just timepoints)
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
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
function [indir, hp, lp, RERBT]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
hp = OPTIONS.hp;
lp = OPTIONS.lp;
RERBT = OPTIONS.RERBT;

end

%--------------------------------------------------------------
% FUNCTION that check if that sets of param exist 
%--------------------------------------------------------------
function [does_exist, count] = check_exist_set_params(filename, subject, OPTIONS)

% Reads all folders that are in indir 
d = dir(fullfile(OPTIONS.indir,subject, strcat(subject,'RERBT',num2str(OPTIONS.RERBT),'_filtered_FFR_F*')));
% No file exists 
if isempty(d) ; does_exist=0 ; count =1 ; return ; end

for ff=1:length(d) 
    
    EEG = pop_loadset('filepath',fullfile(OPTIONS.indir,subject, d(ff).name),'loadmode','info') ;  
    
    % Check if the set of param correspond to the current file
    if isequal(EEG.history_f,OPTIONS)
        does_exist = 1 ; 
        tmp = regexp(d(ff).name,'F\d*','Match');
        tmp2 =  regexp(tmp,'\d*','Match');
        count = cell2mat(tmp2{:}); 
        return ; 
    end
    
end

% At this point the set of params does not exist and a new file needs to be
% crated (with a count increment)
does_exist = 0 ; 
count = length(d) +1;

end
