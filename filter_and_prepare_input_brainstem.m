function [out_filenames] = filter_and_prepare_input_brainstem(ALLEEG, OPTIONS,tube_length, propag_sound,flag_sub_to_create, count,suffix)
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
[indir, hp, lp, RERBT, BT_toolbox]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Inititalize output parameter
out_filenames = [] ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

%Check if RERBT(number).set files exist for all subjects
% for jj=1:length(subjects)
%     setname = dir(fullfile(indir,subjects{jj},strcat(subjects{jj},'_reref_epoched_FFR_RERBT',num2str(RERBT),'.set')));
%     if isempty(setname) ; error('_reref_epoched_FFR_RERBT%s file does not exist for subject %s', num2str(RERBT),subjects{jj}); end
% end

%Loop through subjects
for ii=1:length(subjects) 

     % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    %Set rerbt file to work on
    file_rerbt = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_reref_epoched_FFR_RERBT',num2str(RERBT),'.set'))) ;
    
    %Get filepath
    filepath = file_rerbt.folder ;
    
    % Creates resulting filename
    out_filenames{ii} = fullfile(indir,subjects{ii}, strcat(subjects{ii},suffix,num2str(count),'.set')) ; 
    
    %Load the RERBT .set file to work on
    EEG = pop_loadset(strcat(subjects{ii},'_reref_epoched_FFR_RERBT',num2str(RERBT),'.set'), filepath) ; 
    
    %Filter data
    EEG  = pop_basicfilter(EEG,  1 , 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 ); % GUI: 11-Apr-2022 12:47:48
    
    %Extract mean activity (erp) and replace data
    abr = mean(EEG.data(1,:,:),3);
    EEG.data = abr;

    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    abr_shifted = circshift(abr,nsample_delay) ;
    
    % Create a custom history variable to keep track of OPTIONS 
    EEG.history_f = OPTIONS ;
    
    %% SAVE DATASET
    pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_filtered_FFR','_RERBT',num2str(RERBT),'_F',num2str(count))),'gui','off');

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
function [indir, hp, lp, RERBT, BT_toolbox]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
hp = OPTIONS.hp;
lp = OPTIONS.lp;
RERBT = OPTIONS.RERBT;
BT_toolbox = OPTIONS.bt_toolbox;
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
