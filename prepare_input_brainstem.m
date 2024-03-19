function [out_filenames] = prepare_input_brainstem(ALLEEG, OPTIONS,tube_length, propag_sound,flag_sub_to_create, count,suffix, RFE)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%%%%%%%%% TO UPDATE %%%% this comes from old code 

% This function mainly do : 
% Reject bad epoch based on window of rejection in OPTIONS
% Extract positive and negative FFR separately
% Compute mean activity (FFR) : 
% - mean added FFR ((positive+negative) /2)
% - mean subtracted FFR ((negative-positive) /2)
% - mean positive only
% - mean negative only
% Add tube delay to means
% Export FFR data into .txt file
% Convert .txt file into .avg for BT_toolbox
% Export timepoints from last subject

% Get OPTIONS
[indir, rej_low, rej_high, BT_toolbox, bloc, win_of_interest]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Inititalize output parameter
out_filenames = [] ; 

suffix_stepA = strrep(RFE,'_','') ; 

% Only keeps subjects to process
subjects = subjects(flag_sub_to_create) ; 

%Loop through subjects
for ii=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Input filename
    fname = strcat(subjects{ii},'_',OPTIONS.analysis,'_',suffix_stepA,suffix,num2str(count),'.set'); %dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_FFR_stepA',num2str(stepA),'.set'))) ;
    
     % Error if rfe file does not exist
    if ~exist(fullfile(OPTIONS.indir, subjects{ii},fname),'file') ; error('File %s does not exist for subject %s', fname, subjects{ii}); end
  
    %Load the stepA .set file to work on
    EEG = pop_loadset(fname, fullfile(OPTIONS.indir, subjects{ii})) ; 
    
    % Select only ABR elec
    EEG = pop_select(EEG, 'channel',{'ABR'});
    
    % Add tube delay (27 cm x 340 m/s ) 
    nsample_delay = fix(EEG.srate * (tube_length / propag_sound) ) ; 

    % Get and squeeze data 
    abr_trials = squeeze(EEG.data) ; 

    % Shift
    abr_trials = circshift(abr_trials,nsample_delay,1);

    % Computes the average of ABR signal 
    abr_average = mean(abr_trials,2) ; 
    
    % Creates a vector of alternating flags
    idx_flip1 = ones(1,size(EEG.data,3));
    idx_flip1(2:2:size(EEG.data,3))=0; 
    abr_average1 = mean(abr_trials(:,idx_flip1==1),2);
    df1 = sum(idx_flip1);

    idx_flip2 = ones(1,size(EEG.data,3));
    idx_flip2(1:2:size(EEG.data,3))=0; 
    abr_average2 = mean(abr_trials(:,idx_flip2==1),2);
    df2 = sum(idx_flip2);

    % Compute mean FFR by subtracting 
    abr_subtracted = (abr_average2 - abr_average1) /2 ;
    
    abr_types = {'avg','sub','neg','pos'};
    abr = cat(2,abr_average,abr_subtracted,abr_average2,abr_average1)' ; 

    for ff=1:length(abr_types) 
        fname_out = fullfile(OPTIONS.indir,strcat(subjects{ii},'_',num2str(suffix_stepA),suffix, num2str(count),'_abr_',abr_types{ff},'_shifted_data_HF.txt')) ;
        fid = fopen(fname_out,'w');
        fprintf(fid,'%c\n',abr(ff,:));
        fclose(fid);
       
        % Converts the output file into BT_Toolbpx compatible data 
        bt_txt2avg(fname_out, EEG.srate, EEG.history_stepA.win_of_interest(1)*1000, EEG.history_stepA.win_of_interest(2)*1000);
  
    end
    % % SAVE DATASETS
    % EEG.data = abr_averaged_add_shifted ;
    % EEG.averaged_sub = abr_average ; 
    % EEG.positive_pol = abr_average1;
    % EEG.negative_pol = abr_average2;
    % pop_newset(ALLEEG, EEG, 1, 'setname', strcat(subjects{ii},'_filtered_FFR'),'savenew', fullfile(filepath, strcat(subjects{ii},'_stepA',num2str(stepA),suffix,num2str(count))),'gui','off');

end

%% Export times (from any subject : just timepoints)
fname_out = fullfile(indir,'ABR_timepoints.txt') ;
fid = fopen(fname_out,'w');
fprintf(fid,'%f\n',EEG.times);
fclose(fid);
end


%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir, rej_low, rej_high, BT_toolbox, bloc, win_of_interest]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
rej_low = OPTIONS.rej_low ;                           
rej_high = OPTIONS.rej_high ; 
%varhistory = OPTIONS.varhistory ;
bloc = OPTIONS.bloc ;
BT_toolbox = OPTIONS.bt_toolbox;
win_of_interest = OPTIONS.win_of_interest ;
end
