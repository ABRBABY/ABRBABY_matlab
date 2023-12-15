function [out_filenames] = reject_bad_trials_manual(ALLEEG, OPTIONS, opt_balance, subjects)

% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project
%INPUTS:
%- ALLEEG = EEGLAB whole structure
%- OPTIONS = rbt options
%- opt_balance = 'balanced' for balancing number of standard with number of
%deviants, or 'unbalanced' to keep all standards that are not rejected
%after trial rejection

%Get options
[indir, RFE_num,REJ_num]= get_OPTIONS(OPTIONS) ;
% cond = {'DEV1', 'DEV2'} ;

% Inititalize output parameter
out_filenames = [] ; 

if strcmp(opt_balance, 'balanced') == 1

    cond = {'DEV1', 'DEV2', 'STD1', 'STD2'} ;

% Loop through subjects
for ii=1:length(subjects)
    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Loop through conditions
    for cc = 1:length(cond)
        % Check that 'rejman' .txt file exist for each condition : if yes,
        % proceed to manual rejection
        if  exist(fullfile(indir,subjects{ii}, strcat(subjects{ii},'_rejman_',cond{cc},'.txt')), 'file') == 2
            if exist(fullfile(indir,subjects{ii}, strcat(subjects{ii},'_', cond{cc}, '_thresh_', opt_balance, RFE_num, REJ_num, '_rman.set')), 'file') == 0
            %Get txt and set filenames
            txtname = fullfile(indir,subjects{ii}, strcat(subjects{ii},'_rejman_',cond{cc},'.txt')) ;
            setname = strcat(subjects{ii},'_', cond{cc}, '_thresh_', opt_balance, RFE_num, REJ_num, '.set') ;
            % Read txt file to extract indices of trials to reject
            events_detected = readtext(txtname) ;
            events_to_rej = cell2mat(events_detected) ;
            %Read .set
            EEG = pop_loadset(setname,fullfile(indir, subjects{ii})) ;
            % Remove bad trials and save a new .set file
            EEG = pop_rejepoch( EEG, events_to_rej ,0);
            pop_newset(ALLEEG, EEG, 1,'setname', strcat(subjects{ii}, '_', cond{cc}, '_rman'),'savenew',fullfile(indir, subjects{ii}, strrep(setname, '.set', '_rman')),'gui','off');         
            end
        end
    end

end

end

if strcmp(opt_balance, 'unbalanced') == 1

cond = {'DEV1', 'DEV2', 'STDD'} ;

% Loop through subjects
for ii=1:length(subjects)
    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));
    
    % Loop through conditions
    for cc = 1:length(cond)
        % Check that 'rejman' .txt file exist for each condition : if yes,
        % proceed to manual rejection
        if  exist(fullfile(indir,subjects{ii}, strcat(subjects{ii},'_rejman_',cond{cc},'.txt')), 'file') == 2
            if exist(fullfile(indir,subjects{ii}, strcat(subjects{ii},'_', cond{cc}, '_thresh_', opt_balance, RFE_num, REJ_num, '_rman.set')), 'file') == 0
            %Get txt and set filenames
            txtname = fullfile(indir,subjects{ii}, strcat(subjects{ii},'_rejman_',cond{cc},'.txt')) ;
            if cc==3
                setname = strcat(subjects{ii},'_STD1', '_thresh_', opt_balance, RFE_num, REJ_num, '.set') ;
            else
                setname = strcat(subjects{ii},'_', cond{cc}, '_thresh_', opt_balance, RFE_num, REJ_num, '.set') ;
            end
            % Read txt file to extract indices of trials to reject
            events_detected = readtext(txtname) ;
            events_to_rej = cell2mat(events_detected) ;
            %Read .set
            EEG = pop_loadset(setname,fullfile(indir, subjects{ii})) ;
            % Remove bad trials and save a new .set file
            EEG = pop_rejepoch( EEG, events_to_rej ,0);
            pop_newset(ALLEEG, EEG, 1,'setname', strcat(subjects{ii}, '_', cond{cc}, '_rman'),'savenew',fullfile(indir, subjects{ii}, strrep(setname, '.set', '_rman')),'gui','off');         
            end
        end
    end

end

end

end

%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir,RFE_num,REJ_num]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;                     
%suffix_rman = OPTIONS.suffix_rman ;
RFE_num = OPTIONS.RFE_num ;
REJ_num = OPTIONS.REJ_num ;

end

