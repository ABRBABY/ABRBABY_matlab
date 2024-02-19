function [] = add_manual_bad_trial_detection(ALLEEG,OPTIONS)
% A.-Sophie Dubarry 2024
% Function whch reads a directory containing all files manually marked and
% update the trial_file_description from the indir accordingly

% Reads all folders that are in INDIR 
d = dir(OPTIONS.manualdir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Intialize a flag for bad trials 
flag = ones(6000,1);
                   
%Loop through subjects
for jj=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{jj}, '...\n'));
    
    % Update the (.set) file 
    fname_trial_desc = fullfile(OPTIONS.indir,subjects{jj},strcat(subjects{jj},'_trials_description.txt'));
    
    % Test if trial_description.txt exist 
    if ~exist(fname_trial_desc,'file')
        error('\nABRBABY --------- File _trials_description.txt does not exist. You must run the first part of the analysis');
    else 
        % Read _trials_description file 
        T1 = readtable(fname_trial_desc); 
        
        % Find if column exists skip  subject
        if sum(strcmp('Manual_rejection',T1.Properties.VariableNames))
             fprintf(sprintf('\nABRBABY --------- %s Process already done because Manual_rejection colum exist in trial_description \n',subjects{jj}));
             continue ; 
        end 
        
    end

    % Get the manual marked .set file
    fname= dir(fullfile(OPTIONS.manualdir,subjects{jj},'*_256.set'));

    % For all files detected for this subject
    for ff=1:length(fname) 

        % Reads manually marked data file (.set)
        EEG = pop_loadset(fname(ff).name, fname(ff).folder) ;
    
        % Reads corresponding initial file 
        fname_orig= fullfile(OPTIONS.indir,subjects{jj},strrep(fname(ff).name,'_256.set','*'));
        
        if ~isempty(dir(fname_orig))
            fname_orig = dir(fname_orig) ; 
            EEGorig = pop_loadset(fname_orig, fname(ff).folder) ;
        else
            fname_orig = strrep(fname_orig,'_RFE','_stepA');
            fname_orig = strrep(fname_orig,'_REJ','_stepB');
            fname_orig = strrep(fname_orig,'_thresh','');
            fname_orig = strrep(fname_orig,'STD1','STD*');
            
            [pathname,tmpname,~] = fileparts(fname_orig) ; 
            fname_orig = fullfile(pathname,strcat(subjects{jj},'*',extractAfter(tmpname,subjects{jj}),'.set'));
            if isempty(dir(fname_orig))
                fprintf(strcat('No original file found for ',fname(ff).name,'...SKIP \n'));
                continue
             end
            fname_orig = dir(fname_orig) ; 
            EEGorig = pop_loadset(fname_orig.name, fname_orig.folder) ;
            
       end
          
       idx_poprej = strfind(EEG.history,'pop_rejepoch('); 
       bad_trials = [];
       for ss=1:length(idx_poprej)       
            % Get indices of manually rejected trials
            bad_trials = cat(2,bad_trials,str2num(extractBefore(extractAfter(EEG.history(idx_poprej(ss):end),'pop_rejepoch( EEG, '),',0);')));
       
       end

        % Here the number of trials in the original file is different
        % than the corrected + nb trials
        if EEG.trials+length(bad_trials)~=EEGorig.trials 
            error('ERROR : %s discrepency between file used to mark bad trials and automatically generated file\n',subjects{jj});
        end 
        
        % Remove bad trials and save a new .set file
        EEGorig = pop_rejepoch( EEGorig, bad_trials ,0);
        
        conditions = {'STD','DEV1','DEV2'} ; 
        
        for cc=1:length(conditions)
            if contains(fname(ff).name,conditions{cc})
                
                idx_cond = contains(T1.condition,conditions{cc}) ; 
                idx_column = contains(T1.Properties.VariableNames,extractBefore(extractAfter(fname_orig.name,'_step'),'.set'));
                idx_autorej = T1{:,idx_column} ; 
            
                idx_before_manual_rej = find(idx_cond.*idx_autorej); 
                flag(idx_before_manual_rej(bad_trials)) = 0 ; 
               
            end
        end
        
        % Add history to the EEG file 
        EEGorig.history_stepManualRej.rejectedtrials = bad_trials ; 
        
        % Save datasets 
        pop_newset(ALLEEG, EEGorig, 1, 'setname',fname_orig.name,'savenew', fullfile(fname_orig.folder,fname_orig.name),'gui','off');
        
    end

    % Write a new column in the trial_description file
    add_flag_column_trials_description(fname_trial_desc, 'Manual_rejection',flag') ; 

end

