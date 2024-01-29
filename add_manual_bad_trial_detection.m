function [] = add_manual_bad_trial_detection(OPTIONS)
% A.-Sophie Dubarry 2024
% Function whch reads a directory containing all files manually marked and
% update the trial_file_description from the indir accordingly

% Reads all folders that are in INDIR 
d = dir(OPTIONS.manualdir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Loop through subjects
for jj=1:length(subjects) 

    % Printout the id of the subject in console
    fprintf(strcat(subjects{jj}, '...\n'));

    % Get the manual marked .set file
    fname= dir(fullfile(OPTIONS.manualdir,subjects{jj},'*_256.set'));

    % For all files detected for this subject
    for ff=1:length(fname) 

        % Reads manually marked data file (.set)
        EEG = pop_loadset(fname(ff).name, fname(ff).folder) ;
            
        % Get indices of manually rejected trials
        bad_trials = str2num(extractBefore(extractAfter(EEG.history,'pop_rejepoch( EEG, '),' ,0);'));
        
        % Reads the automatically generated file (.set) 
        fname_orig = fullfile(OPTIONS.indir,subjects{jj},fname(ff).name);

        % Update the (.set) file 
        
        % fname_trial_desc = 

        % Reads trial_description file 
        
        % Update trial_description file


    end


    [~,filename,~] = fileparts(fname.name);    


end
