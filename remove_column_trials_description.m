custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data'; 
% custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data');

% Reads all folders that are in indir
d = dir(indir);
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..te

% selected_subj = {'DVL_013_T24','DVL_011_T10','DVL_044_T8'};
selected_subj = {'DVL_003_T6'};
% selected_subj = [];

if ~isempty(selected_subj)
    % Get selected subjects
    subjects = subjects(ismember(subjects,selected_subj)) ; 
end
%Loop through subjects
for jj=1:length(subjects)

    % Update the (.set) file 
    fname_trial_desc = fullfile(indir,subjects{jj},strcat(subjects{jj},'_trials_description.txt'));
  
   % Read _trials_description file 
   T1 = readtable(fname_trial_desc);
   
    T1.Manual_rejection = [];

    writetable(T1,fname_trial_desc,'WriteVariableNames', true);

end
