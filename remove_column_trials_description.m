% custom_path = '/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/ABRBABY_cfrancois/data'; 
custom_path = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data';

indir = fullfile(custom_path,'DEVLANG_data');

% Select column to remove (0 = keep, 1 = remove):
remove_autorej_erp = 1;
remove_autorej_ffr = 1;
remove_manual_rejection = 1;

% Reads all folders that are in indir
d = dir(indir);
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..te

% selected_subj = {'DVL_013_T24','DVL_011_T10','DVL_044_T8'};
% selected_subj = {'DVL_003_T6', 'DVL_003_T18', 'DVL_003_T10'};
selected_subj = [];

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

   %Remove manual rejection column if selected
   if remove_autorej_ffr==1 && sum(contains(T1.Properties.VariableNames,'ERP_autorej'))>0
       autorej_erp_variableIndex = find(contains(T1.Properties.VariableNames,'ERP_autorej'));
       T1(:,autorej_erp_variableIndex) = [] ;
   end

   if remove_autorej_ffr==1 && sum(contains(T1.Properties.VariableNames,'FFR_autorej'))>0
       autorej_ffr_variableIndex = find(contains(T1.Properties.VariableNames,'FFR_autorej'));
       T1(:,autorej_ffr_variableIndex) = [] ;
   end

   if remove_manual_rejection==1 && sum(strcmp('Manual_rejection',T1.Properties.VariableNames))>0
       T1.Manual_rejection = [];
   end

   writetable(T1,fname_trial_desc,'WriteVariableNames', true);
end
