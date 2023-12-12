function[] = export_rejection_infos_FFR(subjects, OPTIONS)

%% Rejection rates : export into .csv
% Get rejection option name
part_param = split(OPTIONS.params, '_') ;
rejoption = part_param{2} ;
% Initialize variables
total_num_of_trials = [] ;
num_of_trials_rejected = [] ;
percent_of_rej = [] ;
n_trials_kept = [] ;

% Loop through subjects
for loopnum = 1:length(subjects) %for each subject
%for loopnum=find(ismember(subjects,'DVL_003_T10')) ; 
    rej_file = dir(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'_infos_trials_*',rejoption,'.csv'))) ; 
    filepath = fullfile(rej_file.folder, rej_file.name) ; 
    if exist(filepath,'file')==0 ; error(['File does not exist, please run FFR preprocessing steps on ',subjects{loopnum}]); end
    info_rej = readtable(filepath) ;
    total_num_of_trials(loopnum) = size(info_rej,1) ;
    num_of_trials_rejected(loopnum) = sum(info_rej.rejected) ;
    percent_of_rej(loopnum) = num_of_trials_rejected(loopnum)/total_num_of_trials(loopnum)*100 ;
    n_trials_kept(loopnum) = total_num_of_trials(loopnum) -  num_of_trials_rejected(loopnum) ;
end

rej_table = table(subjects,total_num_of_trials',num_of_trials_rejected',percent_of_rej', n_trials_kept','VariableNames', {'subjects', 'n_trial_total', 'n_trials_rejected', 'percent_trials_rejected', 'n_trials_kept'}) ;
writetable(rej_table, fullfile(OPTIONS.indir, strcat('FFR_rejection_info_all_subj_', OPTIONS.params,'.csv'))) ;


disp('Rejection info exported.')