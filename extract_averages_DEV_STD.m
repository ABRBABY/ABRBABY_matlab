function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj, timepoints, labels, xlim] = extract_averages_DEV_STD(subjects,OPTIONS)

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject

    %dimensions XXX.data input (for 1 subject) = channels x timepoints x trials
    %dimensions datasets (.set) output (average) =  subjects x channels x timepoints
    
    
    [STD1_subj(loopnum,:,:), timepoints, labels, xlim]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*STD1*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 
    
%     [STD2_subj(loopnum,:,:),~,~,~]   = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*STD2*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 

    [STD2_subj(loopnum,:,:),~,~,~]   = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*STD1*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 

    [DEV1_subj(loopnum,:,:),~,~,~]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*DEV1*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 
    
    [DEV2_subj(loopnum,:,:),~,~,~]  = get_data(fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'*DEV2*','_',OPTIONS.balance_STD,'*',OPTIONS.params,'*.set'))) ; 

end

end


function [data, timepoints, labels, xlim] = get_data(datafile) 

    FileToLoad = dir(datafile) ; 
    if size(FileToLoad,1) > 1
        rman = find(contains({FileToLoad.name}, '_gd_avg')) ;
        FileToLoad = FileToLoad(rman,1) ;
    end
    if size(FileToLoad,1) > 1
        rman = find(contains({FileToLoad.name}, 'rman')) ;
        FileToLoad = FileToLoad(rman,1) ;
    end
    EEG = pop_loadset(FileToLoad.name,FileToLoad.folder) ; 
    data = mean(EEG.data,3) ;
    
    timepoints = EEG.times;
    labels = {EEG.chanlocs.labels};
    xlim = [EEG.xmin, EEG.xmax]*1000 ; 

end

