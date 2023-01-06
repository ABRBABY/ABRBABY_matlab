function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj, timepoints, labels] = extract_averages_DEV_STD(INDIR,subjects,STD_setfile)

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject

    %dimensions XXX.data input (for 1 subject) = channels x timepoints x trials
    %dimensions datasets (.set) output (average) =  subjects x channels x timepoints
    
    STD1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},STD_setfile{1})) ; 
    EEG_std1 = pop_loadset(STD1File) ; 
    STD1_subj(loopnum,:,:)  = mean(EEG_std1.data,3) ;

    STD2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},STD_setfile{2})) ; 
    EEG_std2 = pop_loadset(STD2File) ; 
    STD2_subj(loopnum,:,:)  = mean(EEG_std2.data,3) ;

    DEV1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV1.set')) ; 
    EEG_dev1 = pop_loadset(DEV1File) ; 
    DEV1_subj(loopnum,:,:)  = mean(EEG_dev1.data,3) ; 

    DEV2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV2.set')) ; 
    EEG_dev2 = pop_loadset(DEV2File) ; 
    DEV2_subj(loopnum,:,:)  = mean(EEG_dev2.data,3) ; 

end

timepoints = EEG_std1.times;
labels = {EEG_std1.chanlocs.labels};

end
