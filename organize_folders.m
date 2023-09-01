function [] = organize_folders(indir, directory, files_num)

% Function that organize destination folder as a function of a source
% folder and moves files into corresponding new folders.
% indir = source folder, contains list of folders to reproduce in
% destination folder
% directory = destination folder that contains files to be organized.
% WARNING : files to be organized should contain the name of the folder in
% which it is going to be moved

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Reads folders that are in directory
inside = dir(directory) ;

% Create new folders in directory based on folder in indir
for ss = 1:length(subjects)
    files = contains({inside.name}, subjects{ss}) ;
     mkdir(fullfile(directory,subjects{ss})) ;
% Move files into corresponding folders    
    if sum(files)>0
        files_subj = inside(files) ;
        for ii = 1:length(files_subj)
            movefile(fullfile(directory,files_subj(ii).name),fullfile(directory,subjects{ss}))
        end
    end
end

% Check that folders contains right number of files; display text otherwise
isub = [inside(:).isdir]; % returns logical vector if is folder
subjects = {inside(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

for ss = 1:length(subjects)
    if length(dir(fullfile(directory,subjects{ss})))<(files_num+2)
        fprintf('File %s is not complete. %n', subjects{ss})
    end
end

end
