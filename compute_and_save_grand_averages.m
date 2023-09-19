function [] = compute_and_save_grand_averages(ALLEEG,OPTIONS)

% Get options
[indir, opt_balance, param, conditions]= get_OPTIONS(OPTIONS) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Get files to average according to options
for ii = 1:length(subjects)
    for cc = 1:length(conditions)
        filename = dir(fullfile(indir, subjects{ii},strcat(subjects{ii},'_', conditions{cc},'_thresh_',opt_balance,'_',param,'*.set'))) ;
        if size(filename,1) > 1
            rman = contains({filename.name},'rman') ;
            filename = filename(rman,1) ;
        end
        fname = filename.name ;
        filepath = filename.folder ;
        EEG = pop_loadset(fname, filepath) ;
        gd_avg = mean(EEG.data,3) ;
        EEG.data = gd_avg ;
        % Resample data to 256 Hz
        EEG = eeg_checkset( EEG );
        EEG = pop_resample( EEG, 256);
        pop_newset(ALLEEG, EEG, 1, 'setname',strrep(fname,'thresh','gd_avg'),'savenew', fullfile(filepath, strrep(fname,'thresh','gd_avg')),'gui','off');
    end
end

%--------------------------------------------------------------
% FUNCTION that get OPTIONS values
%--------------------------------------------------------------
function [indir, opt_balance, param, conditions]= get_OPTIONS(OPTIONS) 

indir = OPTIONS.indir ;
opt_balance = OPTIONS.opt_balance ;
conditions = OPTIONS.conditions ;
param = OPTIONS.param ;

end

end
