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
       
        tmp = conditions{cc}  ; 
        if contains(conditions{cc},'STD') ;  
               conditions{cc} = 'STD*'; 
        end 
        
        filename = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_',OPTIONS.analysis,'_',conditions{cc},'_',opt_balance,'_',param,'*.set'))) ; 
        outname = fullfile(indir,subjects{ii},strcat(subjects{ii},'_',OPTIONS.analysis,'_',tmp,'_',opt_balance,'_',param,'_',OPTIONS.keyword,'.set')) ; 

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
%         EEG = pop_resample( EEG, OPTIONS.srate);
       
        [~,outfname,~] = fileparts(outname);

        pop_newset(ALLEEG, EEG, 1, 'setname',outfname,'savenew', outname,'gui','off');
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
