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
    
    % Printout the id of the subject in console
    fprintf(strcat(subjects{ii}, '...\n'));

    for cc = 1:length(conditions)
       
        % Trick to handles either STD1, STD2 or STDD
        if contains(conditions{cc},'STD') 
               conditions{cc} = 'STD*'; 
        end 
        
        % Sets the filename on which we compute grand average
        filename = dir(fullfile(indir,subjects{ii},strcat(subjects{ii},'_',OPTIONS.analysis,'_',conditions{cc},'_',opt_balance,'_',param,'*.set'))) ; 
           
        % In case the grand average for this subject and this condition
        % already exists
        if size(filename,1) >1
            if contains([filename.name], OPTIONS.keyword)
                %Skip file if gd avg already done
                continue
            else
                error('More than one .set file for participant %s, condition %s.', subjects{ii}, conditions{cc})
            end
        end
        
        % Sets the resulting filename (with grand average)
        outname = fullfile(filename.folder, strrep(filename.name, '.set', strcat('_', OPTIONS.keyword, '.set'))) ;
     
        % Loads data to average
        EEG = pop_loadset(filename.name, filename.folder) ;
        
        % Computes one subject/one condition average 
        gd_avg = mean(EEG.data,3) ;
        EEG.data = gd_avg ;

        % Resample data to 256 Hz if needed
        EEG = eeg_checkset( EEG );
        if EEG.srate ~=256
            EEG = pop_resample( EEG, OPTIONS.srate);
        end
       
        [~,outfname,~] = fileparts(outname);

        % Save new EEGLAB dataset
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
