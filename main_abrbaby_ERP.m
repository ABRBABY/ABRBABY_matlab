% ERPs Pipeline analysisfor ABRBABY
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%% ------------------- Set environment 
% Variables to enter manually before running the code
name = 'EH';
% name = 'ASD';

% This function sets custom path (either for Estelle or AnneSo)
[eeglab_path, biosig_installer_path, erplab_path,indir,plot_dir,~] = get_custom_path(name);

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% ------------------- Preprocess : filter, reref, epoch, set chan positions
OPTIONS_rfe.indir = indir ;
OPTIONS_rfe.hp = 1;                                         % high-pass (Hz) (APICE)
OPTIONS_rfe.lp = 30;                                        % low-pass (Hz) (APICE) 
OPTIONS_rfe.mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; 
OPTIONS_rfe.trig = {'Erg1'};                                % Ref and trigger channels 
OPTIONS_rfe.baseline = [-99, 0] ; 
OPTIONS_rfe.win_of_interest = [-0.1, 0.5] ; 
%OPTIONS_rfe.baseline = [-199, 0] ; 
%OPTIONS_rfe.win_of_interest = [-0.2, 0.5] ; 
OPTIONS_rfe.conditions = {'STD','DEV1','DEV2'} ; 
OPTIONS_rfe.eeg_elec = 1:16 ; 
OPTIONS_rfe.chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc') ; 
suffix_rfe = '_reref_filtered_epoched_RFE' ; 

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_rfe, count_rfe]= test_existance_of_params_in_db(OPTIONS_rfe, suffix_rfe) ; 

% If one element of the vector is 1 (=a subject to create)
if sum(flag_sub_to_create_rfe) ~= 0
    
    fname_txt = fullfile(OPTIONS_rfe.indir,strcat('PARAM',suffix_rfe,num2str(count_rfe),'.txt')) ;
   
    % If txt file does not exist -> creates it 
    if ~exist(fname_txt)
        writetable(struct2table(OPTIONS_rfe), fname_txt, 'Delimiter',',') ;
%         writestruct(OPTIONS_rfe, fname_txt,"FileType", 'auto') ;
%         tabletest = struct2table(OPTIONS_rfe);
%         table2 = rows2vars(tabletest);
%   
%         fileID = fopen(fname_txt,'w+');
%         fwrite(fileID,struct2table(OPTIONS_rfe));
%         fclose(fileID);
    end
    
end
    
% Reref filter epoch erp : only apply to subjects which were not already
% computed with this set of parameters (as defined by flag_sub_to_create) ;
[preproc_filenames] = reref_filter_epoch_erp(ALLEEG, OPTIONS_rfe,flag_sub_to_create_rfe, count_rfe, suffix_rfe) ;

%% ------------------- Preprocess : Select trials per condition and reject BAD trials 
OPTIONS_st.indir = indir ;                             % options st for 'select trials'
OPTIONS_st.rej_low = -150 ;                            % 150 infants; 120 adults
OPTIONS_st.rej_high = 150 ;                            % 150 infants; 120 adults
OPTIONS_st.RFE = 1;                                    % indicates index of rfe set of parameters to use
bloc = repelem(1:30,30) ;                              % creates a vector of [1 1 1 1 (30 times) 2 2 2 2 (30 times) etc. up to 30]
suffix_st = '_select_trials_ST' ;

% Test if this set of params exists and returns the files to process and
% counter to use to name the saved files
[flag_sub_to_create_st, count_st]= test_existance_of_params_in_db(OPTIONS_st, suffix_st) ; 

% % If one element of the vector is 1 (=a subject to create)
% if sum(flag_sub_to_create_st) ~= 0
%     
%     fname_txt = fullfile(OPTIONS_rst.indir,strcat('PARAM',suffix,num2str(count_st),'.txt')) ;
%    
%     % If txt file does not exist -> creates it 
%     if ~exist(fname_txt)
%         writetable(struct2table(OPTIONS_rfe), fname_txt, 'Delimiter',',') ;
%     end
%     
% end
    
% Reject bad trials and save new .set file
select_and_save_trials_per_condition(ALLEEG, OPTIONS_st, 'balanced', flag_sub_to_create_st, count_st, suffix_st) ; 
select_and_save_trials_per_condition(ALLEEG, OPTIONS_st, 'unbalanced', flag_sub_to_create_st, count_st, suffix_st) ; 

% Write csv file directly into the subject dir
[~] = reject_trials_produce_report(preproc_filenames, eeg_elec, bloc, win_of_interest, rej_low, rej_high,'') ; 

%% ------------------- Display : 
elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};
%plot_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder' ; 
display_timeseries_by_condition(preproc_filenames, elec_subset, 'balanced',plot_dir) ; 
