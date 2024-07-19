function [] = ABRBABY_populate_BST_DB_averages(set_params, opt_balance, INDIR)
% ========================================================================
% This file is part of ABRBABY project
% 
% Free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2023 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% ------------------------------------------------------------------------
%
% DESCRIPTION : 
%   Function aims to populate a brainstorm database from a directory
%   containging one folder per participant (name of the folder will be used
%   as the name of the subject). (Brainstorm must be up and running at exec)
%   The files must have been preprocessed with main_abrbay_ERP.m 
%    
%   INPUT : Each subject folder must include : 
%       * REJ_REF files : set_params
%   OUTPUT : a Brainstorm database 
%

%% Main 
% Input directory 
% INDIR = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data';
% INDIR = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data' ;
% INDIR = 'D:\preprocessed_data_EEG\RFE1_REJ1' ;
% INDIR = 'E:\EEG_ANALYSES\EEGdata_CF_revised_byparticipant_all' ; 

OPTIONS.indir = INDIR ;
OPTIONS.params = set_params ;
OPTIONS.opt_balance = opt_balance ;
OPTIONS.writecsv = 0 ;
conditionsMMN = { 'DEV1-STD', 'DEV2-STD'} ; 
vTime = [-0.200012207, 0.4999389648] ; 
sGroups{1} = ["_T6","_T8","_T10"] ; 
sGroups{2} = ["_T18","_T24"] ; 
 
% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects_list = {d(isub).name}';
subjects_list(ismember(subjects_list,{'.','..'})) = []; % Removes . and ..

% Choose subjects to process
% choose_subj = readtable('\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data\ffr_participants_ok.csv', 'Delimiter', 'Space', 'ReadVariableNames',0) ;
% subjects = table2array(choose_subj) ;
subjects = get_subjects(INDIR, '') ;           % get subjects in OPTIONS_rej.file


% % Remove subjects based on number of trial rejected 
% thresh = 0.33;
% subjects = filter_subjects_based_rejection(subjects, thresh, OPTIONS) ;

% Loop through all subjects
for jj=1:length(subjects) 

    % Process individual analysis : Feed the BST database with subjects data 
    process_pipeline(INDIR, subjects{jj}, set_params, opt_balance);

end

% Process contrast by groups : t-test
process_group_analysis(conditionsMMN, vTime, sGroups) ; 
 
%-----------------------------------------------------------------------------
% Process individual analysis (import, rej bad trials, average)
%-----------------------------------------------------------------------------
function [] = process_pipeline(INDIR, SubjectName, set_params,opt_balance)

Conditions = {'DEV1', 'DEV2', 'STD'};

% Input files
sFiles = [];
sFilesAvg = [];

% Start a new report
bst_report('Start', sFiles);

% Creates a subject in database using the default anatomy and default
% channels
[~, ~] = db_add_subject(strcat(set_params,'_',opt_balance, '_',SubjectName) , [], 1, 1);
panel_protocols('UpdateTree'); % Update the Protocol in GUI

% Loop through conditions 
for cc=1:length(Conditions)
    
    RawFile = dir(fullfile(INDIR,SubjectName,strcat(SubjectName,'_ERP_',Conditions{cc}, '*_',opt_balance, '_',set_params,'_gd_avg.set'))) ; 
   if size(RawFile,1) > 1
        rman = find(contains({RawFile.name}, 'rman')) ;
        RawFile = RawFile(rman,1) ;
    end
   RawFile = fullfile(RawFile.folder, RawFile.name);

   % Process: Import MEG/EEG: Time
    sFiles = bst_process('CallProcess', 'process_import_data_time', sFiles, [], ...
        'subjectname',   strcat(set_params,'_',opt_balance, '_', SubjectName) , ...
        'condition',     Conditions{cc}, ...
        'datafile',      {RawFile, 'EEG-EEGLAB'}, ...
        'timewindow',    [], ...
        'split',         0, ...
        'ignoreshort',   1, ...
        'channelalign',  1, ...
        'usectfcomp',    1, ...
        'usessp',        1, ...
        'freq',          [], ...
        'baseline',      [], ...
        'blsensortypes', 'EEG');
   
    % Process: Add EEG positions
    sFiles = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
        'channelfile', {'', ''}, ...
        'usedefault',  'ICBM152: BioSemi 16', ...  % ICBM152: BioSemi 16
        'fixunits',    1, ...
        'vox2ras',     0, ...
        'mrifile',     {'', ''}, ...
        'fiducials',   []);

    sFilesAvg = [sFilesAvg, sFiles] ; 

end

% Subtract DEV1-STD1
sFiles = bst_process('CallProcess', 'process_diff_ab', sFilesAvg(contains({sFilesAvg.FileName}, 'DEV1')), sFilesAvg(contains({sFilesAvg.FileName}, 'STD')));

% Subtract DEV2-STD2
sFiles = bst_process('CallProcess', 'process_diff_ab', sFilesAvg(contains({sFilesAvg.FileName}, 'DEV2')), sFilesAvg(contains({sFilesAvg.FileName}, 'STD')));

% % Save and display report
% ReportFile = bst_report('Save', sFiles);
% bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

end

%-----------------------------------------------------------------------------
% Process Group contrast 
%-----------------------------------------------------------------------------
function [] = process_group_analysis(conditionsMMN, vTime, sGroups)

cond{1}= bst_process('CallProcess', 'process_select_files_data', [], [], 'condition',     conditionsMMN{1}) ; 
cond{2}= bst_process('CallProcess', 'process_select_files_data', [], [], 'condition',     conditionsMMN{2}) ; 

% Process: Average: By folder (grand average)
sFiles = bst_process('CallProcess', 'process_average', [cond{1},cond{2}], [], ...
    'avgtype',       4, ...  % By folder (grand average)
    'avg_func',      1, ...  % Arithmetic average:  mean(x)
    'weighted',      0, ...
    'keepevents',    0);

for ccMMN=1:length(conditionsMMN)
   
    sFiles = cond{ccMMN} ; 
    
    sFilesG1 = sFiles(contains({sFiles.FileName},sGroups{1})); 
    sFilesG2 = sFiles(contains({sFiles.FileName},sGroups{2})); 
    
    % Process: t-test equal [-200ms,500ms]          H0:(A=B), H1:(A<>B)
    sFiles = bst_process('CallProcess', 'process_test_parametric2', sFilesG1, sFilesG2, ...
        'timewindow',    [], ...
        'sensortypes',   '', ...
        'isabs',         0, ...
        'avgtime',       0, ...
        'avgrow',        0, ...
        'Comment',       '', ...
        'test_type',     'ttest_equal', ...  % Student's t-test   (equal variance)        A,B~N(m,v)t = (mean(A)-mean(B)) / (Sx * sqrt(1/nA + 1/nB))Sx = sqrt(((nA-1)*var(A) + (nB-1)*var(B)) / (nA+nB-2)) df = nA + nB - 2
        'tail',          'two');  % Two-tailed
    
    % Process: Average: Everything
    sFilesAvgG1 = bst_process('CallProcess', 'process_average', sFilesG1, [], ...
        'avgtype',       1, ...  % Everything
        'avg_func',      1, ...  % Arithmetic average:  mean(x)
        'weighted',      0, ...
        'keepevents',    0);
    
    % Process: Add tag: Group1
    bst_process('CallProcess', 'process_add_tag', sFilesAvgG1, [], ...
        'tag',           'Group1', ...
        'output',        1);  % Add to file name
    
    % Process: Average: Everything
    sFilesAvgG2 = bst_process('CallProcess', 'process_average', sFilesG2, [], ...
        'avgtype',       1, ...  % Everything
        'avg_func',      1, ...  % Arithmetic average:  mean(x)
        'weighted',      0, ...
        'keepevents',    0);
    
    % Process: Add tag: Group2
    bst_process('CallProcess', 'process_add_tag', sFilesAvgG2, [], ...
        'tag',           'Group2', ...
        'output',        1);  % Add to file name

end

end


end
