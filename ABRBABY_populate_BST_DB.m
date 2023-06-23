function [] = ABRBABY_populate_BST_DB(set_params, opt_balance)
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
%INDIR = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data';
INDIR = '\\Filer\home\Invites\herve\Mes documents\These\EEG\Data\DEVLANG_data' ;
OPTIONS.indir = INDIR ;
OPTIONS.params = set_params ;
OPTIONS.opt_balance = opt_balance ;
OPTIONS.writecsv = 0 ;

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Remove subjects based on number of trial rejected 
thresh = 0.33;
subjects = filter_subjects_based_rejection(subjects, thresh, OPTIONS) ;

% Loop through all subjects
for jj=1:length(subjects) 

    % Call BST functions
    process_pipeline(INDIR, subjects{jj}, set_params, opt_balance)

end
 
%% Function the link data file to BST (Review Raw)
function [] = process_pipeline(INDIR, SubjectName, set_params,opt_balance)

Conditions = {'DEV1', 'DEV2', 'STD1', 'STD2'};

% Input files
sFiles = [];
sFilesAvg = [];

% Start a new report
bst_report('Start', sFiles);

% Creates a subject in database using the dfault anatomy and default
% channels
[~, ~] = db_add_subject(strcat(set_params,'_',opt_balance, '_',SubjectName) , [], 1, 1);
panel_protocols('UpdateTree'); % Update the Protocol in GUI

% Loop through conditions 
for cc=1:length(Conditions)
    
   RawFile = dir(fullfile(INDIR,SubjectName,strcat(SubjectName,'_',Conditions{cc},'*_', opt_balance, '_*',set_params,'*.set'))) ; 
   if size(RawFile,1) > 1
        rman = find(contains({RawFile.name}, 'rman')) ;
        RawFile = RawFile(rman,1) ;
    end
   RawFile = fullfile(RawFile.folder, RawFile.name);
   
    % Process: Import MEG/EEG: Existing epochs
    sFiles = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
        'subjectname',   strcat(set_params,'_',opt_balance, '_', SubjectName) , ...
        'condition',     Conditions{cc}, ...
        'datafile',      {RawFile, 'EEG-EEGLAB'}, ...
        'iepochs',       [], ...
        'eventtypes',    '', ...
        'createcond',    0, ...
        'channelalign',  0, ...
        'usectfcomp',    0, ...
        'usessp',        0, ...
        'freq',          256, ...
        'baseline',      [], ...
        'blsensortypes', 'MEG, EEG');
    
    % Process: Add EEG positions
    sFiles = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
        'channelfile', {'', ''}, ...
        'usedefault',  'ICBM152: BioSemi 16', ...  % ICBM152: BioSemi 16
        'fixunits',    1, ...
        'vox2ras',     0, ...
        'mrifile',     {'', ''}, ...
        'fiducials',   []);

    % Process: Average: By condition 
     avg_sFiles = bst_process(...
        'CallProcess', 'process_average', ...
        sFiles, [], ...
        'avgtype', 4, ...
        'avg_func', 1, ...  % <HTML>Arithmetic average: <FONT color="#777777">mean(x)</FONT>
        'keepevents', 0);
    
     sFilesAvg = [sFilesAvg, avg_sFiles] ; 

end

% Subtract DEV1-STD1
sFiles = bst_process('CallProcess', 'process_diff_ab', sFilesAvg(contains({sFilesAvg.FileName}, 'DEV1')), sFilesAvg(contains({sFilesAvg.FileName}, 'STD1')));

% Subtract DEV2-STD2
sFiles = bst_process('CallProcess', 'process_diff_ab', sFilesAvg(contains({sFilesAvg.FileName}, 'DEV2')), sFilesAvg(contains({sFilesAvg.FileName}, 'STD2')));

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

end

cond1= bst_process('CallProcess', 'process_select_files_data', [], [], 'condition',     'DEV1-STD1') ; 
cond2= bst_process('CallProcess', 'process_select_files_data', [], [], 'condition',     'DEV2-STD2') ; 

% Process: Average: By folder (grand average)
sFiles = bst_process('CallProcess', 'process_average', [cond1,cond2], [], ...
    'avgtype',       4, ...  % By folder (grand average)
    'avg_func',      1, ...  % Arithmetic average:  mean(x)
    'weighted',      0, ...
    'keepevents',    0);

conditionsMMN = { 'DEV1-STD1', 'DEV2-STD2'} ; 

for ccMMN=1:length(conditionsMMN)

    % Process: Select data files in: */DEV1-STD1/_T8
    sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   [], ...
        'condition',    conditionsMMN(ccMMN), ...
        'tag',           '', ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0);
    
    sFilesG1 = sFiles(contains({sFiles.FileName},["_T6","_T8","_T10"])); 
    sFilesG2 = sFiles(contains({sFiles.FileName},["_T18","_T24"])); 
    
    % Process: t-test equal [-200ms,500ms]          H0:(A=B), H1:(A<>B)
    sFiles = bst_process('CallProcess', 'process_test_parametric2', sFilesG1, sFilesG2, ...
        'timewindow',    [-0.200012207, 0.4999389648], ...
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
    
    % Process: Add tag: toto
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
