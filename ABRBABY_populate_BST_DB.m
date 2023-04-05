function [] = ABRBABY_populate_BST_DB(set_params)
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
INDIR = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data';

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Loop through all subjects
for jj=1:length(subjects) 

    % Call BST functions
    process_pipeline(INDIR, subjects{jj}, set_params)

end
 
%% Function the link data file to BST (Review Raw)
function [] = process_pipeline(INDIR, SubjectName, set_params)

Conditions = {'DEV1', 'DEV2', 'STD1', 'STD2'};

% Input files
sFiles = [];
sFilesAvg = [];

% Start a new report
bst_report('Start', sFiles);

% Creates a subject in database using the dfault anatomy and default
% channels
[~, ~] = db_add_subject(SubjectName, [], 1, 1);
panel_protocols('UpdateTree'); % Update the Protocol in GUI

% Loop through conditions 
for cc=1:length(Conditions)
    
   RawFile = dir(fullfile(INDIR,SubjectName,strcat(SubjectName,'_',Conditions{cc},'*_balanced_',set_params,'*.set'))) ; 
   RawFile = fullfile(RawFile.folder, RawFile.name);
   
    % Process: Import MEG/EEG: Existing epochs
    sFiles = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
        'subjectname',   SubjectName, ...
        'condition',     Conditions{cc}, ...
        'datafile',      {RawFile, 'EEG-EEGLAB'}, ...
        'iepochs',       [], ...
        'eventtypes',    '', ...
        'createcond',    0, ...
        'channelalign',  0, ...
        'usectfcomp',    0, ...
        'usessp',        0, ...
        'freq',          [], ...
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

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

end
end
