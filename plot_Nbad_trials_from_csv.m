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
function [] = plot_Nbad_trials_from_csv(indir)

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

for iSubj=1:length(subjects) %for each subject
    
    T = readtable(fullfile(indir,subjects{iSubj}, strcat(subjects{iSubj},'_trials_description.txt'))); 

    nb_trials(iSubj,1) = table2array(sum(T(:,contains(T.Properties.VariableNames,'autorej_low_25_high_25_stepA1_stepB1'))==0));
    nb_trials(iSubj,2) = table2array(sum(T(:,contains(T.Properties.VariableNames,'autorej_low_30_high_30_stepA2_stepB2'))==0));
   
end

figure('Units','normalized','position',[0,0.5,1,0.8]) ; 

subplot(2,2,[1 2]) ; plot(nb_trials(:,1),'k+','MarkerSize',12,'linewidth',2); hold on ; plot(nb_trials(:,2),'m+','MarkerSize',12,'linewidth',2) ;  grid on ; title('Number trials rejected (autorej) with different parameters');
set(gca,'Fontsize',14); ylabel('Number of trial  rejected');

legend('autorej=low25-high25','autorej=low30-high30');

set(gca,'XTick',1:length(subjects), 'XTickLabels',strrep(subjects,'_','\_'),'XTickLabelRotation',45); 

subplot(2,2,3); 
[phH1,phH2] = hist([nb_trials(:,1) nb_trials(:,2)]);
bh = bar(phH2,phH1);
set(bh(1),'FaceColor','k');
set(bh(2),'FaceColor','m');
set(gca,'Fontsize',14); xlabel('Number of trials rejected');legend('autorej=low25-high25','autorej=low30-high30');
title('Histogram of rejected trials accross subjects');

subplot(2,2,4); 
plot(nb_trials(nb_trials(:,2)>3000,1),'k+','MarkerSize',12,'linewidth',2); hold on ; plot(nb_trials(nb_trials(:,2)>3000,2),'m+','MarkerSize',12,'linewidth',2) ;  grid on ; title('Number trials rejected (autorej) with different parameters');
set(gca,'Fontsize',14); ylabel('Number of trial  rejected');

legend('autorej=low25-high25','autorej=low30-high30');
set(gca,'XTick',1:sum(nb_trials(:,2)>3000), 'XTickLabels',strrep(subjects(nb_trials(:,2)>3000),'_','\_'),'XTickLabelRotation',45); title('Subject with more than 3000 trials rejected')

% subplot(2,2,3); hist(nb_trials(:,1),50) ; legend('autorej=low25-high25'); title('Histogram of rejected trials accross subjects'); set(gca,'Fontsize',14); xlabel('Number of trial  rejected');
% subplot(2,2,4); hist(nb_trials(:,2),50) ; legend('autorej=low30-high30');title('Histogram of rejected trials accross subjects'); set(gca,'Fontsize',14); xlabel('Number of trial  rejected');

