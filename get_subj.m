function [subjects] = get_subj(indir)
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads a .xslx file that contains list of subjects to process*
d = readtable(indir) ;
subjects = table2array(d) ;