function [data_avg, vTimes] = compute_grand_average_allcond(subjects, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project


%(dimensions = subjects x channels x timepoints)
[dev1, dev2, std1, ~ , timepoints, labels, xlim] = extract_averages_DEV_STD(subjects,OPTIONS);

vTimes = timepoints ; 

% Compute MMN average across subjects and conditions
data_avg = squeeze(mean((dev1-std1) + (dev2-std1) /2,1) );

% Compute average across channels
data_avg = mean(data_avg(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 