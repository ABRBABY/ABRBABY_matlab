function [lat, amp, auc] = search_for_mmn(subjects, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%(dimensions = subjects x channels x timepoints)
[dev1, dev2, std1, ~ , timepoints, labels, xlim] = extract_averages_DEV_STD(subjects,OPTIONS);

vTimes = timepoints ; 

% Compute MMN average across subjects and conditions
data_avg = squeeze(mean((dev1-std1) + (dev2-std1) /2,1) );

% Compute average across channels
data_avg = mean(data_avg(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 

[lat_gd, ~, ~] = search_for_local_peak(data_avg,vTimes,OPTIONS.win_gd_mmn,OPTIONS) ; 

% Get files to average according to options
for ii = 1:length(subjects)

    % Compute MMN for this participant 
    mmn1 = squeeze(dev1(ii,:,:) - std1(ii,:,:)) ; 
    mmn2 = squeeze(dev2(ii,:,:) - std1(ii,:,:)) ; 
    mmn1_subset = mean(mmn1(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 
    mmn2_subset = mean(mmn2(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 

    [lat(ii,1), amp(ii,1), auc(ii,1)] = search_for_local_peak(mmn1_subset,timepoints,lat_gd+OPTIONS.win_mmn,OPTIONS) ; 
    [lat(ii,2), amp(ii,2), auc(ii,2)] = search_for_local_peak(mmn2_subset,timepoints,lat_gd+OPTIONS.win_mmn,OPTIONS) ; 

end
 