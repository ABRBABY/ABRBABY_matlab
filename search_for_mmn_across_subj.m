function [lat, amp, auc] = search_for_mmn_across_subj(subjects, timeW, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% Get files to average according to options
for ii = 1:length(subjects)

    %(dimensions = subjects x channels x timepoints)
    [dev1, dev2, std1, ~ , timepoints, labels, xlim] = extract_averages_DEV_STD(subjects(ii),OPTIONS);

    % Compute MMN for this participant 
    mmn1 = squeeze(dev1 - std1) ; 
    mmn2 = squeeze(dev2 - std1) ; 
    mmn1_subset = mean(mmn1(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 
    mmn2_subset = mean(mmn2(ismember(labels,OPTIONS.elec_subset(:)'),:),1) ; 

    [lat(ii,1), amp(ii,1), auc(ii,1)] = search_for_local_peak(mmn1_subset,timepoints,timeW ,OPTIONS) ; 
    [lat(ii,2), amp(ii,2), auc(ii,2)] = search_for_local_peak(mmn2_subset,timepoints,timeW ,OPTIONS) ; 

end
 