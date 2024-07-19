function visualize_rejection_rates_FFR(subjects,OPTIONS)
% Load .csv with rejection information by participant to merge all rejection information in one single .csv. 
% This new table is then visualized
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

%Export trial rejection information for all participants
HF_rejected = zeros(size(subjects));

%Loop through subjects to extract number of trials rejected
for jj=1:length(subjects)
    %Get individual .csv file with rejected trials information
    csv_name= fullfile(OPTIONS.indir,subjects{jj},strcat(subjects{jj}, OPTIONS.suffix_csv, OPTIONS.param, '.csv'));

    %Open for reading
    temp= readtable(csv_name);
    HF_rejected(jj) = sum((temp.rejected==1)&strcmp(temp.condition,'HF'));
end

% Create table to store these information
rej_info_all = table(subjects,HF_rejected);

%Create histograms with percentages of trials rejected by condition
table_HF = sortrows(rej_info_all,2);
figure ;
h1 = histogram(table_HF.HF_rejected/5100*100); h1.NumBins = 10 ; hold on ; 
% Add legend
legend('%of HF rejected') ; 

%Create barplots with percentages of trials rejected by condition
figure;
bar(categorical(table_HF.subjects,table_HF.subjects),table_HF.HF_rejected/5100*100);
grid on;
title('Percentage of HF trial rejected');
legend('% HF rejected');
savefig('percentage_HF_trials_rej.fig');

end