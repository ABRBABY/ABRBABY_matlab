function [] = rejected_trials_visu_ERPs(rej_info_all,n_std,n_dev)
%% Script for rejection rates visualization (avoids running all ERPs sanity check)

%Sort data by number of std, dev1 and dev2
table_std = sortrows(rej_info_all,2);
table_dev1 = sortrows(rej_info_all,3);
table_dev2 = sortrows(rej_info_all,4);

%==========================================================================
%Create histograms with percentages of trials rejected by condition

figure ;
h1=  histogram(table_std.std_rejected/n_std*100); h1.NumBins = 10 ; hold on ; 
h2 = histogram( table_std.dev1_rejected/n_dev*100,10);h2.NumBins = 10 ;
h3 = histogram(table_std.dev2_rejected/n_dev*100,10) ; h3.NumBins = 10 ;

% h1.Normalization = 'probability' ; 
h1.BinWidth = 0.5 ; 
% h2.Normalization = 'probability' ; 
h2.BinWidth = 0.5 ; 
% h3.Normalization = 'probability' ; 
h3.BinWidth = 0.5 ; 

legend('STD','DEV1','DEV2') ; 

%==========================================================================
%Create barplots with percentages of trials rejected by condition

%Plot number of trial rejected sorted by number of STD rejected
figure;
bar(categorical(table_std.subjects, table_std.subjects),[table_std.std_rejected/n_std*100, table_std.dev1_rejected/n_dev*100, table_std.dev2_rejected/n_dev*100]);
grid on;
title('Percentage of trial rejected, sorted by number of STD rejected');
legend('% STD rejected','% DEV1 rejected','% DEV2 rejected');
savefig('percentage_trials_rej_bySTD.fig');

%Plot number of trial rejected sorted by number of DEV1 rejected 
figure;
bar(categorical(table_dev1.subjects, table_dev1.subjects),[table_dev1.std_rejected/n_std*100, table_dev1.dev1_rejected/n_dev*100, table_dev1.dev2_rejected/n_dev*100]);
grid on;
title('Percentage of trial rejected, sorted by number of DEV1 rejected');
legend('% STD rejected','% DEV1 rejected','% DEV2 rejected');
savefig('percentage_trials_rej_byDEV1.fig');

%Plot number of trial rejected sorted by number of DEV2 rejected 
figure;
bar(categorical(table_dev2.subjects, table_dev2.subjects),[table_dev2.std_rejected/n_std*100, table_dev2.dev1_rejected/n_dev*100, table_dev2.dev2_rejected/n_dev*100]);
grid on;
title('Percentage of trial rejected, sorted by number of DEV2 rejected');
legend('% STD rejected','% DEV1 rejected','% DEV2 rejected');
savefig('percentage_trials_rej_byDEV2.fig');

end
