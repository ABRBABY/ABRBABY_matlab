function [] = FFR_analysis(subjects, OPTIONS) 
% ERPs analysis script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

grpA = OPTIONS.groups{1} ;
grpB = OPTIONS.groups{2} ;

%Color for plot
FFR_color = [0 0 0]; %black
grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.8902 0 0]; %red

% Get timepoints
FFR_times_file = fullfile(OPTIONS.indir,'ABR_timepoints.txt') ; 
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects) %for each subject
%for loopnum=find(ismember(subjects,'DVL_003_T10')) ; 
    FFR_file = fullfile(OPTIONS.indir,subjects{loopnum},strcat(subjects{loopnum},'_',OPTIONS.params,'_abr_',OPTIONS.ffr_polarity,'_shifted_data_HF.txt')) ; 
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR preprocessing steps on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

 
%% Time domain : FFR visualization

% Compute grand average and plot
grd_FFR = mean(all_subj,2);
figure ; 
plot(timepoints,grd_FFR,'Color',FFR_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR, 6-24 mo');

% Compute grand average by group
grpA.subj = subjects(contains(subjects,grpA.suffix));
grpB.subj = subjects(contains(subjects,grpB.suffix));

grpA.data = zeros(size(timepoints,1),1);
grpB.data = zeros(size(timepoints,1),1);

groups = {grpA.subj, grpB.subj};
data_groups = {grpA.data, grpB.data};
for k = 1:length(groups)
    grp_subj = zeros(size(timepoints,1),1);
    indices = find(ismember(subjects,groups{k})==1);
    for l = 1:length(indices)
        grp_subj(:,l) = all_subj(:,indices(l));
    end
    data_groups{k} = grp_subj;
end

FFR_avg_grpA = mean(data_groups{1,1},2);
FFR_avg_grpB = mean(data_groups{1,2},2);

% Visualization by group

% Temporal
fig = figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');

%Print plot into .svg, png and fig files
print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grp_comparison_', OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grp_comparison_', OPTIONS.polarity,'_FFR_temporal_', OPTIONS.params,'.png')));
saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grp_comparison_', OPTIONS.polarity,'_FFR_temporal_', OPTIONS.params,'.fig')));

% Visualization of each group
fig = figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 6-10mo');

print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grpA_', OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grpA_', OPTIONS.polarity,'_', OPTIONS.params,'.png')));
saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grpA_', OPTIONS.polarity,'_', OPTIONS.params,'.fig')));

fig = figure ; 
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR 18-24mo');

print('-dsvg',fullfile(OPTIONS.indir,strcat('mean_FFR_grpB_', OPTIONS.params,'.svg')));
print('-dpng',fullfile(OPTIONS.plot_dir,strcat('mean_FFR_grpB_', OPTIONS.polarity,'_', OPTIONS.params,'.png')));
saveas(fig, fullfile(strrep(OPTIONS.plot_dir,'png','fig'),strcat('mean_FFR_grpB_', OPTIONS.polarity,'_', OPTIONS.params,'.fig')));

%% Export mean FFRs into .txt files and convert into .avg files
fname_out_grpA = fullfile(OPTIONS.indir,strcat('mean_FFR_grpA_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_grpA,'w');
fprintf(fid,'%c\n',FFR_avg_grpA);
fclose(fid);

fname_out_grpB = fullfile(OPTIONS.indir,strcat('mean_FFR_grpB_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_grpB,'w');
fprintf(fid,'%c\n',FFR_avg_grpB);
fclose(fid);

fname_out_all = fullfile(OPTIONS.indir,strcat('mean_FFR_all_', OPTIONS.params,'.txt')) ;
fid = fopen(fname_out_all,'w');
fprintf(fid,'%c\n',grd_FFR);
fclose(fid);

BT_toolbox_path = fullfile(pwd, strcat('ToolBox_Brainstem\BT_2013\programFiles')) ;
addpath(BT_toolbox_path) ;
bt_txt2avg(fname_out_grpA, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_grpB, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);
bt_txt2avg(fname_out_all, OPTIONS.srate, OPTIONS.win_of_interest(1)*1000, OPTIONS.win_of_interest(2)*1000);

%% Root mean square and SNR on averaged FFRs

% Adapted from bt_rms.m from bt_gui toolbox (Kraus & Skoe, 2010)
% Note: std(X,1) is a shortcut RMS calc. It is equivalent to RMS
% on a baselined (demeaned to 0) waveform.  This is what we want here.

%note : the mean FFR includes the [-40;0] ms prestimulus period. In consequence:
%[0;40] = prestimulus period (->[-40;0]);
%[40;210] = whole FFR(->[0;170]);
%[50-95] = consonant-vowel transition period(->[10;55]);
%[95-210] = steady-state vowel period(->[55;170]);

startPrestim = 0.0610; % (in ms) to obtain 1 in samples
endPrestim = 40; % (in ms)

startCVtransition = 50;
endCVtransition = 95;

startVowel = 95;
endVowel = 210;

% Prestim period
Segment_A = FFR_avg_grpA(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000));
Segment_B = FFR_avg_grpB(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000));
RMSprestimA = std(Segment_A,1);
RMSprestimB = std(Segment_B,1);
%RMSprestimA = std(FFR_avg_grpA(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000)),1);
%RMSprestimB = std(FFR_avg_grpB(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000)),1);

% Whole response period
RMS_whole_grpA = std(FFR_avg_grpA(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
RMS_whole_grpB = std(FFR_avg_grpB(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
%RMS_grpA = std(FFR_avg_grpA,1); 
%RMS_grpB = std(FFR_avg_grpB,1); 

%CV transition period
RMS_CV_grpA = std(FFR_avg_grpA(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000)),1);
RMS_CV_grpB = std(FFR_avg_grpB(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000)),1);

%Vowel period
RMS_vowel_grpA = std(FFR_avg_grpA(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);
RMS_vowel_grpB = std(FFR_avg_grpB(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000)),1);

% Signal-to-noise on whole response
%SNR is obtained by dividing the RMS of a specific time window of the
%response by the RMS of the prestimulus response
% SNR_grpA = RMS_grpA./RMSprestimA;
% SNR_grpB = RMS_grpB./RMSprestimB;
%SNR_whole_grpA = RMS_whole_grpA./RMSprestimA;
%SNR_whole_grpB = RMS_whole_grpB./RMSprestimB;
SNR_whole_grpA = RMS_whole_grpA./RMSprestimA;
SNR_whole_grpB = RMS_whole_grpB./RMSprestimB;

% Signal-to-noise on CV transition period
SNR_CV_grpA = RMS_CV_grpA./RMSprestimA;
SNR_CV_grpB = RMS_CV_grpB./RMSprestimB;

% Signal-to-noise on vowel period
SNR_vowel_grpA = RMS_vowel_grpA./RMSprestimA;
SNR_vowel_grpB = RMS_vowel_grpB./RMSprestimB;

% display results
%disp(["Group A: RMS total = ",RMS_grpA, "RMS prestim = ", RMSprestimA, "SNR = ", SNR_grpA]);
%disp(["Group B: RMS total = ",RMS_grpB, "RMS prestim = ", RMSprestimB, "SNR = ", SNR_grpB]);

Group = ["6-10 mo";"18-24 mo"];
RMS_Prestim = [RMSprestimA;RMSprestimB];
RMS_Total = [RMS_whole_grpA;RMS_whole_grpB];
RMS_CV = [RMS_CV_grpA;RMS_CV_grpB];
RMS_Vowel = [RMS_vowel_grpA;RMS_vowel_grpB];
SNR_Total = [SNR_whole_grpA;SNR_whole_grpB];
SNR_CV = [SNR_CV_grpA;SNR_CV_grpB];
SNR_Vowel = [SNR_vowel_grpA;SNR_vowel_grpB];

RMS_and_SNR = table(Group,RMS_Prestim,RMS_Total,RMS_CV,RMS_Vowel,SNR_Total,SNR_CV,SNR_Vowel);
RMS_and_SNR

% Display distribution of RMS and SNR
%figure ; scatter(age_grA, neural_lags_grA) ; hold on ; scatter(age_grB, neural_lags_grB) ; legend({'6-10 mo', '18-24mo'}) ;
%figure ; boxplot(all_info.Var2, all_info.Var3, 'Notch','on','Labels',{'6-10 mo', '18-24mo'}) ;

% Save table in a .csv file
fname = fullfile(OPTIONS.indir, strcat('RMS_and_SNR_group_comparison_', OPTIONS.params,'.csv'));
writetable(RMS_and_SNR,fname, 'WriteVariableNames', true) ;

%% Root mean square and SNR on indivudual FFRs

Subj = [subjects(:)];

RMS_Prestim_all = zeros(1,length(subjects));
RMS_Total_all = zeros(1,length(subjects));
RMS_CV_all = zeros(1,length(subjects));
RMS_Vowel_all = zeros(1,length(subjects));
SNR_Total_all = zeros(1,length(subjects));
SNR_CV_all = zeros(1,length(subjects));
SNR_Vowel_all = zeros(1,length(subjects));

for s =1:length(subjects)
    RMS_Prestim_all(s) = std(all_subj(round(startPrestim*OPTIONS.srate/1000):round(endPrestim*OPTIONS.srate/1000),s),1);
    RMS_Total_all(s) = std(all_subj(round(endPrestim*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000),s),1);
    RMS_CV_all(s) = std(all_subj(round(startCVtransition*OPTIONS.srate/1000):round(endCVtransition*OPTIONS.srate/1000),s),1);
    RMS_Vowel_all(s) = std(all_subj(round(startVowel*OPTIONS.srate/1000):round(endVowel*OPTIONS.srate/1000),s),1);
    SNR_Total_all(s) = RMS_Total_all(s)/RMS_Prestim_all(s);
    SNR_CV_all(s) = RMS_CV_all(s)/RMS_Prestim_all(s);
    SNR_Vowel_all(s) = RMS_Vowel_all(s)/RMS_Prestim_all(s);
end

RMS_Prestim_all = RMS_Prestim_all';
RMS_Total_all = RMS_Total_all';
RMS_CV_all = RMS_CV_all';
RMS_Vowel_all = RMS_Vowel_all';
SNR_Total_all = SNR_Total_all';
SNR_CV_all = SNR_CV_all';
SNR_Vowel_all = SNR_Vowel_all';

RMS_and_SNR_all = table(Subj,RMS_Prestim_all,RMS_Total_all,RMS_CV_all,RMS_Vowel_all,SNR_Total_all,SNR_CV_all,SNR_Vowel_all);
RMS_and_SNR_all

% Save table in a .csv file
fname = fullfile(OPTIONS.indir, strcat('RMS_and_SNR_all_participants_', OPTIONS.params,'.csv'));
writetable(RMS_and_SNR_all,fname, 'WriteVariableNames', true) ;


% %todo : organiser data par groupe d'âge
% 

%% Neural lag

% Estimation of the transmission delay between stimulus and response. 
% Calculated from the time lag that produces the maximum stimulus-to-response cross-correlation magnitude

% Read files that contains neural lag and age information
neural_lags = readtable(fullfile(OPTIONS.indir, OPTIONS.nlag_filename)) ;
age_in_days = readtable(fullfile(OPTIONS.indir, 'age_in_days.csv'), 'Delimiter',';') ;

% Keep only subjects of interest (not rejected)
neural_lags = neural_lags(contains(neural_lags.suject_ID, subjects),:) ;
age_in_days = age_in_days(contains(age_in_days.subjects, subjects),:) ;

%Get variables of interest: neural lags and ages 
IDlist_grA = neural_lags.suject_ID(contains(neural_lags.group,'A')) ;
IDlist_grB = neural_lags.suject_ID(contains(neural_lags.group,'B')) ;
neural_lags_grA = neural_lags.neural_lag(contains(neural_lags.group,'A')) ;
neural_lags_grB = neural_lags.neural_lag(contains(neural_lags.group,'B')) ;
age_grA = age_in_days.age_in_days(contains(age_in_days.subjects,IDlist_grA)) ;
age_grB = age_in_days.age_in_days(contains(age_in_days.subjects,IDlist_grB)) ;

all_info = table(neural_lags.suject_ID, neural_lags.neural_lag, neural_lags.group, age_in_days.age_in_days) ;

 % Display neural lag distribution as a function of age
figure ; scatter(age_grA, neural_lags_grA) ; hold on ; scatter(age_grB, neural_lags_grB) ; legend({'6-10 mo', '18-24mo'}) ;
figure ; boxplot(all_info.Var2, all_info.Var3, 'Notch','on','Labels',{'6-10 mo', '18-24mo'}) ;

