
function [] = compute_neural_lag_report(subjects_to_process, OPTIONS)

path_stim = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY_matlab/ToolBox_BrainStem/BT_2013/da_170_kraus_16384_LP3000_HP80.avg' ; 
start = 0;
stop = 40;
lagstart = 0;
lagstop = 9.6000;
polarity = 'POSITIVE' ; 
chan =1;
chancomp =1;

for ss = 1:length(subjects_to_process) %for each subject
    
    
    FFR_file = fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_',OPTIONS.params,'_abr_shifted_data_HF.avg')) ;
    
    [LAG_atmaxmin, maxmincor, all_corrs, all_lags] = bt_xcorrelation2(FFR_file, path_stim, start, stop, lagstart, lagstop, polarity, chan, chancomp) ; 
    
    subject_lag(ss) =  LAG_atmaxmin ; 
    
    
end


% HERE Write a table with all lags 
