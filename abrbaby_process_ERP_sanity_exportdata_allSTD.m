function [] = abrbaby_process_ERP_sanity_exportdata_allSTD(eeglab_path, biosig_installer_path,indir,plots_dir,chan_dir,STD_number) 
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%----------------------------------------------------------------------------------------------------
% %WHEN NOT IN 'FUNCTION' MODE:
% eeglab_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1';
% biosig_installer_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1\plugins\Biosig3.7.9\biosig_installer.m';
% chan_dir ='\\Filer\\home\\Invites\\hervé\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc';
% plots_dir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\png_plots_eeg_data';
% indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';
% %indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';
% STD_number = 1;

%----------------------------------------------------------------------------------------------------

% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
tmp = pwd ; 
cd(eeglab_path) ; 
% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
run(biosig_installer_path) ; 
cd(tmp) ;

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%% Variables to enter manually before running the code

% Set variables for filtering
% hp = 0.1; %value for high-pass filter (Hz) (APICE)
% lp = 40; %value for low-pass filter (Hz) (APICE) 
hp = 1; %value for high-pass filter (Hz) (APICE)
lp = 30; %value for low-pass filter (Hz) (APICE) 

elec_to_disp_labels = {'F3','Fz','F4';'C3','Cz','C4'};

% Rejection treshold for bad epochs
rej_low = -150; %150 infants; 120 adults
rej_high = 150; %150 infants; 120 adults

% List of channel labels to reref with 
mastos = {'Lmon','Rmon','MASTOG','MASTOD'};
trig = {'Erg1'};

baseline = [-99, 0] ; 
win_of_interest = [-0.1, 0.5] ; 
conditions = {'STD','DEV1','DEV2'} ; 
cond_sylab = {'BA','GA'} ; 

elec = 1:16 ; 

%Colors for plots
STD_color = [0.4941 0.1019 0.8863]; %purple
DEV1_color = [1 0.7686 0]; %light orange
DEV2_color = [1 0.4 0]; %dark orange
DEV_colors = {DEV1_color, DEV2_color};
DIFF_color = [0 0 0]; %black

% FOR SANITY CHECK
%for jj=find(ismember(subjects,'DVL_010_T24'))

    %Loop through subjects
    for jj=1:length(subjects) 
        
     fprintf(strcat(subjects{jj}, '...\n'));
     
    %% IMPORT
    % Get BDF file
    %[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
    fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
 
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

    % Find REF electrodes indices by labels 
    ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 
    
    % Save a first dataset in EEGLAB 
    [~,filename,~] = fileparts(fname.name);    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');
 
    %% RE-REF (excluding trig channel)
    % Find TRIG electrodes indices by labels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 

    % Re-reference data and rename new file
    EEG = pop_reref(EEG, ref_elec, 'exclude',trigg_elec, 'keepref','on');
    
    %% EVENTS 
    % Extract event from trigger channel (Erg1)
    EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);
    
    % Identifies outliers events (e.g. boundaries) or too close events 
    %idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate)+1,... % minimum intretrial duration = 220 ms
     %                   find(diff([EEG.event.latency])>2*EEG.srate) ];        
    idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate),... % minimum intretrial duration = 220 ms
                        find(diff([EEG.event.latency])>2*EEG.srate) ];  
    % Removes outliers events
    EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ; 
    
    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;
    
    %% FILTERS the data with ERPLab
    EEG  = pop_basicfilter(EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    %% SAVE DATASET BEFORE EPOCHING
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_filtered'),'savenew', fullfile(indir,subjects{jj}, strcat(filename,'_filtered')),'gui','off');
    CURR_FILTERED = CURRENTSET ; 
    
    % Extract ALL conditions epochs
    EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');

    % Remove baseline
    EEG = pop_rmbase( EEG, baseline,[] );
    
    % Add channels information
    %EEG=pop_chanedit(EEG, 'lookup','/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'lookup',chan_dir);

    % Select DEV
    [EEG_DEV1,target_indices1] = pop_selectevent(EEG,'type','DEV1');
    [EEG_DEV2,target_indices2] = pop_selectevent(EEG,'type','DEV2');
    [EEG_STD,target_indices_std] = pop_selectevent(EEG,'type','STD');
    
    if STD_number == 1
        idx_std1 = target_indices_std(ismember(target_indices_std,target_indices1-1));
        idx_std2 = target_indices_std(ismember(target_indices_std,target_indices2-1));
    elseif STD_number == 2
        begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
        idx_std1 = setdiff(target_indices_std,begining_of_block);
        idx_std2 = setdiff(target_indices_std,begining_of_block);
    end
    
    [EEG_STD1,target_indices_std1] = pop_selectevent(EEG,'event',idx_std1);
    [EEG_STD2,target_indices_std2] = pop_selectevent(EEG,'event',idx_std2);
   
    [EEG_STD1_thresh,idx_std1_rej] = pop_eegthresh(EEG_STD1,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_STD2_thresh,idx_std2_rej] = pop_eegthresh(EEG_STD2,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV1_thresh,idx_dev1_rej] = pop_eegthresh(EEG_DEV1,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV2_thresh,idx_dev2_rej] = pop_eegthresh(EEG_DEV2,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
   
    if STD_number == 1
        [EEG_STD_thresh, idx_removed] = pop_eegthresh(EEG_STD,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
        std_good = setdiff(1:900,target_indices_std(idx_removed)); 
        begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
    end
    
     %Export information about rejected trials
      % Get indices of the trials which were rejected (without messing around with the relative indices)
    [EEG_all,idx_rejected_all] = pop_eegthresh(EEG,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
      % Extract variables of interest
    trial_index = 1:EEG.trials;
    trial_num = [EEG.event.urevent];
    condition = {EEG.event.type} ;
    latency = [EEG.event.latency]/EEG.srate; 
    %latency = [EEG.event.latency]; 
    rejected = ismember(trial_index,idx_rejected_all) ;
    bloc = zeros(900,1);
    tr = 1;
    for bc = 1:30
        bloc(tr:tr+29,1)= repmat(bc,30,1);
        bc = bc+1;
        tr = tr+30;
    end
    trial_index = trial_index';
    condition = condition';
    latency = latency';
    trial_num = trial_num';
    rejected = rejected';
      % Create table to store these information
    list_trial_infos = table(trial_index,condition,latency, trial_num,rejected, bloc) ;
      %  Save this table into a csv file (use function writetable)
    writetable(list_trial_infos,fullfile(indir,subjects{jj},strcat(filename,'_low_',num2str(rej_low),'_high_',num2str(rej_high),'infos_trials.csv'))) ; 

    if STD_number == 1
        
        % If nubmber of STD1 < number of DEV1 : randomly select other STD
        if length(EEG_DEV1_thresh.event)>length(EEG_STD1_thresh.event)
            % Find number of trial to add in STD
            ntrials =  length(EEG_DEV1_thresh.event)-length(EEG_STD1_thresh.event) ;
            %Apply function to balance number of STD trial to reach the number of DEV trials
            EEG_STD1_thresh = balance_number_of_STD(EEG,ntrials,std_good,target_indices1,begining_of_block,target_indices_std1,idx_std1_rej,idx_std1) ;
            
        end
        
        % If nubmber of STD2 < number of DEV2 : randomly select other STD
        if length(EEG_DEV2_thresh.event)>length(EEG_STD2_thresh.event)
            % Find number of trial to add in STD
            ntrials =  length(EEG_DEV2_thresh.event)-length(EEG_STD2_thresh.event) ;
            %Apply function to balance number of STD trial to reach the number of DEV trials
            EEG_STD2_thresh = balance_number_of_STD(EEG,ntrials,std_good,target_indices2,begining_of_block,target_indices_std2,idx_std2_rej,idx_std2) ;
            
        end
        
        % If nubmber of STD1 > number of DEV1 : delete random STD preceding missing DEV1
        if length(EEG_DEV1_thresh.event)<length(EEG_STD1_thresh.event)
            %Find indices of rejected DEV1 in the 900 trials referential
            list_of_rej_dev1 = target_indices1(idx_dev1_rej);
            %Select random indices amont the rejected DEV1
            random_std1_to_remove = randsample(list_of_rej_dev1-1,size(idx_dev1_rej,2)-size(idx_std1_rej,2));
            %Remove STD corresponding to rejected DEV1
            idx_std_already_included = setdiff(target_indices_std1, target_indices_std1(idx_std1_rej)) ;
            new_list_of_std1 = setdiff(idx_std_already_included,random_std1_to_remove);
            %Modify EEG struct
            [EEG_STD1_thresh,~] = pop_selectevent(EEG,'event',[new_list_of_std1]);
        end
        
        % If nubmber of STD2 > number of DEV2 : delete random STD preceding
        % missing DEV2
        if length(EEG_DEV2_thresh.event)<length(EEG_STD2_thresh.event)
            %Find indices of rejected DEV1 in the 900 trials referential
            list_of_rej_dev2 = target_indices2(idx_dev2_rej);
            %Select random indices amont the rejected DEV1
            random_std2_to_remove = randsample(list_of_rej_dev2-1,size(idx_dev2_rej,2)-size(idx_std2_rej,2));
            %Remove STD corresponding to rejected DEV1
            idx_std_already_included = setdiff(target_indices_std2, target_indices_std2(idx_std2_rej)) ;
            new_list_of_std2 = setdiff(idx_std_already_included,random_std2_to_remove);
            %Modify EEG struct
            [EEG_STD2_thresh,~] = pop_selectevent(EEG,'event',[new_list_of_std2]);
        end
        
        
    end
   
    for cc=1:2
        
        if STD_number == 1
            figure('Name',strcat('Subject :',subjects{jj},'Condition :',conditions{cc+1}),'Units','normalized','Position',[0,0,1,1]);
        elseif STD_number == 2
            figure('Name',strcat('Subject :',subjects{jj},'Condition :',conditions{cc+1},'_allSTD'),'Units','normalized','Position',[0,0,1,1]);
        end
        
        EEG_DEV = eval(sprintf('EEG_DEV%d_thresh',cc)) ;
        EEG_STD = eval(sprintf('EEG_STD%d_thresh',cc)) ;
        nfig = 1 ;
        
        for elec_letter=1:size(elec_to_disp_labels,1)
            
            for elec_numb=1:size(elec_to_disp_labels,2)
                
                hAxes = subplot(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2),nfig) ;
                nfig = nfig +1 ;
                
                idx_elec = find(ismember({EEG_DEV.chanlocs.labels},elec_to_disp_labels(elec_letter,elec_numb))) ;
                
                % Compute grand average over one electrode
                grd_STD = squeeze(mean(EEG_STD.data(idx_elec,:,:),3)) ;
                grd_DEV = squeeze(mean(EEG_DEV.data(idx_elec,:,:),3)) ;
                grd_DIFF = grd_DEV - grd_STD ;
                
                % Plot timeseries
                plot(EEG_STD.times,grd_STD,'Color', STD_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DEV,'Color',DEV_colors{cc},'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DIFF,'Color',DIFF_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                
                % Plot transparent halo (+-mad)
                %plotHaloPatchMAD(hAxes, EEG_STD.times, squeeze(EEG_STD.data(idx_elec,:,:)), STD_color*255) ;
                %plotHaloPatchMAD(hAxes, EEG_DEV.times, squeeze(EEG_DEV.data(idx_elec,:,:)), DEV_colors{cc}*255);
                
                % Adjust graphics
                xlim([EEG_STD.xmin, EEG_STD.xmax]*1000); grid on ;
                title(elec_to_disp_labels(elec_letter,elec_numb));
                xlabel('Times (ms)'); ylabel('uV');
                set(hAxes,'Fontsize',12);
                
            end
            
            %Add a single legend for 6 plots
            fig = gcf;
            fig.Position(3) = fig.Position(3) + 250;
            Lgnd = legend('STD (/DA/)',sprintf('DEV (/%s/)',cond_sylab{cc}),sprintf('DEV-STD (/%s/)',cond_sylab{cc}),'Location','bestoutside');
            Lgnd.Position(1) = 0.06;
            Lgnd.Position(2) = 0.8;
            
            %Add a single title for 6 plots
            if STD_number == 1
                sgtitle([filename,' balanced number of STDs'],'Interpreter', 'None', 'Fontsize', 16, 'FontWeight', 'bold');
            elseif STD_number == 2
                sgtitle([filename,' unbalanced number of STDs (all)'],'Interpreter', 'None', 'Fontsize', 16, 'FontWeight', 'bold')
            end
            
            % To save data in vectoriel and png
            if STD_number == 1
                print('-dsvg',fullfile(indir,subjects{jj},strcat(subjects{jj},'_',conditions{cc+1},'_balancedSTD.svg')));
                %print('-dpng',fullfile(indir,subjects{jj},strcat(subjects{jj},'_',conditions{cc+1},'_balancedSTD.png')));
                print('-dpng',fullfile(plots_dir,strcat(subjects{jj},'_',conditions{cc+1},'_balancedSTD.png')));
            elseif STD_number == 2
                print('-dsvg',fullfile(indir,subjects{jj},strcat(subjects{jj},'_',conditions{cc+1},'_allSTD.svg')));
                %print('-dpng',fullfile(indir,subjects{jj},strcat(subjects{jj},'_',conditions{cc+1},'_allSTD.png')));
                print('-dpng',fullfile(plots_dir,strcat(subjects{jj},'_',conditions{cc+1},'_allSTD.png')));
            end
            
        end
        
        %Export data
        [ALLEEG, EEG_DEV, CURRENTSET] = pop_newset(ALLEEG, EEG_DEV, CURRENTSET, 'setname',strcat(filename,'_','DEV',num2str(cc)),'savenew', fullfile(indir,subjects{jj},strcat(filename,'_','DEV',num2str(cc))),'gui','off');
        if STD_number == 1
            [ALLEEG, EEG_STD, CURRENTSET] = pop_newset(ALLEEG, EEG_STD, CURRENTSET, 'setname',strcat(filename,'_','STD',num2str(cc)),'savenew', fullfile(indir,subjects{jj},strcat(filename,'_','STD',num2str(cc))),'gui','off');
        elseif STD_number == 2
            [ALLEEG, EEG_STD, CURRENTSET] = pop_newset(ALLEEG, EEG_STD, CURRENTSET, 'setname',strcat(filename,'_','all_STD',num2str(cc)),'savenew', fullfile(indir,subjects{jj},strcat(filename,'_','all_STD',num2str(cc))),'gui','off');
        end
        
    end

    end   
end

    
%--------------------------------------------------------------
% FUNCTION that select from EEG_STD ntrial with no repetition with exisitng
% STD 
%--------------------------------------------------------------
function [EEG_STD_ALL] = balance_number_of_STD(EEG,ntrials,std_good,target_indices,begining_of_block,target_indices_std,idx_std_rej, idx_std)
   
        % Pool of STD without those rejected by threshold detection
        pool_std = setdiff(std_good,target_indices);   
        
        % Pool of STD without beginners in block (3 first trials) 
        pool_std_w_no_beginners = setdiff(pool_std,begining_of_block);
        
        % Trials which were already selected 
        idx_std_already_included = setdiff(target_indices_std, target_indices_std(idx_std_rej)) ; 
     
        % Trials to add to balance the number of trial to the same number
        % as DEV
        idx_to_add = pool_std_w_no_beginners(randperm(length(pool_std_w_no_beginners),ntrials));
        
        % Select trial : 1) randomly a number = ntrial and 2) those which
        % were already selected and 'good'
        [EEG_STD_ALL,~] = pop_selectevent(EEG,'event',[idx_std_already_included idx_to_add]);
end
    
%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = read_custom_events(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 0);

% Insert info from .txt into EEG.event
my_events = table2array(my_events);

out_event = struct('latency', {in_event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {in_event(:).urevent});

end