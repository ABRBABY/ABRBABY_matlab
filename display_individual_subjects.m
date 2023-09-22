function [] = display_individual_subjects(subjects_to_process, OPTIONS)
% ERPs visualization script 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% OPTIONS is a structure containing:
% params = 'RFE1_REJ1';                            % option of preprocess to consider
% elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
% indir = indir ;                                  % directory path of files to process
% diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
% plot_dir = plot_dir ;                            % path to save png files of plots
% balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
% ylim = [-20,20] ;                                % limits of y axis


%Colors for plots
STD_color = [0.4941 0.1019 0.8863]; %purple
DEV1_color = [1 0.7686 0]; %light orange
DEV2_color = [1 0.4 0]; %dark orange
DEV_colors = {DEV1_color, DEV2_color};
DIFF_color = [0 0 0]; %black

cond_sylab = {'BA','GA'} ; 

for ss=1:length(subjects_to_process)

    % Gets files 
    for cc=1:length(cond_sylab)
        
        fname_DEV = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_DEV',num2str(cc),'_thresh_',OPTIONS.balance_STD,'_',OPTIONS.params,'*.set'))) ; 
        if contains(OPTIONS.balance_STD,'unbalanced') == 1
            fname_STD = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_STDD_thresh_',OPTIONS.balance_STD,'_',OPTIONS.params,'*.set'))) ; 
            if isempty(fname_STD) == 1
                fname_STD = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_STD1_thresh_',OPTIONS.balance_STD,'_',OPTIONS.params,'*.set'))) ; 
            end
        else
            fname_STD = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_STD',num2str(cc),'_thresh_',OPTIONS.balance_STD,'_',OPTIONS.params,'*.set'))) ; 
        end 

        if size(fname_STD,1) > 1
        rman = find(contains({fname_STD.name}, 'rman')) ;
        fname_STD = fname_STD(rman,1) ;
        end

        if size(fname_DEV,1) > 1
        rman = find(contains({fname_DEV.name}, 'rman')) ;
        fname_DEV = fname_DEV(rman,1) ;
        end
        
        % Loads DEV trials
        EEG_DEV = pop_loadset(fname_DEV.name,fullfile(OPTIONS.indir,subjects_to_process{ss})) ;

        % Loads STD trials
        EEG_STD = pop_loadset(fname_STD.name,fullfile(OPTIONS.indir,subjects_to_process{ss})) ;
        
        nfig =1 ; 
        
        % Create figure for one condition (e.g. DEV1)
        figure('Name',strcat('Subject :',subjects_to_process{ss},' | Condition :',strcat('DEV',num2str(cc))),'Units','normalized','Position',[0,0,1,1]);

        
        % Compute grand average over one electrode
        signals =       {EEG_STD.data,  EEG_DEV.data,   mean(EEG_DEV.data,3)- mean(EEG_STD.data,3)};
        OPTIONS.color = {STD_color,     DEV_colors{cc}, DIFF_color} ; 
        
        % Get indices of electrodes in subset
        [sharedvals,idxA,idxB] = intersect({EEG_STD.chanlocs.labels}, OPTIONS.elec_subset(:),'stable') ; 
        [~,id1,id2] = intersect(OPTIONS.elec_subset(:), sharedvals,'stable') ;
        OPTIONS.elec_indices = reshape(idxA(id2),size(OPTIONS.elec_subset)) ;
        
        OPTIONS.vTime = EEG_STD.times; 
        OPTIONS.xlim = [EEG_STD.xmin, EEG_STD.xmax]*1000 ;  
        OPTIONS.legend{1} = 'STD (/DA/)' ; 
        OPTIONS.legend{2} = sprintf('DEV (/%s/)',cond_sylab{cc}) ;
        OPTIONS.legend{3} = sprintf('DEV-STD (/%s/)',cond_sylab{cc}) ; 
        OPTIONS.title = [strcat('Subject -> ',subjects_to_process{ss},' | Condition ->',strcat('DEV/STD',num2str(cc))),' (' ,OPTIONS.balance_STD,' number of STDs)'] ; 
        
        % Call visualisation function (grids with electrode subset) 
        [fig] = plot_electrodes_subset(signals,OPTIONS) ; 
        
        % Save figures 
        if OPTIONS.savefigs == 1 
            % Save data in vectoriel in subject folder
            out_fname = fullfile(OPTIONS.indir,subjects_to_process{ss},strrep(fname_DEV.name,'.set','.svg'));
            print('-dsvg', out_fname);
            
            % Save data in png (with same filename as vectoriel) but different directory
            print('-dpng',fullfile(OPTIONS.plot_dir,strrep(fname_DEV.name,'.set','.png')));
    
            % Save data in fig (with same filename as vectoriel) but different directory
            saveas(fig, fullfile(strrep(OPTIONS.plot_dir, 'png', 'fig'),strrep(fname_DEV.name,'.set','.fig')));
        end

    end
end
