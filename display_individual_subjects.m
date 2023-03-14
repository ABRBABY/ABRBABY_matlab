function [] = display_individual_subjects(subjects_to_process, OPTIONS)
% ERPs visualization script 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%OPTIONS is a structure containing:
%params = 'RFE1_REJ1';                            % option of preprocess to consider
%elec_subset = {'F3','Fz','F4';'C3','Cz','C4'};   % electrodes to display
%indir = indir ;                                  % directory path of files to process
%diff_display = 1 ;                               % 1 to display difference wave (MMN), 0 to not display
%plot_dir = plot_dir ;                            % path to save png files of plots
%balance_STD = 'unbalanced';                      % 'balanced' or 'unbalanced' number of STD
%ylim = [-20,20] ;                                % limits of y axis


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
        
        fname_DEV = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_DEV',num2str(cc),'*_',OPTIONS.balance_STD,'_',OPTIONS.params,'.set'))) ; 
        fname_STD = dir(fullfile(OPTIONS.indir,subjects_to_process{ss},strcat(subjects_to_process{ss},'_STD',num2str(cc),'*_',OPTIONS.balance_STD,'_',OPTIONS.params,'.set'))) ; 
        
        % Loads DEV trials
        EEG_DEV = pop_loadset(fname_DEV.name,fullfile(OPTIONS.indir,subjects_to_process{ss})) ;

        % Loads STD trials
        EEG_STD = pop_loadset(fname_STD.name,fullfile(OPTIONS.indir,subjects_to_process{ss})) ;
        
        nfig =1 ; 
        
        % Create figure for one condition (e.g. DEV1)
        figure('Name',strcat('Subject :',subjects_to_process{ss},' | Condition :',strcat('DEV',num2str(cc))),'Units','normalized','Position',[0,0,1,1]);

        % Compute grand average over one electrode
        grd_STD = mean(EEG_STD.data,3) ;
        grd_DEV = mean(EEG_DEV.data,3) ;
        grd_DIFF = grd_DEV - grd_STD ;
                
        for elec_letter=1:size(OPTIONS.elec_subset,1)
        
            for elec_numb=1:size(OPTIONS.elec_subset,2)
                
                hAxes = subplot(size(OPTIONS.elec_subset,1),size(OPTIONS.elec_subset,2),nfig) ;
                nfig = nfig +1 ;
                
                % Find index of electrode to display
                idx_elec = find(ismember({EEG_DEV.chanlocs.labels},OPTIONS.elec_subset(elec_letter,elec_numb))) ;
               
                % Plot timeseries
                plot(EEG_STD.times,grd_STD(idx_elec,:),'Color', STD_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DEV(idx_elec,:),'Color',DEV_colors{cc},'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DIFF(idx_elec,:),'Color',DIFF_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                
                % Plot transparent halo (+-mad)
                plotHaloPatchSEM(hAxes, EEG_STD.times, squeeze(EEG_STD.data(idx_elec,:,:)), STD_color*255) ;
                plotHaloPatchSEM(hAxes, EEG_DEV.times, squeeze(EEG_DEV.data(idx_elec,:,:)), DEV_colors{cc}*255);
                
                % Adjust scales (y-axis and x-axis) (transform in milliseconds)
                xlim([EEG_STD.xmin, EEG_STD.xmax]*1000); ylim(OPTIONS.ylim) ; grid on ;
                
                % Add label of electrode in title 
                title(OPTIONS.elec_subset(elec_letter,elec_numb));
                
                % Display labels
                xlabel('Times (ms)'); ylabel('uV'); set(hAxes,'Fontsize',12);
                
            end
            
            % Legend : Add one single legend for 6 plots
            fig = gcf; fig.Position(3) = fig.Position(3) + 250;
            Lgnd = legend('STD (/DA/)',sprintf('DEV (/%s/)',cond_sylab{cc}),sprintf('DEV-STD (/%s/)',cond_sylab{cc}),'Location','bestoutside');
            Lgnd.Position(1) = 0.06; Lgnd.Position(2) = 0.8;
            
            % Title : Add a single title for 6 plots
            sgtitle([strcat('Subject -> ',subjects_to_process{ss},' | Condition ->',strcat('DEV/STD',num2str(cc))),' (' ,OPTIONS.balance_STD,' number of STDs)'],'Interpreter', 'None', 'Fontsize', 16, 'FontWeight', 'bold');
              
        end

        % Save data in vectoriel in subject folder
        out_fname = fullfile(OPTIONS.indir,subjects_to_process{ss},strrep(fname_DEV.name,'.set','.svg'));
        print('-dsvg', out_fname);
        
        % Save data in png (with same filename as vectoriel) but different directory
        print('-dpng',fullfile(OPTIONS.plot_dir,strrep(fname_DEV.name,'.set','.png')));

    end
end
