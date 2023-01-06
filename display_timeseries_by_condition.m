function [] = display_timeseries_by_condition(preproc_filenames, elec_subset, opt_balance, plots_dir)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

%Colors for plots
STD_color = [0.4941 0.1019 0.8863]; %purple
DEV1_color = [1 0.7686 0]; %light orange
DEV2_color = [1 0.4 0]; %dark orange
DEV_colors = {DEV1_color, DEV2_color};
DIFF_color = [0 0 0]; %black

YLIM_SCALE = [-20, 20] ;

cond_sylab = {'BA','GA'} ; 

for subject=1:length(preproc_filenames)

    % Gets files 
    for cond=1:length(cond_sylab)

        fname_DEV = strrep(preproc_filenames{subject},'reref_filtered_epoched',strcat('reref_filtered_epoched_DEV',num2str(cond),'_thresh_',opt_balance)) ; 
        fname_STD = strrep(preproc_filenames{subject},'reref_filtered_epoched',strcat('reref_filtered_epoched_STD',num2str(cond),'_thresh_',opt_balance)) ; 
    
        % Loads DEV trials
        [filepath,filename,ext] = fileparts(fname_DEV) ;
        EEG_DEV = pop_loadset(strcat(filename,ext),filepath) ;

        % Loads STD trials
        [filepath,filename,ext] = fileparts(fname_STD) ;
        EEG_STD = pop_loadset(strcat(filename,ext),filepath) ;
        
        nfig =1 ; 
        [~,subject_id,ext] = fileparts(fileparts(fname_DEV)) ; 

        figure('Name',strcat('Subject :',subject_id,'Condition :',strcat('DEV',num2str(cond))),'Units','normalized','Position',[0,0,1,1]);

        for elec_letter=1:size(elec_subset,1)
        
            for elec_numb=1:size(elec_subset,2)
                
                hAxes = subplot(size(elec_subset,1),size(elec_subset,2),nfig) ;
                nfig = nfig +1 ;
                
                idx_elec = find(ismember({EEG_DEV.chanlocs.labels},elec_subset(elec_letter,elec_numb))) ;
                
                % Compute grand average over one electrode
                grd_STD = squeeze(mean(EEG_STD.data(idx_elec,:,:),3)) ;
                grd_DEV = squeeze(mean(EEG_DEV.data(idx_elec,:,:),3)) ;
                grd_DIFF = grd_DEV - grd_STD ;
                
                % Plot timeseries
                plot(EEG_STD.times,grd_STD,'Color', STD_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DEV,'Color',DEV_colors{cond},'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DIFF,'Color',DIFF_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                
                % Plot transparent halo (+-mad)
                plotHaloPatchSEM(hAxes, EEG_STD.times, squeeze(EEG_STD.data(idx_elec,:,:)), STD_color*255) ;
                plotHaloPatchSEM(hAxes, EEG_DEV.times, squeeze(EEG_DEV.data(idx_elec,:,:)), DEV_colors{cond}*255);
                
                % Adjust graphics

                xlim([EEG_STD.xmin, EEG_STD.xmax]*1000); ylim(YLIM_SCALE) ; grid on ;
                title(elec_subset(elec_letter,elec_numb));
                xlabel('Times (ms)'); ylabel('uV');
                set(hAxes,'Fontsize',12);
                
            end
            
            %Add a single legend for 6 plots
            fig = gcf;
            fig.Position(3) = fig.Position(3) + 250;
            Lgnd = legend('STD (/DA/)',sprintf('DEV (/%s/)',cond_sylab{cond}),sprintf('DEV-STD (/%s/)',cond_sylab{cond}),'Location','bestoutside');
            Lgnd.Position(1) = 0.06;
            Lgnd.Position(2) = 0.8;
            
            %Add a single title for 6 plots
            sgtitle([strcat('Subject -> ',subject_id,' | Condition ->',strcat('DEV/STD',num2str(cond))),' (' ,opt_balance,' number of STDs)'],'Interpreter', 'None', 'Fontsize', 16, 'FontWeight', 'bold');
              
        end

        % Save data in vectoriel in subject folder
        out_fname = strrep(fname_DEV,strcat('reref_filtered_epoched_DEV',num2str(cond),'_thresh_balanced'),strcat('timeseries_condition',num2str(cond))); 
        print('-dsvg', strrep(out_fname,'.set','svg'));
        
        % Save data in png (with same filename as vectoriel) but different directory
        [~,fname,~] = fileparts(out_fname) ; 
        print('-dpng',fullfile(plots_dir,fname));

    end
end
