function visualisation = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,mean_grpA, mean_grpB,timepoints,fig_name,legend1,legend2,STD_number, grpA, grpB)

visualisation = figure('Name', fig_name,'Units','normalized','Position',[0,0,1,1]);

span = 1000;
%filter_disp = 10; %(Hz)

grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.2 0.7765 0.2]; %green

nfig = 1 ; 

for elec_letter=1:size(elec_to_disp_labels,1)
    
    for elec_numb=1:size(elec_to_disp_labels,2)
        hAxes = subplot(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2),nfig) ; 
        nfig = nfig +1 ;
        idx_elec = elec_indices(elec_letter,elec_numb);
        
        % Compute grand average over one electrode
        grd_avg_grpA = squeeze(mean(mean_grpA(idx_elec,:),1)) ; 
        grd_avg_grpB = squeeze(mean(mean_grpB(idx_elec,:),1)) ; 
                
        %Add 10 Hz lowpass filter for display
        %functions to try :
        %smooth(timepoints, grd_avg_grpA,1000); %1000 points moving average (ex from pop_smootherp in ERPLAB)
        %lowpass(grd_avg_grpA,10,16384);
        grd_avg_grpA_smooth = smoothdata(grd_avg_grpA, 'movmean', span);
        grd_avg_grpB_smooth = smoothdata(grd_avg_grpB, 'movmean', span);
            
        % Plot timeseries smoothed
        plot(timepoints,grd_avg_grpA_smooth,'Color', grpA_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
        plot(timepoints,grd_avg_grpB_smooth,'Color', grpB_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
       
        % Plot timeseries (dotted line)
        plot(timepoints,grd_avg_grpA,':','Color', grpA_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
        plot(timepoints,grd_avg_grpB,':','Color', grpB_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
    
        % Plot transparetn halo (+-mad)
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpA(:,idx_elec,:)+raw_avg2_grpA(:,idx_elec,:))/2,1)), [0,255,0]) ; 
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpB(:,idx_elec,:)+raw_avg2_grpB(:,idx_elec,:))/2,1)), [255,0,0]) ; 
        %plotHaloPatchMAD(hAxes, timepoints, grd_avg_grpA, grpA_color*255)
        %plotHaloPatchMAD(hAxes, timepoints, grd_avg_grpA, grpA_color*255)
        std_dev_patch(grd_avg_grpA, grpA_color);
        display_std_error(grd_avg_grpB, grpB_color);
        
        % Adjust graphics
         %xlim([grd_STD_grpA.xmin, grd_STD_grpB.xmax]*1000); grid on ; 
         grid on;
         ylim([-4 5]);
         title(elec_to_disp_labels(elec_letter,elec_numb));
         xlabel('Times (ms)'); ylabel('uV');
         set(hAxes,'Fontsize',12);
                 
    end
    
    %Add a single legend for 6 plots
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend({legend1,legend2},'Location','bestoutside');
    Lgnd.Position(1) = 0.06;
    Lgnd.Position(2) = 0.8;
    
    %Add a single title for 6 plots
    if STD_number == 1
        sgtitle('Balanced number of STDs', 'Fontsize', 16, 'FontWeight', 'bold');
    elseif STD_number == 2
        sgtitle('Unbalanced number of STDs (all)', 'Fontsize', 16, 'FontWeight', 'bold')
    end
    
    
end


%%
% %function [] = display_std_error(y);
% y = rand(1,10); % your mean vector;
% x = 1:numel(y);
% %stderror= std(data) / sqrt( length( data ))
% %std_dev = 1;
% std_dev =std(y)/sqrt(length(y));
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% %patch(curve1, fliplr(curve2), 'g');
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);
% %end

end