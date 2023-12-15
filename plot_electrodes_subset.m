function [fig] = plot_electrodes_subset(signals, OPTIONS)
% ERPs visualization script 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Init number of figure(s)
nfig = 1;

% TODO : Check that all signals are the same length 

% Display in grid mode (following elect_subset format) 
for elec_letter=1:size(OPTIONS.elec_subset,1)

    for elec_numb=1:size(OPTIONS.elec_subset,2)

        hAxes = subplot(size(OPTIONS.elec_subset,1),size(OPTIONS.elec_subset,2),nfig) ;
        nfig = nfig +1 ;

        % Index of electrode to display
        idx_elec = OPTIONS.elec_indices(elec_letter,elec_numb) ; 

        for ss=1:length(signals)
            
            % if the matrix is 3D (plot variance)
            if numel(size(signals{ss}))>=3 
                % Plot timeseries
                plot(OPTIONS.vTime,mean(signals{ss}(idx_elec,:,:),3),'Color', OPTIONS.color{ss},'Linewidth',1.5) ; hold on ; set(gca,'YDir','reverse') ;

                % Plot transparent halo (+-SEM)
                plotHaloPatchSEM(hAxes, OPTIONS.vTime, squeeze(signals{ss}(idx_elec,:,:)), OPTIONS.color{ss}*255) ;
            else 
                % Plot timeseries
                plot(OPTIONS.vTime,signals{ss}(idx_elec,:),'Color', OPTIONS.color{ss},'Linewidth',1.5) ; hold on ; set(gca,'YDir','reverse') ;

            end

        end

        % Adjust scales (y-axis and x-axis) (transform in milliseconds)
        xlim(OPTIONS.xlim); ylim(OPTIONS.ylim) ; grid on ;

        % Add label of electrode in title 
        title(OPTIONS.elec_subset(elec_letter,elec_numb));

        % Display labels
        xlabel('Times (ms)'); ylabel('uV'); set(hAxes,'Fontsize',12);

    end

end
    % Legend : Add one single legend for 6 plots
    fig = gcf; fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend(OPTIONS.legend,'Location','bestoutside');
    Lgnd.Position(1) = 0.06; Lgnd.Position(2) = 0.8;

    % Title : Add a single title for 6 plots
    sgtitle(OPTIONS.title,'Interpreter', 'None', 'Fontsize', 16, 'FontWeight', 'bold');
end
