%%draft_export_stats_bst

%Example: create figure of temporal course + t values for condition 2 (group
%comparison) at F3 electrode
% 1. Export data to matlab (from Brainstorm interface) 
% for 2 conditions : tt1 and tt2
% for 2 groups : g1 and g2
% 2. Plot data
f3 = strcmp({chan.Channel.Name},'F3') ;
tthresh = (tt1.tmap(f3,:).*tt1.pmap(f3,:)<0.05)*-2 ; figure ; plot(tthresh,'g'); hold on ; plot(g1.F(f3,:)*10^6,'k') ; hold on ;  plot(g2.F(f3,:)*10^6,'r') ;  xticks(linspace(0,180,length(xlab))) ; xticklabels(xlab) ; grid on ; ylim([-2;2]) ; set(gca,'ydir','reverse') ; xlabel('Time (ms)') ; ylabel('uV') ; title('F3 - Condition 2');