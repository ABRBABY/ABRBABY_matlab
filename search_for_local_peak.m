function [lat,amp, auc] = search_for_local_peak(signal, vTimes, timeW, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

if size(signal,1)~= 1; error('search_for_local_peak : wrong size of signal') ; end

[pks,locs] = findpeaks(-signal(:,vTimes>timeW(1)&vTimes<timeW(2)));

if isempty(pks) % No local min 
    [amp, slat] = min(signal(:,vTimes>timeW(1)&vTimes<timeW(2)),[],2);
elseif size(pks,2)>1 % several local min
    [~,idx]=min(-pks) ; amp  = -pks(idx) ; slat = locs(idx) ; 
else amp=-pks ; slat = locs;

end

if OPTIONS.disp ; figure ; plot(vTimes(:,vTimes>timeW(1)&vTimes<timeW(2)), signal(:,vTimes>timeW(1)&vTimes<timeW(2))) ; hold on ; plot(vTimes(locs-1)+find(vTimes>timeW(1),1),-pks,'r*') ; end

% Get latency in the right time referential 
lat =vTimes(slat)+find(vTimes>timeW(1),1) ; 

% Compute auc
auc = mean(signal(:,vTimes>lat-OPTIONS.auc_delta&vTimes<lat+OPTIONS.auc_delta)) ; 
