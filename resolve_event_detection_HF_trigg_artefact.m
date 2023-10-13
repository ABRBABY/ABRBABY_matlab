function [idx_to_remove] = resolve_event_detection_HF_trigg_artefact(EEG)

% When issue occurs during event detection because of noisy HF triggers,
% run this code to change idx_to_remove variable and detect the right
% outliers to delete

n = 1;
m = 1;
ii = 1 ;
idx_to_remove = [];

while m<size(EEG.event,2)
    m = n+1 ;
    while EEG.event(m).latency-EEG.event(n).latency < 3604
        idx_to_remove(ii) = m ;
        ii = ii+1 ;
        m = m+1;
    end
    n=m ;
end


end