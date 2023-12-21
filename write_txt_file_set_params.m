function [] = write_txt_file_set_params(flag, count,suffix, OPTIONS)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2023 - %80PRIME Project

% If one element of the vector is 1 (=a subject to create)
if sum(flag) ~= 0
    
    fname_txt = fullfile(OPTIONS.indir,strcat('PARAM','_',OPTIONS.analysis,suffix,num2str(count),'.txt')) ;
   
    % If txt file does not exist -> creates it 
    if ~exist(fname_txt)
        t = table(repmat({'OPTIONS'}, [height(fieldnames(OPTIONS)), 1]), fieldnames(OPTIONS), struct2cell(OPTIONS));
        writetable(t, fname_txt, 'WriteVariableNames',0, 'Delimiter', ' ');
    end

end
