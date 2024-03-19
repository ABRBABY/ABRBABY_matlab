function [] = add_flag_column_trials_description(fname, header,flag, overwrite)

    % Read _trials_description file 
    T1 = readtable(fname); 
    
    % Find if column exists
    exist_column = strcmp(header,T1.Properties.VariableNames) ; 

     % % If the variable does not exost in the table create a new 
    if sum(exist_column)==0
        T2 = table(flag','VariableNames',{header});
        T1= [T1 T2]; 
        
    elseif overwrite % If column exist AND is different from the existing display a WARNINIG message
        if any(T1{:,exist_column}~=flag')
            fprintf('WARNING : trial_decription was replaced by new values');
        end
        % Add or replace column values
        T1{:,exist_column}=flag';
    else
        T1= T1{:,exist_column}+flag'; 
    end
    writetable(T1,fname,'WriteVariableNames', true);