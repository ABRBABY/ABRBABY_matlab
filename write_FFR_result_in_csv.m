function [] = write_FFR_result_in_csv(OPTIONS, flag, value, fname)

% Check if values length correspond to number of subjects (given by flag)
if length(value)~=sum(flag)
   error("Something went wrong, there is not the same number of values than subjects which were computed (as indicated per flag)");
end

% Get list of all subjects present in the database
all_sub = get_subjects(OPTIONS.indir,[]);

% Create a table with values of selected (by flag) subjects
myVal = table(all_sub(flag), value', 'VariableNames', {'suject_ID', 'value'}) ;

% Write a table containing values
writetable(myVal,fullfile(OPTIONS.indir,fname), 'WriteVariableNames', true) ;
fprintf(sprintf('%s  created...\n',fullfile(OPTIONS.indir,fname)));
  