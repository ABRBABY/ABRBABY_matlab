interpol.subj = {'DVL_004_T8','DVL_006_T10','DVL_006_T18','DVL_007_T18','DVL_008_T10','DVL_008_T18','DVL_010_T18','DVL_012_T10','DVL_013_T10','DVL_013_T8','DVL_018_T10','DVL_037_T8','DVL_021_T18','DVL_045_T8', 'DVL_031_T24'} ;
interpol.chan = {{'T7';'T8'},  {'T7';'T8'},    {'T8'},       {'T7'},      {'Cz'},       {'T8'},        {'T7'},      {'T8'},     {'T7';'T8'},   {'T7';'T8'},  {'T7'},      {'T8'},      {'F3';'O2'},    {'Cz'},     {'T7','T8'}   } ;


for jj=1:length(subj) 
 
     fid = fopen(strcat(subj{jj},'_interp.txt'),'wt');
     fprintf(fid, cell2mat(strcat(chan{jj},',')));
     fclose(fid);
 end
 