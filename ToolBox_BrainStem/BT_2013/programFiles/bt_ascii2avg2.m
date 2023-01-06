function bt_ascii2avg(ASCIIfile)
% this function takes an ascii file and converts its contents to 
%  1) a Neuroscan-format .avg file, and
%  2) a text file containing the marked latencies
% Both resulting files will retain the same base filename as the ASCII file
%  with extensions of .avg and .mrk, respectively.
% Usage: bt_ascii2avg('filename.txt')

% Dependencies: eeg_load eeg_write ms2row2

BaseFilename = ASCIIfile(1:end-4);

% (1) LOCATE EACH HEADER & DATASET IN ASCII FILE
% Locate "marker report" using the fact 'Chan 1 Bank 1',  appears in column
%         1 of the marker report
% open file
fid = fopen(ASCIIfile);
       
x = 1;  % index for the number of replications sets found
y = 1;  % index for the number of headers found
i = 0;  % index for line number, increments before each line was read.
        % wasn't able to find function that returned line number
        % automatically. 
while(~feof(fid));   % go line by line until the end of the file (feof) is reached.
    i = i+1;   % increment line number
    lyne{i} = fgetl(fid);  % read "line"
    % find marker report
    if (strfind(lyne{i},'Chan 1 Bank 1'));  % search for all instances of 'Chan1 Bank1';
        LatencyLine = i ;                   % returns line number if a match is found
        add_start = i+2;                    % match+2 = beginning of data
    end
end
fclose(fid); % done using textscan, can now close file.

% Extract latencies and save .mrk file
Latencies = lyne(1,LatencyLine);
Latencies = cell2mat(Latencies);
Latencies(Latencies == '"') = [];
Latencies = str2double(Latencies)';
Latencies = Latencies(5:end);
dlmwrite([BaseFilename '.mrk'],Latencies'); % comma delimited

clear lyne i x y fid LatencyLine header_start

% Save .avg file
templateFILE = eeg_load('BioMARK-Template.avg'); % load template file
pnts = templateFILE.header.pnts;  % number of points to extract
add = dlmread(ASCIIfile, '\t', [add_start-1 0 (add_start+pnts-2) 0]);

% prestim baseline "ADD file"
baselined = add - mean(add(ms2row2(templateFILE, ...
    templateFILE.header.xmin*1000):ms2row2(templateFILE, 0)));
templateFILE.data.samples = baselined;
%templateFILE.data.samples = add;
% write baselined data to file.
eeg_write(templateFILE, [BaseFilename '.avg']);  % generate avg file
