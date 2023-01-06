function ihsascii2avg(Filename,TemplateFilename);
% Converting IHS ASCII-exported data files to .avg format.
% Assumes a one-channel recording in the following format:
%
% Bunch of header lines
% Marked points, if applicable
% Then, the data arranged in eight columns, with column 2 being
%  latency (in ms), and columns 5 and 7 being the two polarity waveforms
%
% Requirements:
% A template .avg file called IHS-template.avg
% eeg_load.m, eeg_write.m
% 2nd argument is optional.  Will default to 'IHS-template.avg' in current
% folder if omitted.

if nargin < 2
    TemplateFilename = 'IHS-Template.avg';
end

% Last header line begins with "Data" so we read in headers until found.
% While at it, we pick up line numbers of sweeps and gain.
line = 0;        
fid = fopen(Filename);
while(~feof(fid));
    line = line+1;
    Header{line} = fgetl(fid);
    if (strfind(Header{line},'Sweep'))
        Sweepline = line;
    end
    if (strfind(Header{line},'Amp'))
        Gainline = line;
    end
    if (strfind(Header{line},'Data'))
        break
    end
end
fclose(fid);

% Load in data from ascii file:
% Note, the line number to start import is "line" not "line+1"
% because dlmread is 0-indexed
Data = dlmread(Filename,',',line,0);
clear ans fid line 

% Load in template file, with four pre-labeled channels
T = eeg_load(TemplateFilename);

% insert two polarites in 1, 2.  We are picking up the A/D sample columns
% because the uV conversion in the ascii files are very rough.  We will
% convert to uV below. (channels 3 % 4, add/sub, will come later after uV
% conversion)
T.data(1,1).samples = Data(:,5);
T.data(1,2).samples = Data(:,7);
T.variance(1,1).samples = zeros(length(Data),1);
T.variance(1,2).samples = zeros(length(Data),1);
T.variance(1,3).samples = zeros(length(Data),1);
T.variance(1,4).samples = zeros(length(Data),1);

% Now take care of sampling rate, start/stop points, and other things of
% that nature.  Mostly relying on "Data" and math rather than Header info
% except for #sweeps and gain
T.header.pnts = length(Data);
T.header.rate = round(length(Data)./((Data(end,2)-Data(1,2))./1000));
T.header.xmin = Data(1,2)/1000;
T.header.automin = Data(1,2)/1000;
T.header.displayxmin = Data(1,2)/1000;
T.header.xmax = Data(end,2)/1000;
T.header.automax = Data(end,2)/1000;
T.header.displayxmax = Data(end,2)/1000;

% Pull # of sweeps out of text line using reg exp (no I don't understand
% reg exps)
[junk junk junk NumSwps] = regexp(Header{Sweepline}, '(\d)*(\.)?(\d)*$');
NumSwps = str2double(NumSwps);
T.header.N = NumSwps;
T.header.compsweeps = NumSwps;
T.header.acceptcnt = NumSwps;

% Pull gain out of text line using reg exp (no I don't understand
% reg exps)
[junk junk junk Gain] = regexp(Header{Gainline}, '(\d).(\d)');
Gain = str2double(Gain).*1000; % expressed in "k" must mult by 1000.

% Now, convert to µV.
% Some constants
VFS=10; % volts full scale
ADRange = 2^16/2; % 16 bit sampling (signed, so half)

% perform A/D sampling integer to uV conversion.  Note that for channels 1
%  and 2 (the 2 polarities) the total #sweeps must be divided by 2 for the
%  #sweeps of each polarity
T.data(1,1).samples = (T.data(1,1).samples.*1e6.*VFS)./(ADRange.*(NumSwps./2).*Gain);
T.data(1,2).samples = (T.data(1,2).samples.*1e6.*VFS)./(ADRange.*(NumSwps./2).*Gain);
%create add/sub channels
T.data(1,3).samples = T.data(1,1).samples + T.data(1,2).samples;
T.data(1,4).samples = T.data(1,1).samples - T.data(1,2).samples;

% save to .avg file
[Pathy Rooty] = fileparts(Filename);
if isempty(Pathy)==1
    Outfile = [Rooty '.avg'];
else
    Outfile = [Pathy '\' Rooty '.avg'];
end
eeg_write(T, Outfile);