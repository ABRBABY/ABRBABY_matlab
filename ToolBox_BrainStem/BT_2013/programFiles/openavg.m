function [avg] = openavg(FILENAME)
% openavg: Reads NeuroScan AVG File (version 4.2 and 4.3, and earlier)

% UPDATE    version remark
% 1062001   0.1     primitive version based on a C program 
% 1092001   1.0     fully working version based on loadegg.m that I programmed 
% 1102001   1.1     adding channel names loading
% 03/2002   1.2     Darren.Weber@flinders.edu.au
%                   - S_nsweeps_offset from 362 to 364;
%                     so that it finds ACCEPTED not TOTAL SWEEPS, which
%                     has a great impact on conversion to uV.
%                   - modified xmin/xmax to msec from sec.
%                   - modified output to structure with fields
%                   - Changed the name of the function (from loadavg).
% 02/2004   1.3     eeskoe
%                   nsweeps_offset from 364 back to 362;
% 05/2004   1.4     eeskoe
%                   msweeps is initially read from line 362. If 362
%                   returns 1, the program continues. Else, nsweeps is
%                   re-read from 364. 
% 12/2004  1.5      nsweeps is initally read from line 362. If 362 is greater than 1
%                   then check value of line 361. If 362 equals 361, then S_nsweeps is read from line 362.
%                   If 362 and 361 are different then nsweeps is re-read from line 364. 
%                   For grand average files line 361 = 362, and this value equals the number of files that went
%                   into the grandaverage.
% 1/10/2005 1.6     In version 1.4, avg.nsweeps was being read from the wrong line.
%                   Corrected to:  avg.nsweeps is now being read from line
%                   364. This is the value of the accepted sweeps.
%
%                   In scaling the signal, three values are important:
%                   sweeps (line 360), compsweeps (line 362), acceptcnt
%                   (line 364).
%
%                   If line 360 is greater than 1 then check value of line 362. 
%                   If 360 equals 362, scale signal using line 362. If 360
%                   is not equal to 362, scale using line 364.
%                   

%%% Daren Webber's Modifications: 
% This is a modified version of the cumbersomely titled 'eeg_load_scan3avg' function.
% 
% Other modifications:
% 1. chan_names is now a character array
% 2. signal and variance are now formatted (points,channels) rather than (channels,points)
%    that is, equivalent to the Neuroscan's rows=points exporting option

%%% eeskoe's modifications
%%%%  version 1.3 %%%%%
% 1. The number of sweeps in the average is read from line 262 (modified
%    from 264).  This corrects previous signal scaling problems between scan4.2 files, and earlier versions. (see
%    more detailed explanation, line 90).
%%%% version 1.4 %%%%%%%
%~ .  The corrects signal scaling problems for Scan 4.1,
%     4.2, 4.3 files.  The trend so far seems to be:
                % 4.1 or earlier.  For some 4.1 files, the TOTAL nsweeps is
                % reported as 1 (line 362). For other 4.1 files it is not,
                % and in these cases 362 and 364 return the same value. 
       
                % 4.2 and 4.3 files are scaled based on the ACCEPTED number
                % of sweeps (line 364).  For the case of 4.2 files, line
                % 362 and 364 (i.e. accepted and total nsweeps), are
                % usually very similar in value. Whereas, for 4.3 they are
                % not. 
%%%% version 1.5 %%%%%%%
% Modification require to correctly scale grand average files.

%%%% version 1.6
%  Modified 2/28/05 Erika Skoe, eeskoe@northwestern.edu 
% 	(1) Corrected xmax so that it corresponds to what Neuroscan) and Biologic
% 	are reports.  Original openavg reported xmax+1 sample point. This problem was corrected in the original version of openavgdBioMARK in 2005
% 	but was never incorporated into openavg
% 	(2) Biologic uses non-integer sample rates.  openavg reads in only integer values for sample rates. Consequently, for BioMARK files the sample rate
% 	is recalculated based on the time window and number of points. This
% 	will only happen if there is one channel, and that channel is named "BioMARK" , or has "BioMARK in the name".
% 	This recalculation happens after (1).  (total time window = ms of last point-ms of first point+1 sample point).
%   (3) GFP AND REF channels scaled by either compsweeps (if grand
%   	average) or acceptsweeps (if single average). However,  Neuroscan doesn't
%   	actually allow you to create average GFPs. 
%	(4) openavg and openavgBioMARK are now identical.
     


% USEAGE:  [avg] = openavg(FILENAME)
%
%   FILENAME     input Neuroscan .avg file
%   avg          output data structure, with fields:
%   
%   avg.signal      - ERP signal (uV, Npnts x Mchan)
%   avg.variance    - variance of the signal (Npnts x Mchan)
%   avg.chan_names  - electrode labels
%   avg.pnts        - number of points in ERP waveform
%   avg.rate        - sample rate (Hz)
%   avg.xmin        - prestimulus epoch start (e.g., -100 msec)
%   avg.xmax        - poststimulus epoch end (e.g., 900 msec)
%   avg.nsweeps     - number of accepted trials/sweeps in avg
%   
%   e.g.
%   avg = openavg( 'test.avg' );
%   plot( avg.signal );
%
% This program is distributed under the GNU GPL; you can redistribute 
% it and/or modify it. This program is distributed in the hope that it 
% will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% Version 1.1, arno_delorme@salk.edu
% Version 1.2, Darren.Weber@flinders.edu.au
%
% Average data is stored as 4-byte floats in vectored format for each
% channel. Each channel has a 5-byte header that is no longer used. Thus,
% after the main file header, there is an unused 5-byte header followed by
% erp.pnts of 4-byte floating point numbers for the first channel; then a
% 5-byte header for channel two followed by erp.pnts*sizeof(float) bytes,
% etc. Therefore, the total number of bytes after the main header is:
% erp.nchannels * (5 + erp.pnts*sizeof(float)). 

% To scale a data point to microvolts, multiply by the channel-specific calibration factor (i.e., for
% electrode j: channel[j]->calib) and divide by the number of sweeps in the
% average (i.e., channel[j]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1,
help openavg; return;
end;
fid = fopen(FILENAME,'r','ieee-le');
if fid<0,
    msg = sprintf('openavg: Cannot find:\n... %s\n', FILENAME);
    error(msg);
end;
BOOL='int16';
ULONG='int32'; 
FLOAT='float32';
% read # of channels, # of samples, variance flag, and time bounds
% ----------------------------------------------------------------
S_sweeps_offset    = 360; 
S_compsweeps_offset = 362;  % added 12/17/2004
S_acceptcnt_offset =  364;
S_pnts_offset       = 368;
S_nchans_offset     = 370;
S_variance_offset   = 375;
S_rate_offset       = 376;
S_xmin_offset       = 505;
S_xmax_offset       = 509;
packed_sizeof_SETUP = 900;

% the variable names 'nsweeps', 'compsweeps', and 'acceptcnt' are taken
% from eeg_load.  nsweeps is not equivalent to avg.nsweeps.  avg.nsweeps is
% the total number of accepted sweeps (acceptcnt).  'nsweeps' is also NOT
% the total number of sweeps (accepted + rejected).  
% The signal is scaled based on the values of nsweeps, compsweeps, and
% acceptcnt.


fseek(fid, S_sweeps_offset, 'bof');    sweeps = fread(fid, 1, 'ushort');
fseek(fid, S_compsweeps_offset, 'bof'); compsweeps = fread(fid, 1, 'ushort');
fseek(fid, S_acceptcnt_offset, 'bof');  acceptcnt = fread(fid, 1, 'ushort');


fseek(fid, S_pnts_offset, 'bof');       avg.pnts = fread(fid, 1, 'ushort');
fseek(fid, S_nchans_offset, 'bof');     chan = fread(fid, 1, 'ushort');
fseek(fid, S_variance_offset, 'bof');   variance_flag = fread(fid, 1, 'uchar');
fseek(fid, S_rate_offset, 'bof');       avg.rate = fread(fid, 1, 'ushort');
fseek(fid, S_xmin_offset, 'bof');       avg.xmin = fread(fid, 1, 'float32') * 1000;
fseek(fid, S_xmax_offset, 'bof');       avg.xmax = fread(fid, 1, 'float32') * 1000;
fseek(fid, packed_sizeof_SETUP, 'bof');

% Erika's modification: The function fprintf writes text to the screen. 
% To speed up execution, and to free up memory, all lines containing fprint have 
% commented out.

% fprintf('number of channels : %d\n', chan);
% fprintf('number of points   : %d\n', avg.pnts);
% fprintf('sampling rate (Hz) : %f\n', avg.rate);
% fprintf('xmin (msec)        : %f\n', avg.xmin);
% fprintf('xmax (msec)        : %f\n', avg.xmax);
% fprintf('Accepted sweeps    : %d\n', avg.nsweeps);

% read electrode configuration
% ----------------------------

% fprintf('Electrode configuration\n');
for elec = 1:chan,
    channel_label_tmp = fread(fid, 10, 'uchar');
    avg.chan_names(elec,:) = channel_label_tmp';
    for index = 2:9,
        if avg.chan_names(elec,index) == 0,
            avg.chan_names(elec,index) = ' ';
        end;
    end;
    erp = fread(fid, 47-10, 'uchar');
    baseline(elec) = fread(fid, 1, 'ushort');
    erp = fread(fid, 10, 'uchar');
    sensitivity(elec) = fread(fid, 1, 'float32');
    erp = fread(fid, 8, 'uchar');
    calib(elec) = fread(fid, 1, 'float32');
    %fprintf('%s: baseline: %d\tsensitivity: %f\tcalibration: %f\n', avg.chan_names(elec,1:4), baseline(elec), sensitivity(elec), calib(elec));
    factor(elec) = calib(elec) * sensitivity(elec) / 204.8;
end;
% Read signal data (amplifier units)
signal = zeros(avg.pnts, chan);
for elec = 1:chan,
    fseek(fid, 5, 'cof'); % skip sweeps header
    signal(:, elec) = fread(fid, avg.pnts, 'float32');
end;
if variance_flag,
    variance = zeros(avg.pnts, chan);
    for elec = 1:chan,
        variance(:, elec) = fread(fid, avg.pnts, 'float32');
    end;
    avg.variance = variance;
else
    avg.variance = [];
end;

% Convert signal to microvolts
baseline = repmat(baseline,avg.pnts,1);
calib    = repmat(calib,   avg.pnts,1);


% eeskoe modification 12/16/2004
% if nsweeps is not equal to one, then check whether nsweeps = compsweeps. 
if compsweeps ~= 1
    
    if compsweeps ==  sweeps               % grand avg files
 
        signal  = (signal - baseline) .* calib/compsweeps;
       
    else  

        signal = (signal - baseline).* calib/acceptcnt;
       
    end
end

% scale GFP and REF files.
avg.chan_names=char(avg.chan_names);
for i = 1:size(avg.chan_names,1);
    if findstr('GFP', avg.chan_names(i,:))

        matchGFP = i;
        if compsweeps ~= 1
            if compsweeps == sweeps               % grand avg files
               
                signal(:, matchGFP) = signal(:,matchGFP).*compsweeps;
               

            else
               
               signal(:,matchGFP) = signal(:,matchGFP).*acceptcnt;
               

            end
        end
        break
    end
end

for i = 1:size(avg.chan_names,1);

    if findstr('REF', avg.chan_names(i,:))
        matchREF = i;
      
         if compsweeps ~= 1
            if compsweeps == sweeps               % grand avg files

                signal(:,matchREF) = signal(:,matchREF).*compsweeps;
               

            else

               signal(:,matchREF) = signal(:,matchREF).*acceptcnt;
               

            end
            
         end
         break
         end
end
 
      





% number of sweeps is read from line 364. This value is what gets reported.
fseek(fid, S_acceptcnt_offset, 'bof');  avg.nsweeps = fread(fid, 1, 'ushort');

avg.signal = signal;

% need to correct xmax so that it corresponds to what Biologic and
% Neuroscan are reporting.
if size(avg.chan_names,1)==1
    if strfind(avg.chan_names, 'BioMARK');
        avg.rate = avg.pnts/((avg.xmax-avg.xmin)/1000);
    end
    
    if strfind(avg.chan_names, 'BioMAP');
        avg.rate = avg.pnts/((avg.xmax-avg.xmin)/1000);
    end
end
    
 avg.xmax = avg.xmin+((avg.pnts-1)*(1/avg.rate)*1000);     



fclose(fid);

% assignin('base','avg', avg);  
% return;
