function [RMS RMSprestim SNR] = bt_rms(FILENAME,start,stop, chan)
% Computes RMS amplitude of response period (selectable), of prestim
% period, and the resulting SNR.
%
% Usage: [RMS RMSprestim SNR] = bt_rms('filename.avg',start_latency,stop_latency);
%    Note, if start,stop omitted, BioMARK default latencies (11.38, 40.58) are used.
% Three variables returned to workspace are 
% RMS: RMS amplitude of selected latency range
% RMSprestim: RMS amplitude of prestimulus activity
% SNR: The quotient of RMS/RMSprestim

% accessory m-files needed: ms2row, openavg

if nargin < 3
    disp('Using default FFR period of 11.38-40.58 ms.');
    start = 11.38;
    stop = 40.58;
end

avg = openavg(FILENAME);
data = avg.signal(:,chan);

% Note: std(X,1) is a shortcut RMS calc. It is equivalent to RMS
%   on a baselined (demeaned to 0) waveform.  This is what we want here.

% response period
Segment = data(ms2row(avg,start):ms2row(avg,stop));
RMS = std(Segment,1); 

% Prestim period
Segment = data(ms2row(avg,avg.xmin):ms2row(avg,0));
RMSprestim = std(Segment,1);

% signal-to-noise
SNR = RMS./RMSprestim;
