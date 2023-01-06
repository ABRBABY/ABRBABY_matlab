function wav2avg(filename, fs)
% Description:
% Takes sound files (.wav) and converts them to Neuroscan (.avg) files.
%   Output avg file only has one (left) channel if wav was stereo recording
%   Output resampled to desired rate
%
% Usage: wav2avg('filename.wav', fs)
%   fs is the target sampling rate
%
% Example: wav2avg('afile.wav', 10000)
%
% Dependencies wavread, resample, writeavg

[y OldFS] = wavread(filename);
y = y(:,1);
newy = resample(y,fs,OldFS);

%Create for avg file
f.nsweeps = 1;
f.pnts = length(newy);
f.rate = fs;
f.xmin = 0;
f.xmax = (length(newy)/fs)*1000;
f.chan_names = 'stimulus  ';
f.variance = [];
f.signal = newy;

writeavg (f, [filename(1:end-4) '.avg']);