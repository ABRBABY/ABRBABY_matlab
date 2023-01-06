function [PosCorr PosLag NegCorr NegLag] = bt_srcorr(file, start, stop)
% Compares BioMaRK file to the evoking stimulus using cross correlation
% Usage: [PosCorr PosLag NegCorr NegLag] = bt_srcorr('file.avg', start, stop);
%
% A typical BioMARK response, when cross-correlated with the stimulus,
%  exhibits two major correlogram maxima: one positive, one negative. 
%  This function returns both and their respective lags.
%
% INPUT ARGUMENTS:
% 'file.avg': BioMARK data file (in .avg format)
% start, stop: latency, in ms., over which the correlation is run
%
% OUTPUT ARGUMENTS:
% PosCorr, NegCorr: maximum/minimum r-values found at lags corresponding
%   to onset response.
% PosLag, NegLag: lags corresponding to above.
%
% Dependancies: 
%   filtered_stimulus.avg
%   openavg.m
%   xcorrelation.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open stimulus file
stim = 'filtered_stimulus.avg';

% Within selected range find max correlation (positive or negative)
%   over the following 2 ranges:
% (a) positive corr between 6.6 & 9.6 ms
% (b) negative corr between 5.8 & 7.8 ms
[PosLag PosCorr] = xcorrelation(file, stim, start, stop, 6.6, 9.6, 'POSITIVE');
[NegLag NegCorr] = xcorrelation(file, stim, start, stop, 5.8, 7.8, 'NEGATIVE');





    
    
    
