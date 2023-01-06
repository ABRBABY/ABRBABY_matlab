function [StraightCorr CrossCorr Lag]=bt_qncorr(Qfile, Nfile, start, stop, chan)
% bt_qncorr compares two files via correlation.  Typically this is used to compare a response recorded 
% in background noise with a response to the same stimulus in quiet.
%
% Usage: [StraightCorr CrossCorr Lag]=bt_qncorr(Qfile, Nfile, start, stop);
%
% Requires 4 input arguments:
% Qfile, Nfile:  quiet and noise .avg files.
% start, stop: latency, in ms., over which the correlations are run
%
% 3 variables are returned to the workspace: 
% (1) StraightCorr: Simple Pearson's r over desired latency range
% (2) CrossCorr: Maximum correlation, allowing noise response to lag quiet by
%     up to 2 ms.
% (3) Lag: Amount of timeshift (in ms.) to achieve CrossCorr
%       In cases where the best correlation occurs with no shift, Lag will be 0
%       and StraightCorr will equal CrossCorr.
%
% Dependancies: 
%   openavg.m
%   xcorrelation.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perform "straight" correlation and cross-correlation, allowing noise to
% lag quiet by up to two ms.  Using bt_xcorrelation to perform both, though
% the former could be accomplished with corrcoef.
[junk  StraightCorr] = bt_xcorrelation(Nfile, Qfile, start, stop, 0, 0, 'POSITIVE', chan);
[Lag  CrossCorr] = bt_xcorrelation(Nfile, Qfile, start, stop, 0, 2, 'POSITIVE', chan);
