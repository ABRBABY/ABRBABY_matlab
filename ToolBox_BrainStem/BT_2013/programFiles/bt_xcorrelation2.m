function [LAG_atmaxmin, maxmincor, all_corrs, all_lags] = bt_xcorrelation2(filename, comparison, start, stop, lagstart, lagstop, polarity, chan, chancomp)

%% This function calculates the maximum correlation (maxmincor) value and its
%% respective lag (LAG) over an imputted lag range, for a specified portion of a
%% file. The user must specify whether to find the max positive or max
%% negative correlation value

% Description of function arguments:
% 1) filename:  the data of interest
% 2) Comparison:  what you want to correlate the file with 
% 3) start: start latency of the comparison file
% 4) stop: stop latency of the comparison file
% 5) lagstart:   starting lag   (how much the File lags behind the
%                               Comparison)
% 6) lagstop:    stopping lag
% 7) polarity:   If 'POSITIVE', then find max positive correlation value
                % If 'NEGATIVE', then find max negative correlation value

% author: Erika Skoe  eeskoe@yahoo.com
% date:  May 2004

% open file, extract signal, and calculate step size
file = openavg(filename);
filename = file.signal(:,chan);
xmin = file.xmin;
xmax = file.xmax;
pnts = file.pnts;
fs = pnts/(xmax-xmin);
latency = linspace(xmin, xmax, pnts);

% open comparison file
compfile = openavg(comparison);
comparison = compfile.signal(:, chancomp);
latency_comp = linspace(compfile.xmin, compfile.xmax, compfile.pnts);



% (a) Isolate the desired region of the comparison file
startPtComp = ms2row(compfile, start);
startLatencyComp = latency_comp(startPtComp);
stopPtComp = ms2row(compfile, stop);
total = stopPtComp-startPtComp;

% (b) Isolate desired region of the file, this is based on the
% inputted lag ranges.
startLat = lagstart + startLatencyComp;
startPt = ms2row(file, startLat);
startLatency = latency(startPt);
stopPt = startPt + total;
stopLatency = latency(stopPt);
% (c)  find start lag 
REALstartLag = startLatency-startLatencyComp;
msPoints = latency(2)-latency(1);


% (d) Find out how many correlations to perform
totalLagPoints = round((lagstop-lagstart)/msPoints);

% (e) Start doing Correlations
%% Lag increases by msPoints each time through loop
for j = 1:totalLagPoints+1;
    all_corrs(j,1)=nancorrcoef(filename(startPt+(j-1):stopPt+(j-1)), comparison(startPtComp:stopPtComp,1));
    all_lags(j,1) = REALstartLag + (msPoints*(j-1));
end


clear data file stim

% (f)  Find max pos/negative correlation

if strncmpi(polarity,'POSITIVE',3)==1;
    % calculate the max positive correlation 
    maxmincor = max(all_corrs);
    [y, index] = max(all_corrs);
end

if strncmpi(polarity,'NEGATIVE',3)==1;
    maxmincor = min(all_corrs);
    [y, index] = min(all_corrs);
end

%  (g) calculate the LAG
LAG_atmaxmin = all_lags(index,1);


clear data file stim  startLatencyComp startPtComp total stopPtComp startLat startLatency startPt stopPt
clear stopLatency totalLagPoints y polarity start stop lagstart lagstop REALstartLag filename comparison compFile
clear xmin xmax