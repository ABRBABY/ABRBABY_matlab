function [Picked Auto] = bt_peaks3(AvgFile,Markers,chan)
% Find lats, amps, etc of up to seven peaks/troughs. 
%    Inputs are the avg file and an X x 2 input where
%    col1 is the (up to) 10 lats and col2 is polarity (1 = pos).
% Can be used as a stand-alone, but designed to be used with bt_gui.  

% Dependancies: openavg, ms2row, closestrc

% read in data file
y = openavg(AvgFile);
y.signal = y.signal(:,chan);
Baseline = mean(y.signal(ms2row(y,y.xmin):ms2row(y,0)));

% There will be two sets of latency/amplitude/VAcomplex values returned.
% 1) As picked (in .mrk file)  Variable prefix: Picked
% 2) Auto-picked, using .mrk file as seed. Variable prefix: Auto

%%%% AS PICKED %%%%
% read in latency file
PickedLats = Markers(:,1);
PickedAmps = y.signal(ms2row(y,PickedLats));

% Execute only if there any marks to read; allows for null
if sum(sum(PickedLats))
    %%%% AUTOPICKED %%%%
    NumPeaks = length(PickedLats);
    AutoLats = zeros(NumPeaks,1); % preallocate for loop below
    AutoAmps = zeros(NumPeaks,1); % preallocate for loop below
    Step = 1000/y.rate;

    % Find min or max amplitude -2/+2 points from picked
    for x=1:NumPeaks
        points = y.signal(ms2row(y,PickedLats(x))-2:ms2row(y,PickedLats(x))+2);
        
        % seek peak or trough depending on Markers(:,2)
        if Markers(x,2) == 1
            AutoAmps(x) = max(points);
        else
            AutoAmps(x) = min(points);
        end
        adjustment = (find(points==AutoAmps(x)));
        
        % sometimes, more than one point has same amp; choose closest to picked
        if length(adjustment)>1
            [junk row junk] = closestrc(adjustment,3); % find nearest to picked
            adjustment = adjustment(row);
        end
        adjustment = (adjustment-3).*Step;
        AutoLats(x) = PickedLats(x)+adjustment;
    end

    % V/A calcs: assumes V and A are first 2 peaks
  
    % Consolodate into two cell arrays, correcting amps to prestim baseline
    Auto{1} = AutoLats;
    Auto{2} = AutoAmps-Baseline;

    Picked{1} = PickedLats;
    Picked{2} = PickedAmps-Baseline;

else
    Picked = cell(1);
    Auto = cell(0);
end
