function [Picked Auto] = bt_peaks(FilenameBase)
% 

% Dependancies: openavg, ms2row, closestrc

% read in data file
y = openavg([FilenameBase '.avg']);
Baseline = mean(y.signal(ms2row(y,y.xmin):ms2row(y,0)));

% There will be two sets of latency/amplitude/VAcomplex values returned.
% 1) As picked (in .mrk file)  Variable prefix: Picked
% 2) Auto-picked, using .mrk file as seed. Variable prefix: Auto

%%%% AS PICKED %%%%
% read in latency file
PickedLats = dlmread([FilenameBase '.mrk'])';
PickedAmps = y.signal(ms2row(y,PickedLats));

% V/A calcs: assumes V and A are first 2 peaks
Picked_VA_interpeak_lat = PickedLats(2)-PickedLats(1);
Picked_VA_interpeak_amp = PickedAmps(1)-PickedAmps(2);
Picked_VA_slope = -Picked_VA_interpeak_amp./Picked_VA_interpeak_lat;
Picked_VA_area = y.signal-PickedAmps(2); % waveform shifted to "A" = 0 uV
Picked_VA_area = sum(Picked_VA_area(ms2row(y,PickedLats(1)):...
    ms2row(y,PickedLats(2)))); % sum points btw V and A
Picked_VA_area = Picked_VA_area*((y.xmax-y.xmin)/y.pnts); % convert 
        % from µV*pt to µV*ms

%%%% AUTOPICKED %%%%
NumPeaks = length(PickedLats);
AutoLats = zeros(NumPeaks,1); % preallocate for loop below
AutoAmps = zeros(NumPeaks,1); % preallocate for loop below
Step = 1000/y.rate;

% For wave V (lat(1)), find max amplitude -0/+2 points from picked
points = y.signal(ms2row(y,PickedLats(1)):ms2row(y,PickedLats(1))+2);
AutoAmps(1) = max(points);
adjustment = (find(points==AutoAmps(1)));
% sometimes, more than one point has same amp; choose closest to picked
if length(adjustment)>1
    [junk row junk] = closestrc(adjustment,1); % find nearest to picked
    adjustment = adjustment(row);
end
adjustment = (adjustment-1).*Step;
AutoLats(1) = PickedLats(1)+adjustment;

% For all other peaks, find min amplitude -2/+2 points from picked
for x=2:NumPeaks
    points = y.signal(ms2row(y,PickedLats(x))-2:ms2row(y,PickedLats(x))+2);
    AutoAmps(x) = min(points);
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
Auto_VA_interpeak_lat = AutoLats(2)-AutoLats(1);
Auto_VA_interpeak_amp = AutoAmps(1)-AutoAmps(2);
Auto_VA_slope = -Auto_VA_interpeak_amp./Auto_VA_interpeak_lat;
Auto_VA_area = y.signal-AutoAmps(2); % waveform shifted to "A" = 0 uV
Auto_VA_area = sum(Auto_VA_area(ms2row(y,AutoLats(1)):...
    ms2row(y,AutoLats(2)))); % sum points btw V and A
Auto_VA_area = Auto_VA_area*((y.xmax-y.xmin)/y.pnts); % convert 
        % from µV*pt to µV*ms
        
% Consolodate into two cell arrays, correcting amps to prestim baseline
Auto{1} = AutoLats;
Auto{2} = AutoAmps-Baseline;
Auto{3} (1,1) = Auto_VA_interpeak_lat;
Auto{3} (2,1) = Auto_VA_interpeak_amp;
Auto{3} (3,1) = Auto_VA_slope;
Auto{3} (4,1) = Auto_VA_area;

Picked{1} = PickedLats;
Picked{2} = PickedAmps-Baseline;
Picked{3} (1,1) = Picked_VA_interpeak_lat;
Picked{3} (2,1) = Picked_VA_interpeak_amp;
Picked{3} (3,1) = Picked_VA_slope;
Picked{3} (4,1) = Picked_VA_area;

