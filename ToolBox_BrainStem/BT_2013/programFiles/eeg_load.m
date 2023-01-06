function [f,fid] = eeg_load(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As of Feb 2008, eeg_load and eeg_loadBioMARK are identical. At this point, both files will be kept 
%  so that programs continue to run without hiccups. (eeskoe, Feb 28, 2008).
%
%Modified version of eeg_load_scan41
%
% Load a Neuroscan 'scan4.1, scan4.2, scan4.3' average file
%
% Useage:   [f,fid] = eeg_load(filename)
%
% where:    filename is a complete 'path\fileprefix.ext' string;
%           
%
%           f is a structure containing:
%           
%               f.header
%               f.electloc
%               f.data
%               f.variance
%               f.tag
%               f.chan_names
%
%           fid is a file handle to 'filename'
%
% Author:   Sean Fitzgibbon <psspf@id.psy.flinders.edu.au>
% Created:  08/2000 (load scan4.1 file only)
% Modified: 07/2001, Darren.Weber@flinders.edu.au
%                    tab formatted this file and added help

% Modified: 7/23/2004 Erika Skoe eeskoe@yahoo.com
%                     modified to load scan4.2 and scan4.3 files
%                    
% Average data is stored in terms of amplifier units.  To scale a data point to microvolts, divide by the number of sweeps in the
% average which is either h.compsweeps (total sweeps) OR h.acceptcnt (accepted sweeps).  If the value of h.compsweeps equals 1, then the signal (i.e.
% f.data.samples) is already scaled into uV. If the value is something
% other than 1, the signal must be scaled into uV by dividing the signal by
% the the value of h.acceptcnt.

% Modified: 10/01/2004  Erika Skoe eeskoe@northwestern.edu
%                       Added code to extract channel names from file.
%                       added f.chan_names to f.
% Modified: 12/17/2004 Erika Skoe eeskoe@northwestern.edu
%                      Modified to load grand average files.
%
% To produce a properly scaled waveform, the signal must be divided by the value of h.compsweeps and not h.acceptcnt. 
% Unlike, other signals examined thus far, h.nsweeps = h.compsweeps for grand avg files.  For these files, h.nsweeps
% refers to the number of files that were averaged to make the grandavg.
% 
% Modified: 1/10/05  Erika Skoe eeskoe@northwestern.edu
%                    Signal is scaled using the baseline and calibration
%                    values of each electrode. 
%                    signal = ((signal - baseline)*calibration factor)./sweeps
%                    For all of our LLB files survey thus far, baseline = 0, and calib = 1;
%                     
% Modified 2/28/05 Erika Skoe, eeskoe@northwestern.edu 
% 	(1) Corrected xmax so that it corresponds to what Neuroscan and Biologic
% 	reports  Original eeg_load reported xmax+1 sample point. (This problem was corrected in the original version of eeg_loadBioMARK in 2005
% 	but was never incorporated into eeg_load)
%
% 	(2) Biologic uses non-integer sample rates.  eeg_load reads in only integer values for sample rates. Consequently, for BioMARK files the sample rate
% 	is recalculated based on the time window and number of points. This will only happen if there is one channel, and that channel is named "BioMARK".
% 	This recalculation happens after (1).  (total time window = ms of last point-ms of first point+1 sample point).
%
%       (3) GFP AND REF channels are now properly scaled. Unlike the other channels, GFP and REF are not scaled by  h.compsweeps (if grand average) or h.acceptsweeps (if single average).
%
%  	(4) Marker information is extracted,  including latency (in pts and ms), amplitude marker name, and channel.
%   	Added f.marker to f (original author Jade Wang jadewang@gmail.com)
%   	Modified Jade's code so that amplitude is reported

%  	 For example:  
% 		single-channel file, two markers:   f.marker = {[1x1 struct]  [1x1 struct]}
% 		f.marker{x} = 	name: 'New                 '
%     				at_pts: 623
%     				latency: 2.12
%                      		amplitude: -1.5
%                      		channel: 'CZ-A2     '
% 	 			ch_num: 1
% 	 The order of the marker is based on the order that the marks were added.  The marks do not get sorted by channel name, or latency at this point.
%       * If there are no markers, markers:   marker = No Markers'
% 	
%    	%known bug:  Marker name cannot be numeric values.  If numeric, the name will be reported as 'New'.
%
%	(5) eeg_load and eeg_loadBioMARK are now identical.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fid = fopen(filename,'r');

h.rev               = fread(fid,12,'char');
h.nextfile          = fread(fid,1,'long');
h.prevfile          = fread(fid,1,'long');
h.type              = fread(fid,1,'char');
h.id                = fread(fid,20,'char');
h.oper              = fread(fid,20,'char');
h.doctor            = fread(fid,20,'char');
h.referral          = fread(fid,20,'char');
h.hospital          = fread(fid,20,'char');
h.patient           = fread(fid,20,'char');
h.age               = fread(fid,1,'short');
h.sex               = fread(fid,1,'char');
h.hand              = fread(fid,1,'char');
h.med               = fread(fid,20, 'char');
h.category          = fread(fid,20, 'char');
h.state             = fread(fid,20, 'char');
h.label             = fread(fid,20, 'char');
h.date              = fread(fid,10, 'char');
h.time              = fread(fid,12, 'char');
h.mean_age          = fread(fid,1,'float');
h.stdev             = fread(fid,1,'float');
h.n                 = fread(fid,1,'short');
h.compfile          = fread(fid,38,'char');
h.spectwincomp      = fread(fid,1,'float');
h.meanaccuracy      = fread(fid,1,'float');
h.meanlatency       = fread(fid,1,'float');
h.sortfile          = fread(fid,46,'char');
h.numevents         = fread(fid,1,'int');
h.compoper          = fread(fid,1,'char');
h.avgmode           = fread(fid,1,'char');
h.review            = fread(fid,1,'char');
h.nsweeps           = fread(fid,1,'ushort');
h.compsweeps        = fread(fid,1,'ushort');
h.acceptcnt         = fread(fid,1,'ushort');
h.rejectcnt         = fread(fid,1,'ushort');
h.pnts              = fread(fid,1,'ushort');
h.nchannels         = fread(fid,1,'ushort');
h.avgupdate         = fread(fid,1,'ushort');
h.domain            = fread(fid,1,'char');
h.variance          = fread(fid,1,'char');
h.rate              = fread(fid,1,'ushort');
h.scale             = fread(fid,1,'double');
h.veogcorrect       = fread(fid,1,'char');
h.heogcorrect       = fread(fid,1,'char');
h.aux1correct       = fread(fid,1,'char');
h.aux2correct       = fread(fid,1,'char');
h.veogtrig          = fread(fid,1,'float');
h.heogtrig          = fread(fid,1,'float');
h.aux1trig          = fread(fid,1,'float');
h.aux2trig          = fread(fid,1,'float');
h.heogchnl          = fread(fid,1,'short');
h.veogchnl          = fread(fid,1,'short');
h.aux1chnl          = fread(fid,1,'short');
h.aux2chnl          = fread(fid,1,'short');
h.veogdir           = fread(fid,1,'char');
h.heogdir           = fread(fid,1,'char');
h.aux1dir           = fread(fid,1,'char');
h.aux2dir           = fread(fid,1,'char');
h.veog_n            = fread(fid,1,'short');
h.heog_n            = fread(fid,1,'short');
h.aux1_n            = fread(fid,1,'short');
h.aux2_n            = fread(fid,1,'short');
h.veogmaxcnt        = fread(fid,1,'short');
h.heogmaxcnt        = fread(fid,1,'short');
h.aux1maxcnt        = fread(fid,1,'short');
h.aux2maxcnt        = fread(fid,1,'short');
h.veogmethod        = fread(fid,1,'char');
h.heogmethod        = fread(fid,1,'char');
h.aux1method        = fread(fid,1,'char');
h.aux2method        = fread(fid,1,'char');
h.ampsensitivity    = fread(fid,1,'float');
h.lowpass           = fread(fid,1,'char');
h.highpass          = fread(fid,1,'char');
h.notch             = fread(fid,1,'char');
h.autoclipadd       = fread(fid,1,'char');
h.baseline          = fread(fid,1,'char');
h.offstart          = fread(fid,1,'float');
h.offstop           = fread(fid,1,'float');
h.reject            = fread(fid,1,'char');
h.rejstart          = fread(fid,1,'float');
h.rejstop           = fread(fid,1,'float');
h.rejmin            = fread(fid,1,'float');
h.rejmax            = fread(fid,1,'float');
h.trigtype          = fread(fid,1,'char');
h.trigval           = fread(fid,1,'float');
h.trigchnl          = fread(fid,1,'char');
h.trigmask          = fread(fid,1,'short');
h.trigisi           = fread(fid,1,'float');
h.trigmin           = fread(fid,1,'float');
h.trigmax           = fread(fid,1,'float');
h.trigdir           = fread(fid,1,'char');
h.autoscale         = fread(fid,1,'char');
h.n2                = fread(fid,1,'short');
h.dir               = fread(fid,1,'char');
h.dispmin           = fread(fid,1,'float');
h.dispmax           = fread(fid,1,'float');
h.xmin              = fread(fid,1,'float');
h.xmax              = fread(fid,1,'float');
h.automin           = fread(fid,1,'float');
h.automax           = fread(fid,1,'float');
h.zmin              = fread(fid,1,'float');
h.zmax              = fread(fid,1,'float');
h.lowcut            = fread(fid,1,'float');
h.highcut           = fread(fid,1,'float');
h.common            = fread(fid,1,'char');
h.savemode          = fread(fid,1,'char');
h.manmode           = fread(fid,1,'char');
h.ref               = fread(fid,10,'char');
h.rectify           = fread(fid,1,'char');
h.displayxmin       = fread(fid,1,'float');
h.displayxmax       = fread(fid,1,'float');
h.phase             = fread(fid,1,'char');
h.screen            = fread(fid,16,'char');
h.calmode           = fread(fid,1,'short');
h.calmethod         = fread(fid,1,'short');
h.calupdate         = fread(fid,1,'short');
h.calbaseline       = fread(fid,1,'short');
h.calsweeps         = fread(fid,1,'short');
h.calattenuator     = fread(fid,1,'float');
h.calpulsevolt      = fread(fid,1,'float');
h.calpulsestart     = fread(fid,1,'float');
h.calpulsestop      = fread(fid,1,'float');
h.calfreq           = fread(fid,1,'float');
h.taskfile          = fread(fid,34,'char');
h.seqfile           = fread(fid,34,'char');
h.spectmethod       = fread(fid,1,'char');
h.spectscaling      = fread(fid,1,'char');
h.spectwindow       = fread(fid,1,'char');
h.spectwinlength    = fread(fid,1,'float');
h.spectorder        = fread(fid,1,'char');
h.notchfilter       = fread(fid,1,'char');
h.headgain          = fread(fid,1,'short');
h.additionalfiles   = fread(fid,1,'int');
h.unused            = fread(fid,5,'char');
h.fspstopmethod     = fread(fid,1,'short');
h.fspstopmode       = fread(fid,1,'short');
h.fspfvalue         = fread(fid,1,'float');
h.fsppoint          = fread(fid,1,'short');
h.fspblocksize      = fread(fid,1,'short');
h.fspp1             = fread(fid,1,'ushort');
h.fspp2             = fread(fid,1,'ushort');
h.fspalpha          = fread(fid,1,'float');
h.fspnoise          = fread(fid,1,'float');
h.fspv1             = fread(fid,1,'short');
h.montage           = fread(fid,40,'char');
h.eventfile         = fread(fid,40,'char');
h.fratio            = fread(fid,1,'float');
h.minor_rev         = fread(fid,1,'char');
h.eegupdate         = fread(fid,1,'short');
h.compressed        = fread(fid,1,'char');
h.xscale            = fread(fid,1,'float');
h.yscale            = fread(fid,1,'float');
h.xsize             = fread(fid,1,'float');
h.ysize             = fread(fid,1,'float');
h.acmode            = fread(fid,1,'char');
h.commonchnl        = fread(fid,1,'uchar');
h.xtics             = fread(fid,1,'char');
h.xrange            = fread(fid,1,'char');
h.ytics             = fread(fid,1,'char');
h.yrange            = fread(fid,1,'char');
h.xscalevalue       = fread(fid,1,'float');
h.xscaleinterval    = fread(fid,1,'float');
h.yscalevalue       = fread(fid,1,'float');
h.yscaleinterval    = fread(fid,1,'float');
h.scaletoolx1       = fread(fid,1,'float');
h.scaletooly1       = fread(fid,1,'float');
h.scaletoolx2       = fread(fid,1,'float');
h.scaletooly2       = fread(fid,1,'float');
h.port              = fread(fid,1,'short');
h.numsamples        = fread(fid,1,'long');
h.filterflag        = fread(fid,1,'char');
h.lowcutoff         = fread(fid,1,'float');
h.lowpoles          = fread(fid,1,'short');
h.highcutoff        = fread(fid,1,'float');
h.highpoles         = fread(fid,1,'short');
h.filtertype        = fread(fid,1,'char');
h.filterdomain      = fread(fid,1,'char');
h.snrflag           = fread(fid,1,'char');
h.coherenceflag     = fread(fid,1,'char');
h.continuoustype    = fread(fid,1,'char');
h.eventtablepos     = fread(fid,1,'long');
h.continuousseconds = fread(fid,1,'float');
h.channeloffset     = fread(fid,1,'long');
h.autocorrectflag   = fread(fid,1,'char');
h.dcthreshold       = fread(fid,1,'uchar');

for n = 1:h.nchannels
    e(n).lab            = fread(fid,10,'char');
    e(n).reference      = fread(fid,1,'char');
    e(n).skip           = fread(fid,1,'char');
    e(n).reject         = fread(fid,1,'char');
    e(n).display        = fread(fid,1,'char');
    e(n).bad            = fread(fid,1,'char');
    e(n).n              = fread(fid,1,'ushort');
    e(n).avg_reference  = fread(fid,1,'char');
    e(n).clipadd        = fread(fid,1,'char');
    e(n).x_coord        = fread(fid,1,'float');
    e(n).y_coord        = fread(fid,1,'float');
    e(n).veog_wt        = fread(fid,1,'float');
    e(n).veog_std       = fread(fid,1,'float');
    e(n).snr            = fread(fid,1,'float');
    e(n).heog_wt        = fread(fid,1,'float');
    e(n).heog_std       = fread(fid,1,'float');
    e(n).baseline       = fread(fid,1,'short');
    e(n).filtered       = fread(fid,1,'char');
    e(n).fsp            = fread(fid,1,'char');
    e(n).aux1_wt        = fread(fid,1,'float');
    e(n).aux1_std       = fread(fid,1,'float');
    e(n).senstivity     = fread(fid,1,'float');
    e(n).gain           = fread(fid,1,'char');
    e(n).hipass         = fread(fid,1,'char');
    e(n).lopass         = fread(fid,1,'char');
    e(n).page           = fread(fid,1,'uchar');
    e(n).size           = fread(fid,1,'uchar');
    e(n).impedance      = fread(fid,1,'uchar');
    e(n).physicalchnl   = fread(fid,1,'uchar');
    e(n).rectify        = fread(fid,1,'char');
    e(n).calib          = fread(fid,1,'float');
end



for i = 1:h.nchannels
    d(i).header      = fread(fid,5,'char');
    d(i).samples     = fread(fid,h.pnts,'float');
end

%eeskoe modification: 10/01/04 added code to extract channel names
% channel names start on line 900.
startchan_names = 900;
fseek(fid, startchan_names, 'bof');   % 'bof' = beginning of file.

for j = 1:h.nchannels
    channel_label_tmp = fread(fid, 10, 'uchar');
    chan_names(j,:) = channel_label_tmp';
    for index = 2:9,
        if chan_names(j,index) == 0,
            chan_names(j,index) = ' ';
        end
    end
    
    startchan_names = startchan_names + 75;  
    fseek(fid, startchan_names, 'bof');
end

for j = 1:h.nchannels
    v(j).samples     = fread(fid,h.pnts,'float');
end


%The trend so far seems to be:
% 4.1 or earlier.  For some 4.1 files, h.compsweeps is
% reported as 1. For other 4.1 files it is not,
% and in these cases h.compsweeps and h.acceptcnt return the same value. 

% 4.2 and 4.3 files are scaled based h.acceptcnt.
% For the case of 4.2 files, h.compsweep and h.acceptcnt, are
%usually very similar in value. Whereas, for 4.3 they are
% not. 

% Grand average files are scaled based on h.compsweeps.
% For these files, h.compsweeps = h.nsweeps

if h.compsweeps ~= 1
    
    if h.compsweeps == h.nsweeps               % grand avg files
        for i = 1:h.nchannels
            d(i).samples = ((d(i).samples-e(i).baseline).*e(i).calib)./h.compsweeps;
            v(i).samples = ((v(i).samples-e(i).baseline).*e(i).calib)./h.compsweeps;
         
        end
    else  
        for i = 1:h.nchannels% 4.2 and 4.3 files
            
            d(i).samples = ((d(i).samples-e(i).baseline).*e(i).calib)./h.acceptcnt;
            v(i).samples = ((v(i).samples-e(i).baseline).*e(i).calib)./h.acceptcnt;

        end
    end

end

% h.nsweeps          
% h.compsweeps       
% h.acceptcnt       
% h.rejectcnt         

for i = 1:h.nchannels
    if strmatch('GFP      ', chan_names(i,:))

        matchGFP = i;
        if h.compsweeps ~= 1
            if h.compsweeps == h.nsweeps               % grand avg files
                
                d(matchGFP).samples = d(matchGFP).samples.*h.compsweeps;
                v(matchGFP).samples = v(matchGFP).samples.*h.compsweeps;

            else


                d(matchGFP).samples = d(matchGFP).samples.*h.acceptcnt;
                v(matchGFP).samples = d(matchGFP).samples.*h.acceptcnt;


            end
        end
       
    end

   
    
    if strmatch('REF      ', chan_names(i,:))
        matchREF = i;
      
        if h.compsweeps ~= 1
            
            if h.compsweeps == h.nsweeps
                             % grand avg files
                   
                    d( matchREF).samples = d( matchREF).samples.*h.compsweeps;
                    v( matchREF).samples = v( matchREF).samples.*h.compsweeps;
            else
                    
                    d( matchREF).samples = d( matchREF).samples.*h.acceptcnt;
                    v( matchREF).samples = v( matchREF).samples.*h.acceptcnt;

                end
            end
      
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
t = fread(fid, 'uchar');
f.header = h;
f.electloc = e;
f.data = d;

f.variance = v;
f.tag = t;
f.chan_names = char(chan_names);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%%%%%
% Need to correct xmax so that it corresponds to what N-SCAN and Biologic
% are reports.  (eeskoe);
if size(f.chan_names,1)==1
    if strfind(f.chan_names, 'BioMARK');
        f.header.rate = f.header.pnts/((f.header.xmax-f.header.xmin));
    end
end
    
f.header.xmax = f.header.xmin+((f.header.pnts-1)*(1/f.header.rate));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% jqw modification:  added code to get marker information from
%%%% t or tag.  markers are 168 elements in length with start/stop
%%%% sequences that are 20 elements in length.  

% parameters
taglen   = 168;
seqlen   = 20;
latindex = 85;
factor   = 256;
chan_idx = 97;
markerID = [80 195 0 0 5 1 0 0 152 0 0 0 152 0 0 0 255 0 0 0]; 
% this is the marker initiation sequence

% part 1: is there a marker? how many markers?

markerstart  = strfind(t',markerID); % how many times and where does it happen?
n_marker     = length(markerstart);
markerexists = (n_marker >= 1);

% part 2: read markers (name and lat)
latency = linspace(f.header.xmin, f.header.xmax, f.header.pnts); %line added by eeskoe; vector of latency values.
if markerexists
    for x=1:n_marker
        start    = markerstart(x);
        myMarker = t(start:start+taglen -1);
        myChan   = myMarker(chan_idx)+1;     % counts channels from 0
                                             % but arrays from 1. :)
                                             % eeg_write counter-compensates
        % marker is cell array
        marker{x}.name    = char(myMarker(seqlen+1:seqlen*2)');
        marker{x}.lat_pts = myMarker(latindex) + factor* myMarker(latindex+1)+1;  % +1 added by eeskoe, original calculations off by one.
        marker{x}.latency = latency(marker{x}.lat_pts);
        marker{x}.amplitude = f.data(myChan).samples(marker{x}.lat_pts);
        marker{x}.channel = char(chan_names(myChan,:));
        marker{x}.ch_num  = myChan;
        
    end
else
    marker = 'No Markers';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f.marker = marker;
frewind(fid);
fclose(fid);