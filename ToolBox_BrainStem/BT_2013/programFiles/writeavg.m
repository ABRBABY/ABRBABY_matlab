function writeavg(f,filename)

% writeavg - Write a Neuroscan average file
%
% Useage:   writeavg(f,filename)
%
% where:    filename is a complete 'path\fileprefix.ext' string;
%           f is a structure in "openavg" format:
%           
%               f.nsweeps
%               f.pnts
%               f.rate
%               f.xmin
%               f.xmax
%               f.chan_names
%               f.variance
%               f.signal
%
%   See companion function: openavg.m
%
%   Original function: eeg_load_scan41.m
%   See also: eeg_load_scan41.m

%  Revision: 2.0:  companion to openavg.m
%
% Licence:  GNU GPL, no implied or express warranties
% Created:  08/2000, Sean Fitzgibbon <psspf@id.psy.flinders.edu.au>
% Modified: 08/2001, Darren.Weber@flinders.edu.au
%                    distribute under GPL, with Sean's permission
% Modified: 06/2004  Copiously changed to write files from structures
%                    created with "openavg.m" which is a pared down
%                    .avg file reader.  openavg ignores most of the
%                    header information and only brings in the most
%                    crucial bits, as outlined above.
%                    Trent Nicol, tgn@northwestern.edu
% Modified 03/2008   To be read correctly by Neuroscan, f.xmax must be set to xmax (display value) +1;
%                    To accomondate this, f.xmax is directly calculated from f.xmin, pnts and sample rate (first line of code). 
%                    calculated from the sample rate, f.xmin and the number of points.  
%                    If the calculated value does not equal f.xmax+1pt, a warning message is displayed.
%                    eeskoe@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2.xmax = f.xmin+((f.pnts)*(1/f.rate)*1000);   %% added by eeskoe (2/2008)

if f2.xmax ~=f.xmax+((1/f.rate*1000))
    display(['warning: xmax adjusted to ',  num2str(f2.xmax-((1/f.rate)*1000))]);
else
    f.xmax = f2.xmax;
end


fid = fopen(filename,'w');
fwrite(fid,[86 101 114 115 105 111 110 32 51 46 48 0]','char');  % was f.header.rev  this spells "Version 3.0"
fwrite(fid,0,'long');  % was f.header.nextfile
fwrite(fid,0,'long');  % was f.header.prevfile
fwrite(fid,0,'char');  % was f.header.type
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]' ,'char');  % was f.header.id  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.oper  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.doctor  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.referral  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.hospital  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.patient  This spells "Unspecified"
fwrite(fid,0,'short');  % was f.header.age
fwrite(fid,85,'char');  % was f.header.sex  "U" for unspecified
fwrite(fid,85,'char');  % was f.header.hand "U" for unspecified
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.med  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.category  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.state  This spells "Unspecified"
fwrite(fid,[85 110 115 112 101 99 105 102 105 101 100 0 0 0 0 0 0 0 0 0]','char');  % was f.header.label  This spells "Unspecified"
fwrite(fid,[48 48 47 48 48 47 48 48 0 0]','char');  % was f.header.date  This is a bogus 00/00/00 date
fwrite(fid,[48 48 58 48 48 58 48 48 0 0 0 0]','char');  % was f.header.time  This is a bogus 00:00:00 time.
fwrite(fid,0,'float');  % was f.header.mean_age
fwrite(fid,0,'float');  % was f.header.stdev  ************* toggle???? **************
fwrite(fid,f.nsweeps,'short');  % was f.header.n
fwrite(fid,zeros(38,1),'char');  % was f.header.compfile  (blank)
fwrite(fid,0,'float');  % was f.header.spectwincomp
fwrite(fid,0,'float');  % was f.header.meanaccuracy
fwrite(fid,0,'float');  % was f.header.meanlatency
fwrite(fid,zeros(46,1),'char');  % was f.header.sortfile
fwrite(fid,0,'int');  % was f.header.numevents
fwrite(fid,0,'char');  % was f.header.compoper
fwrite(fid,0,'char');  % was f.header.avgmode
fwrite(fid,0,'char');  % was f.header.review
fwrite(fid,5000,'ushort');  % was f.header.nsweeps  ******** does this matter? ******
fwrite(fid,f.nsweeps,'ushort');  % was f.header.compsweeps 
fwrite(fid,f.nsweeps,'ushort');  % was f.header.acceptcnt 
fwrite(fid,0,'ushort');  % was f.header.rejectcnt **** note; this was 792 on 322-f ******
fwrite(fid,f.pnts,'ushort');  % was f.header.pnts
fwrite(fid,size(f.signal,2),'ushort');  % was f.header.nchannels 
fwrite(fid,50,'ushort');  % was f.header.avgupdate ******** does this matter? ******
fwrite(fid,0,'char');  % was f.header.domain
fwrite(fid,int2str(~isempty(f.variance)),'char');  % was f.header.variance; toggles 1/0
fwrite(fid,f.rate,'ushort');  % was f.header.rate
fwrite(fid,17.1875,'double');  % was f.header.scale   ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.veogcorrect
fwrite(fid,0,'char');  % was f.header.heogcorrect
fwrite(fid,0,'char');  % was f.header.aux1correct
fwrite(fid,0,'char');  % was f.header.aux2correct
fwrite(fid,10,'float');  % was f.header.veogtrig ******* does this matter? *******
fwrite(fid,10,'float');  % was f.header.heogtrig ******* does this matter? *******
fwrite(fid,10,'float');  % was f.header.aux1trig ******* does this matter? *******
fwrite(fid,10,'float');  % was f.header.aux2trig ******* does this matter? *******
fwrite(fid,0,'short');  % was f.header.heogchnl
fwrite(fid,0,'short');  % was f.header.veogchnl
fwrite(fid,0,'short');  % was f.header.aux1chnl
fwrite(fid,0,'short');  % was f.header.aux2chnl
fwrite(fid,0,'char');  % was f.header.veogdir
fwrite(fid,0,'char');  % was f.header.heogdir
fwrite(fid,0,'char');  % was f.header.aux1dir
fwrite(fid,0,'char');  % was f.header.aux2dir
fwrite(fid,0,'short');  % was f.header.veog_n
fwrite(fid,0,'short');  % was f.header.heog_n
fwrite(fid,0,'short');  % was f.header.aux1_n
fwrite(fid,0,'short');  % was f.header.aux2_n
fwrite(fid,20,'short');  % was f.header.veogmaxcnt ******* does this matter? *******
fwrite(fid,20,'short');  % was f.header.heogmaxcnt ******* does this matter? *******
fwrite(fid,20,'short');  % was f.header.aux1maxcnt ******* does this matter? *******
fwrite(fid,20,'short');  % was f.header.aux2maxcnt ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.veogmethod
fwrite(fid,0,'char');  % was f.header.heogmethod
fwrite(fid,0,'char');  % was f.header.aux1method
fwrite(fid,0,'char');  % was f.header.aux2method
fwrite(fid,10,'float');  % was f.header.ampsensitivity ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.lowpass
fwrite(fid,0,'char');  % was f.header.highpass
fwrite(fid,0,'char');  % was f.header.notch
fwrite(fid,0,'char');  % was f.header.autoclipadd
fwrite(fid,1,'char');  % was f.header.baseline ******* does this matter? *******
    % For next several lines, openavg's xmin and max are expressed in msec; must be converted to seconds
fwrite(fid,f.xmin/1000,'float');  % was f.header.offstart
fwrite(fid,f.xmax/1000,'float');  % was f.header.offstop
fwrite(fid,1,'char');  % was f.header.reject ******* does this matter? *******
fwrite(fid,f.xmin/1000,'float');  % was f.header.rejstart
fwrite(fid,f.xmax/1000,'float');  % was f.header.rejstop
fwrite(fid,-100,'float');  % was f.header.rejmin ******* does this matter? *******
fwrite(fid,100,'float');  % was f.header.rejmax ******* does this matter? *******
fwrite(fid,2,'char');  % was f.header.trigtype ******* does this matter? *******
fwrite(fid,0,'float');  % was f.header.trigval
fwrite(fid,16,'char');  % was f.header.trigchnl ******* does this matter? *******
fwrite(fid,0,'short');  % was f.header.trigmask 
fwrite(fid,1,'float');  % was f.header.trigisi ******* does this matter? *******
fwrite(fid,0,'float');  % was f.header.trigmin
fwrite(fid,5,'float');  % was f.header.trigmax ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.trigdir
fwrite(fid,0,'char');  % was f.header.autoscale *** always 0 whether auto or not; how to work this?
fwrite(fid,0,'short');  % was f.header.n2
fwrite(fid,0,'char');  % was f.header.dir  0=neg down; 1=neg up
% **** presumably next 2 lines work in conjuction with f.header.autoscale,
% but I can't get autoscaling to work.  Thus, I'm calculating the maximum
% size of signal, and scaling proportionally.
scaling = max(max(abs(f.signal))).*1.1;
fwrite(fid,-scaling,'float');  % was f.header.dispmin
fwrite(fid,scaling,'float');  % was f.header.dispmax
fwrite(fid,f.xmin/1000,'float');  % was f.header.xmin
fwrite(fid,f.xmax/1000,'float');  % was f.header.xmax
fwrite(fid,f.xmin/1000,'float');  % was f.header.automin
fwrite(fid,f.xmax/1000,'float');  % was f.header.automax
fwrite(fid,0,'float');  % was f.header.zmin
fwrite(fid,0.1,'float');  % was f.header.zmax ******* does this matter? *******
fwrite(fid,0.1,'float');  % was f.header.lowcut ******* does this matter? *******
fwrite(fid,70,'float');  % was  ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.common
fwrite(fid,0,'char');  % was f.header.savemode
fwrite(fid,0,'char');  % was f.header.manmode
fwrite(fid,zeros(10,1),'char');  % was f.header.ref **** provisional; not actually all zeros
fwrite(fid,0,'char');  % was f.header.rectify
fwrite(fid,f.xmin/1000,'float');  % was f.header.displayxmin
fwrite(fid,f.xmax/1000,'float');  % was f.header.displayxmax
fwrite(fid,0,'char');  % was f.header.phase
fwrite(fid,zeros(16,1),'char');  % was f.header.screen  **** provisional; not actually all zeros
fwrite(fid,2,'short');  % was f.header.calmode ******* does this matter? *******
fwrite(fid,0,'short');  % was f.header.calmethod
fwrite(fid,1,'short');  % was f.header.calupdate ******* does this matter? *******
fwrite(fid,0,'short');  % was f.header.calbaseline
fwrite(fid,5,'short');  % was f.header.calsweeps ******* does this matter? *******
fwrite(fid,1,'float');  % was f.header.calattenuator ******* does this matter? *******
fwrite(fid,0.2,'float');  % was f.header.calpulsevolt ******* does this matter? *******
fwrite(fid,0,'float');  % was f.header.calpulsestart
fwrite(fid,0,'float');  % was f.header.calpulsestop
fwrite(fid,20,'float');  % was f.header.calfreq  ******* does this matter? *******
fwrite(fid,zeros(34,1),'char');  % was f.header.taskfile  *** provisional; not actually all zeros
fwrite(fid,zeros(34,1),'char');  % was f.header.seqfile  *** provisional; not actually all zeros
fwrite(fid,0,'char');  % was f.header.spectmethod
fwrite(fid,0,'char');  % was f.header.spectscaling
fwrite(fid,0,'char');  % was f.header.spectwindow
fwrite(fid,0.1,'float');  % was f.header.spectwinlength ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.spectorder
fwrite(fid,0,'char');  % was f.header.notchfilter
fwrite(fid,30,'short');  % was f.header.headgain  ***** do we need to keep track of 30/150?
fwrite(fid,0,'int');  % was f.header.additionalfiles
fwrite(fid,zeros(5,1),'char');  % was f.header.unused
fwrite(fid,0,'short');  % was f.header.fspstopmethod
fwrite(fid,0,'short');  % was f.header.fspstopmode
fwrite(fid,2.5,'float');  % was f.header.fspfvalue ******* does this matter? *******
fwrite(fid,0,'short');  % was f.header.fsppoint
fwrite(fid,200,'short');  % was f.header.fspblocksize ******* does this matter? *******
fwrite(fid,0,'ushort');  % was f.header.fspp1
fwrite(fid,0,'ushort');  % was f.header.fspp2
fwrite(fid,0,'float');  % was f.header.fspalpha
fwrite(fid,0,'float');  % was f.header.fspnoise
fwrite(fid,0,'short');  % was f.header.fspv1
fwrite(fid,zeros(40,1),'char');  % was f.header.montage
fwrite(fid,zeros(40,1),'char');  % was f.header.eventfile
fwrite(fid,1.6667,'float');  % was f.header.fratio ******* does this matter? *******
fwrite(fid,9,'char');  % was f.header.minor_rev ******* does this matter? *******
fwrite(fid,1,'short');  % was f.header.eegupdate ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.compressed
fwrite(fid,0,'float');  % was f.header.xscale
fwrite(fid,0,'float');  % was f.header.yscale
fwrite(fid,42.6762,'float');  % was f.header.xsize ******* does this matter? *******
fwrite(fid,19.1567,'float');  % was f.header.ysize ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.acmode ******* toggle? matters? ********
fwrite(fid,0,'uchar');  % was f.header.commonchnl
fwrite(fid,0,'char');  % was f.header.xtics
fwrite(fid,0,'char');  % was f.header.xrange
fwrite(fid,0,'char');  % was f.header.ytics
fwrite(fid,0,'char');  % was f.header.yrange
fwrite(fid,670.2626,'float');  % was f.header.xscalevalue ******* does this matter? *******
fwrite(fid,1280,'float');  % was f.header.xscaleinterval ******* does this matter? *******
fwrite(fid,50,'float');  % was f.header.yscalevalue ******* does this matter? *******
fwrite(fid,50,'float');  % was f.header.yscaleinterval ******* does this matter? *******
fwrite(fid,86.265,'float');  % was f.header.scaletoolx1 ******* does this matter? *******
fwrite(fid,186.4966,'float');  % was f.header.scaletooly1 ******* does this matter? *******
fwrite(fid,134.7468,'float');  % was f.header.scaletoolx2 ******* does this matter? *******
fwrite(fid,181.7075,'float');  % was f.header.scaletooly2 ******* does this matter? *******
fwrite(fid,715,'short');  % was f.header.port ******* does this matter? *******
fwrite(fid,1865200,'long');  % was f.header.numsamples ******* does this matter? *******
fwrite(fid,0,'char');  % was f.header.filterflag 
fwrite(fid,10,'float');  % was f.header.lowcutoff ******* does this matter? *******
fwrite(fid,2,'short');  % was f.header.lowpoles ******* does this matter? *******
fwrite(fid,125,'float');  % was f.header.highcutoff ******* does this matter? *******
fwrite(fid,2,'short');  % was f.header.highpoles ******* does this matter? *******
fwrite(fid,3,'char');  % was f.header.filtertype ******* does this matter? *******
fwrite(fid,1,'char');  % was f.header.filterdomain
fwrite(fid,0,'char');  % was f.header.snrflag
fwrite(fid,0,'char');  % was f.header.coherenceflag
fwrite(fid,3,'char');  % was f.header.continuoustype ******* does this matter? *******
fwrite(fid,37305650,'long');  % was f.header.eventtablepos ******* does this matter? *******
fwrite(fid,0,'float');  % was f.header.continuousseconds
fwrite(fid,80,'long');  % was f.header.channeloffset ******* does this matter? *******
fwrite(fid,1,'char');  % was f.header.autocorrectflag ******* does this matter? *******
fwrite(fid,50,'uchar');  % was f.header.dcthreshold ******* does this matter? *******

% munch channel names (convert from #x10 char array to cell array)
for munch = 1:size(f.signal,2) 
    chan_names{munch} = double(f.chan_names(munch,:));
end

ypos = 18; % seed for y position
xpos = 5;  % seed for x position
for n = 1:size(f.signal,2)
    fwrite(fid,chan_names{n}','char');  % was f.electloc(n).lab this is electrode name
    fwrite(fid,0,'char');  % was f.electloc(n).reference
    fwrite(fid,0,'char');  % was f.electloc(n).skip
    fwrite(fid,1,'char');  % was f.electloc(n).reject
    fwrite(fid,1,'char');  % was f.electloc(n).display
    fwrite(fid,0,'char');  % was f.electloc(n).bad
    fwrite(fid,f.nsweeps,'ushort');  % was f.electloc(n).n
    fwrite(fid,0,'char');  % was f.electloc(n).avg_reference
    fwrite(fid,0,'char');  % was f.electloc(n).clipadd
    fwrite(fid,xpos,'float');  % was f.electloc(n).x_coord *** places them horizontally
    fwrite(fid,ypos,'float');  % was f.electloc(n).y_coord *** places them vertically
    fwrite(fid,0,'float');  % was f.electloc(n).veog_wt
    fwrite(fid,0,'float');  % was f.electloc(n).veog_std
    fwrite(fid,0,'float');  % was f.electloc(n).snr
    fwrite(fid,0,'float');  % was f.electloc(n).heog_wt
    fwrite(fid,0,'float');  % was f.electloc(n).heog_std
    fwrite(fid,0,'short');  % was f.electloc(n).baseline
    fwrite(fid,1,'char');  % was f.electloc(n).filtered  *** does this matter? ****
    fwrite(fid,0,'char');  % was f.electloc(n).fsp
    fwrite(fid,0,'float');  % was f.electloc(n).aux1_wt
    fwrite(fid,0,'float');  % was f.electloc(n).aux1_std
    fwrite(fid,17.1875,'float');  % was f.electloc(n).sensitivity *** does this matter? ****
    fwrite(fid,5,'char');  % was f.electloc(n).gain *** does this matter? ****
    fwrite(fid,0,'char');  % was f.electloc(n).hipass
    fwrite(fid,4,'char');  % was f.electloc(n).lopass *** does this matter? ****
    fwrite(fid,0,'uchar');  % was f.electloc(n).page
    fwrite(fid,1,'uchar');  % was f.electloc(n).size
    fwrite(fid,0,'uchar');  % was f.electloc(n).impedance
    fwrite(fid,n-1,'uchar');  % was f.electloc(n).physicalchnl *** does this matter? ****
    fwrite(fid,0,'char');  % was f.electloc(n).rectify
    fwrite(fid,1,'float');  % was f.electloc(n).calib
    xpos = xpos + 100;  % increments x position
    if n/6 == fix(n/6)  % increments y position (resetting x) when row is full (6-channels).
        ypos = ypos + 45;
        xpos = 5;
    end
end

for i = 1:size(f.signal,2)
   fwrite(fid,chan_names{i}(1:5)','char');  % was f.data(i).header (note 1st 5 chars only)
   fwrite(fid,f.signal(:,i).*f.nsweeps,'float');  % was f.data(i).samples
end

% if f.variance is null, we need to create one...even tho the variance flag
% above is set to 0.
if isempty(f.variance) == 1
     f.variance = zeros(f.pnts,size(f.signal,2));
end
    
for j = 1:size(f.signal,2) 
    fwrite(fid,f.variance(:,j),'float');  % was f.variance(j).samples
end

fwrite(fid,zeros(0,1),'char');  % was f.tag  f.tag is an empty matrix; size 0x1


frewind(fid);

fclose(fid);
