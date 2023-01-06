function eeg_write(f,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified version of eeg_write_scan41
%
% Useage:   used in conjunction with eeg_load. eeg_load opens Neuroscan avg
%           files in the form of a structural array, f.
%
%            f is a structure containing:
%
%               f.header
%               f.electloc
%               f.data
%               f.variance
%               f.tag
%               f.chan_names
%
%           eeg_write(f, filename)
%
% where:    filename is a complete 'path\fileprefix.ext' string;
%
%
% original Author:  unknown
% Created:  unknown
% 
%Modification history:
% (1) 7/23/2004 Erika Skoe eeskoe@northwestern.edu
%   	modified to write avg files that can be read by scan4.1, scan4.2 and scan4.3.
%	eeg_load scales the signal to microvolts (from amplifier units), by dividing by the number of sweeps in the
%	average which is either f.header.compsweeps (total sweeps) OR f.header.acceptcnt (accepted sweeps).
% 	In eeg_load, if the value of h.compsweeps equals 1, then the signal (i.e.
% 	f.data.samples) is already scaled into uV. If the value is something
% 	other than 1, the signal must be scaled into uV by dividing the signal by
% 	the the value of h.acceptcnt.  (see 12/2004 modification)
%
% 	Before the signal can be written to file is must be scaled back to to
%	 amplifier units.
%
% (2) 11/2004  Erika Skoe eeskoe@northwestern.edu
% 	If header.minor_rev is larger than 15, than the data is stored in 32-bit
% 	integers, otherwise it is 16-bit integers.   When eeg_write generated avg
% 	files will minor_rev > 13, the resulting avg could not be opened in
% 	N-SCAN. As a temporary solution, minor_rev is changed to 13, when the minor_rev of the template avg file
% 	is greater than 16.
% (3) 12/2004  Erika Skoe eeskoe@northwestern.edu
% 	Modification to rescale grand-average files back to amplifier units.
%	In eeg_load, if the value of h.compsweeps equals 1, then the signal (i.e.
% 	f.data.samples) is already scaled into uV. If the value is something
%	 other than 1, then first check whether h.nsweeps = h.compsweeps. If
% 	yes, then signal must be scaled by the value of h.compsweeps. If no,
% 	the signal must be scaled into uV by dividing the signal by
%	 the the value of h.acceptcnt.   For grandaverage files, h.nsweeps =
%	 h.compsweeps, and this value is equivalent to the number of files that
%	 were used to make the grandaverage.
% (4) 1/10/05 Erika Skoe eeskoe@northwestern.edu
% 	Signal is scaled using the baseline and calibration values of each electrode.
%       signal = ((signal + baseline)/calibration factor)*sweeps
%       For all files surveyed thus far, baseline = 0, and calib = 1;
% (5) 8/31/07  Jade Wang jadewang@gmail.com
% 	Marker information is written into f.tag.  Existing marker information is erased.
%	 F. marker must take on the following format
% 	f.marker{x} = name: 	'New                 '
%     				lat_pts: 623
%     				latency: 0.0212
%                      		channel: 'CZ-A2     '
% 	 			ch_num: 1
%(6) March 2008:  Erika Skoe eeskoe@northwestern.edu
%	(a) To be read correctly by Neuroscan, f.xmax must be set to xmax (display value) +1;
%        To accomondate this, f.xmax is directly calculated from f.xmin, pnts and sample rate (first line of code). 
%        calculated from the sample rate, f.xmin and the number of points.  
%        If the calculated value does not equal f.xmax+1pt, a warning message is displayed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax = f.header.xmin+((f.header.pnts)*(1/f.header.rate));   %% added by eeskoe (3/2008)

if xmax ~=f.header.xmax+((1/f.header.rate))
    display(['warning: xmax adjusted to ',  num2str(xmax-((1/f.header.rate)))]);
else
    f.header.xmax = xmax;
end


% for i = 1:f.header.nchannels;
%     if strmatch('GFP      ', f.chan_names(i,:))
%         f.electloc(i).lab = double('GFP Matlab');
%     end           
%     if strmatch('REF      ', f.chan_names(i,:))
%         f.electloc(i).lab = double('REF Matlab');
%     end
% end
 

fid = fopen(filename,'w');

fwrite(fid,f.header.rev,'char');
fwrite(fid,f.header.nextfile,'long');
fwrite(fid,f.header.prevfile,'long');
fwrite(fid,f.header.type,'char');
fwrite(fid,f.header.id,'char');
fwrite(fid,f.header.oper,'char');
fwrite(fid,f.header.doctor,'char');
fwrite(fid,f.header.referral,'char');
fwrite(fid,f.header.hospital,'char');
fwrite(fid,f.header.patient,'char');
fwrite(fid,f.header.age,'short');
fwrite(fid,f.header.sex,'char');
fwrite(fid,f.header.hand,'char');
fwrite(fid,f.header.med,'char');
fwrite(fid,f.header.category,'char');
fwrite(fid,f.header.state,'char');
fwrite(fid,f.header.label,'char');
fwrite(fid,f.header.date,'char');
fwrite(fid,f.header.time,'char');
fwrite(fid,f.header.mean_age,'float');
fwrite(fid,f.header.stdev,'float');
fwrite(fid,f.header.n,'short');
fwrite(fid,f.header.compfile,'char');
fwrite(fid,f.header.spectwincomp,'float');
fwrite(fid,f.header.meanaccuracy,'float');
fwrite(fid,f.header.meanlatency,'float');
fwrite(fid,f.header.sortfile,'char');
fwrite(fid,f.header.numevents,'int');
fwrite(fid,f.header.compoper,'char');
fwrite(fid,f.header.avgmode,'char');
fwrite(fid,f.header.review,'char');
fwrite(fid,f.header.nsweeps,'ushort');
fwrite(fid,f.header.compsweeps,'ushort');
fwrite(fid,f.header.acceptcnt,'ushort');
fwrite(fid,f.header.rejectcnt,'ushort');
fwrite(fid,f.header.pnts,'ushort');
fwrite(fid,f.header.nchannels,'ushort');
fwrite(fid,f.header.avgupdate,'ushort');
fwrite(fid,f.header.domain,'char');
fwrite(fid,f.header.variance,'char');
fwrite(fid,f.header.rate,'ushort');
fwrite(fid,f.header.scale,'double');
fwrite(fid,f.header.veogcorrect,'char');
fwrite(fid,f.header.heogcorrect,'char');
fwrite(fid,f.header.aux1correct,'char');
fwrite(fid,f.header.aux2correct,'char');
fwrite(fid,f.header.veogtrig,'float');
fwrite(fid,f.header.heogtrig,'float');
fwrite(fid,f.header.aux1trig,'float');
fwrite(fid,f.header.aux2trig,'float');
fwrite(fid,f.header.heogchnl,'short');
fwrite(fid,f.header.veogchnl,'short');
fwrite(fid,f.header.aux1chnl,'short');
fwrite(fid,f.header.aux2chnl,'short');
fwrite(fid,f.header.veogdir,'char');
fwrite(fid,f.header.heogdir,'char');
fwrite(fid,f.header.aux1dir,'char');
fwrite(fid,f.header.aux2dir,'char');
fwrite(fid,f.header.veog_n,'short');
fwrite(fid,f.header.heog_n,'short');
fwrite(fid,f.header.aux1_n,'short');
fwrite(fid,f.header.aux2_n,'short');
fwrite(fid,f.header.veogmaxcnt,'short');
fwrite(fid,f.header.heogmaxcnt,'short');
fwrite(fid,f.header.aux1maxcnt,'short');
fwrite(fid,f.header.aux2maxcnt,'short');
fwrite(fid,f.header.veogmethod,'char');
fwrite(fid,f.header.heogmethod,'char');
fwrite(fid,f.header.aux1method,'char');
fwrite(fid,f.header.aux2method,'char');
fwrite(fid,f.header.ampsensitivity,'float');
fwrite(fid,f.header.lowpass,'char');
fwrite(fid,f.header.highpass,'char');
fwrite(fid,f.header.notch,'char');
fwrite(fid,f.header.autoclipadd,'char');
fwrite(fid,f.header.baseline,'char');
fwrite(fid,f.header.offstart,'float');
fwrite(fid,f.header.offstop,'float');
fwrite(fid,f.header.reject,'char');
fwrite(fid,f.header.rejstart,'float');
fwrite(fid,f.header.rejstop,'float');
fwrite(fid,f.header.rejmin,'float');
fwrite(fid,f.header.rejmax,'float');
fwrite(fid,f.header.trigtype,'char');
fwrite(fid,f.header.trigval,'float');
fwrite(fid,f.header.trigchnl,'char');
fwrite(fid,f.header.trigmask,'short');
fwrite(fid,f.header.trigisi,'float');
fwrite(fid,f.header.trigmin,'float');
fwrite(fid,f.header.trigmax,'float');
fwrite(fid,f.header.trigdir,'char');
fwrite(fid,f.header.autoscale,'char');
fwrite(fid,f.header.n2,'short');
fwrite(fid,f.header.dir,'char');
fwrite(fid,f.header.dispmin,'float');
fwrite(fid,f.header.dispmax,'float');
fwrite(fid,f.header.xmin,'float');
fwrite(fid,f.header.xmax,'float');
fwrite(fid,f.header.automin,'float');
fwrite(fid,f.header.automax,'float');
fwrite(fid,f.header.zmin,'float');
fwrite(fid,f.header.zmax,'float');
fwrite(fid,f.header.lowcut,'float');
fwrite(fid,f.header.highcut,'float');
fwrite(fid,f.header.common,'char');
fwrite(fid,f.header.savemode,'char');
fwrite(fid,f.header.manmode,'char');
fwrite(fid,f.header.ref,'char');
fwrite(fid,f.header.rectify,'char');
fwrite(fid,f.header.displayxmin,'float');
fwrite(fid,f.header.displayxmax,'float');
fwrite(fid,f.header.phase,'char');
fwrite(fid,f.header.screen,'char');
fwrite(fid,f.header.calmode,'short');
fwrite(fid,f.header.calmethod,'short');
fwrite(fid,f.header.calupdate,'short');
fwrite(fid,f.header.calbaseline,'short');
fwrite(fid,f.header.calsweeps,'short');
fwrite(fid,f.header.calattenuator,'float');
fwrite(fid,f.header.calpulsevolt,'float');
fwrite(fid,f.header.calpulsestart,'float');
fwrite(fid,f.header.calpulsestop,'float');
fwrite(fid,f.header.calfreq,'float');
fwrite(fid,f.header.taskfile,'char');
fwrite(fid,f.header.seqfile,'char');
fwrite(fid,f.header.spectmethod,'char');
fwrite(fid,f.header.spectscaling,'char');
fwrite(fid,f.header.spectwindow,'char');
fwrite(fid,f.header.spectwinlength,'float');
fwrite(fid,f.header.spectorder,'char');
fwrite(fid,f.header.notchfilter,'char');
fwrite(fid,f.header.headgain,'short');
fwrite(fid,f.header.additionalfiles,'int');
fwrite(fid,f.header.unused,'char');
fwrite(fid,f.header.fspstopmethod,'short');
fwrite(fid,f.header.fspstopmode,'short');
fwrite(fid,f.header.fspfvalue,'float');
fwrite(fid,f.header.fsppoint,'short');
fwrite(fid,f.header.fspblocksize,'short');
fwrite(fid,f.header.fspp1,'ushort');
fwrite(fid,f.header.fspp2,'ushort');
fwrite(fid,f.header.fspalpha,'float');
fwrite(fid,f.header.fspnoise,'float');
fwrite(fid,f.header.fspv1,'short');
fwrite(fid,f.header.montage,'char');
fwrite(fid,f.header.eventfile,'char');
fwrite(fid,f.header.fratio,'float');


if f.header.minor_rev == 16
    f.header.minor_rev = 13;
end

fwrite(fid,f.header.minor_rev,'char');
fwrite(fid,f.header.eegupdate,'short');
fwrite(fid,f.header.compressed,'char');
fwrite(fid,f.header.xscale,'float');
fwrite(fid,f.header.yscale,'float');
fwrite(fid,f.header.xsize,'float');
fwrite(fid,f.header.ysize,'float');
fwrite(fid,f.header.acmode,'char');
fwrite(fid,f.header.commonchnl,'uchar');
fwrite(fid,f.header.xtics,'char');
fwrite(fid,f.header.xrange,'char');
fwrite(fid,f.header.ytics,'char');
fwrite(fid,f.header.yrange,'char');
fwrite(fid,f.header.xscalevalue,'float');
fwrite(fid,f.header.xscaleinterval,'float');
fwrite(fid,f.header.yscalevalue,'float');
fwrite(fid,f.header.yscaleinterval,'float');
fwrite(fid,f.header.scaletoolx1,'float');
fwrite(fid,f.header.scaletooly1,'float');
fwrite(fid,f.header.scaletoolx2,'float');
fwrite(fid,f.header.scaletooly2,'float');
fwrite(fid,f.header.port,'short');
fwrite(fid,f.header.numsamples,'long');
fwrite(fid,f.header.filterflag,'char');
fwrite(fid,f.header.lowcutoff,'float');
fwrite(fid,f.header.lowpoles,'short');
fwrite(fid,f.header.highcutoff,'float');
fwrite(fid,f.header.highpoles,'short');
fwrite(fid,f.header.filtertype,'char');
fwrite(fid,f.header.filterdomain,'char');
fwrite(fid,f.header.snrflag,'char');
fwrite(fid,f.header.coherenceflag,'char');
fwrite(fid,f.header.continuoustype,'char');
fwrite(fid,f.header.eventtablepos,'long');
fwrite(fid,f.header.continuousseconds,'float');
fwrite(fid,f.header.channeloffset,'long');
fwrite(fid,f.header.autocorrectflag,'char');
fwrite(fid,f.header.dcthreshold,'uchar');

for n = 1:f.header.nchannels
    fwrite(fid,f.electloc(n).lab,'char');
    fwrite(fid,f.electloc(n).reference,'char');
    fwrite(fid,f.electloc(n).skip,'char');
    fwrite(fid,f.electloc(n).reject,'char');
    fwrite(fid,f.electloc(n).display,'char');
    fwrite(fid,f.electloc(n).bad,'char');
    fwrite(fid,f.electloc(n).n,'ushort');
    fwrite(fid,f.electloc(n).avg_reference,'char');
    fwrite(fid,f.electloc(n).clipadd,'char');
    fwrite(fid,f.electloc(n).x_coord,'float');
    fwrite(fid,f.electloc(n).y_coord,'float');
    fwrite(fid,f.electloc(n).veog_wt,'float');
    fwrite(fid,f.electloc(n).veog_std,'float');
    fwrite(fid,f.electloc(n).snr,'float');
    fwrite(fid,f.electloc(n).heog_wt,'float');
    fwrite(fid,f.electloc(n).heog_std,'float');
    fwrite(fid,f.electloc(n).baseline,'short');
    fwrite(fid,f.electloc(n).filtered,'char');
    fwrite(fid,f.electloc(n).fsp,'char');
    fwrite(fid,f.electloc(n).aux1_wt,'float');
    fwrite(fid,f.electloc(n).aux1_std,'float');
    fwrite(fid,f.electloc(n).senstivity,'float');
    fwrite(fid,f.electloc(n).gain,'char');
    fwrite(fid,f.electloc(n).hipass,'char');
    fwrite(fid,f.electloc(n).lopass,'char');
    fwrite(fid,f.electloc(n).page,'uchar');
    fwrite(fid,f.electloc(n).size,'uchar');
    fwrite(fid,f.electloc(n).impedance,'uchar');
    fwrite(fid,f.electloc(n).physicalchnl,'uchar');
    fwrite(fid,f.electloc(n).rectify,'char');
    fwrite(fid,f.electloc(n).calib,'float');
end

% eeskoe modification:
%eeg_load scales the signal to microvolts (from amplifier units), by dividing by the number of sweeps in the
%average which is either f.header.compsweeps (total sweeps) OR f.header.acceptcnt (accepted sweeps).
% In eeg_load, if the value of h.compsweeps equals 1, then the signal (i.e.
% f.data.samples) is already scaled into uV. If the value is something
% other than 1, the signal must be scaled into uV by dividing the signal by
% the the value of h.acceptcnt.

% Before the signal can be written to file is must be scaled back to to
% amplifier units.

for i = 1:f.header.nchannels;
    fwrite(fid,f.data(i).header,'char');

    if f.header.compsweeps ~= 1
        if f.header.compsweeps == f.header.nsweeps
            if findstr('GFP', f.chan_names(i,:))
                f.data(i).samples =  ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib);  
            elseif findstr('REF', f.chan_names(i,:))
                f.data(i).samples =  ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib); 
            else  
                f.data(i).samples = ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib).*f.header.nsweeps;  % added 12/2004
            end
        else
            
            if findstr('GFP', f.chan_names(i,:))
                f.data(i).samples =  ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib);  
            elseif findstr('REF', f.chan_names(i,:))
                f.data(i).samples =  ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib); 
            else  
                f.data(i).samples = ((f.data(i).samples+f.electloc(i).baseline)./f.electloc(i).calib).*f.header.acceptcnt;  % added 12/2004
            end
            
            
            
            
            
        end
    end
    fwrite(fid,f.data(i).samples,'float');

end





% % correct the signal
% for i = 1:f.header.nchannels;
%    
%     if f.header.compsweeps ~= 1
%         
%         if f.header.compsweeps == f.header.nsweeps
% 
%             if findstr('GFP', f.chan_names(i,:))
% 
%                 f.data(i).samples = ((f.data(i).samples)./f.header.nsweeps);
%             elseif findstr('REF', f.chan_names(i,:))
%                 f.data(i).samples = ((f.data(i).samples)./f.header.nsweeps);
%                 %
%             end
%         else
%             
%             %
%             if findstr('GFP', f.chan_names(i,:))
%                           
%                 f.data(i).samples = ((f.data(i).samples)./f.header.acceptcnt);
%                 figure; plot(f.data(i).samples);
%                 f.variance(i).samples = ((f.variance(i).samples)/f.header.acceptcnt);
%             elseif findstr('REF', f.chan_names(i,:))
%                 f.data(i).samples = ((f.data(i).samples)./f.header.acceptcnt);
%                 f.variance(i).samples = ((f.variance(i).samples)./f.header.acceptcnt);
%                 
%             end
%         end
%     end
% end
% 
%     figure; plot(f.data(34).samples);    






% jqw modification:
% f.marker contains marker information extracted from f.tag
% any modifications to f.marker needs to be reflected in a new f.tag before
% being written to file.

% basic parameters

markerprefix  = [80 195 0 0 5 1 0 0 152 0 0 0 152 0 0 0 255 0 0 0];
tagpostfix    = [112 17 1 0 0 0 0 0 4 0 0 0 4 0 0 0 1 0 0 0];
markerpostfix = [18 0 1 0 0 0 0 0 0 0];
seqlen        = 20;
pad1          = zeros(1,30);
pad2          = [75 0 0 0 0 0 0 0 0 0 1 0 0 0]; % meaningless sequence
pad3          = zeros(1,10);
pad4          = zeros(1,62);
factor        = 256;
markerlen     = 168;

% part 1: does marker exist?  if so, how many?

if isstr(f.marker)             % f.marker is str 'No Markers' or cell array
    n_markers = 0;
else
    n_markers = length(f.marker);
end

% part 2: extract non-marker f.tag info (erase existing markers)

len          = length(f.tag);
markerstart  = strfind(f.tag',markerprefix);
n_marker_old = length(markerstart);

if n_marker_old >=1
    tagprefix = f.tag(1:markerstart(1)-1)';
else
    tagprefix = f.tag(1:end-seqlen)';
end

% part 3: if n_markers > 0, write new markers in

allmarkers = [];

if n_markers > 0
    for x = 1:n_markers
        markername = f.marker{x}.name + 0;   % force it to be numerical
        markerchan = f.marker{x}.ch_num - 1; % channels count from 0
        markerpts  = f.marker{x}.lat_pts -1; % latency in pts trumps all.
        %         To be read properly by eeg_load, need to subtract 1 (eeskoe)
        lat1 = mod(markerpts,factor);
        lat2 = floor(markerpts/factor);
        markerlat  = [lat1 lat2];
        mymarker   = cat(2,markerprefix, markername, pad1, pad2, markerlat, pad3, markerchan, pad4, markerpostfix);
        %         disp(mymarker');
        allmarkers = cat(2,allmarkers,mymarker);
    end
end

% new f.tag
f.tag = cat(2,tagprefix, allmarkers, tagpostfix)';

fwrite(fid,f.tag,'uchar');

frewind(fid);

fclose(fid);
