function bt_txt2avg(txtfile, fs, xmin, xmax);
%
%  bt_txt2avg reads in a selected text file and creates an *.avg file that stores all of the file information including 
%  channel names, and sampling rate. 
%
% Usage: bt_txt2avg('Filename.txt', 20000, -15, 60, 'CZ')
%
%
% Dependancies: writeavg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


txt = load(txtfile);  %read text file

% create
f.pnts = length(txt);
f.rate = fs;
f.xmin = xmin;
f.xmax = xmax;
f.chan_names = 'bt_toolbox';
f.variance = [];
f.nsweeps = 1;
f.signal = txt;

f.signal = f.signal-mean(f.signal(1:ms2row(f,0))); % pre-stimulus baseline

[path name ext]=fileparts(txtfile);

if isempty(path); % if the file is in the current directory.
     writeavg(f, [name, '.avg'])
else
    writeavg(f, fullfile(path, strcat(name, '.avg')));
end



