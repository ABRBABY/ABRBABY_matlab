function varargout = bt_gui(varargin)
    % The brainstem toolbox  is free software: you can redistribute it and/or modify
    %~ it under the terms of the GNU General Public License as published by
    %~ the Free Software Foundation.
%
    %~ The brainstem toolbox is distributed in the hope that it will be useful,
    %~ but WITHOUT ANY WARRANTY; without even the implied warranty of
    %~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    %~ See the GNU General Public License for more details <http://www.gnu.org/licenses/>/
%
%Copyright 2007,2008 Trent Nicol and Erika Skoe
% Dependencies:
% brainstem toolbox m-files: bt_fft2, bt_fftsc, 
%          bt_peaks2, bt_qncorr, bt_rms, bt_xlswrite
% other m-files: closestrc, ms2row, nancorrcoef, openavg, xcorrelation (at
%          least)
% fig-files: bt_gui, bt_gui_out
%
% 2007-2008, Trent Nicol and Erika Skoe.  Based on "Brainstem Toolbox" by Erika Skoe and
% Trent Nicol.
%
%
%
%
%% Begin auto-generated code
addpath(fullfile(cd, 'programFiles'))
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bt_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @bt_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

               
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before bt_bui is made visible.
function bt_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bt_bui (see VARARGIN)

% Choose default command line output for bt_bui
% guifig = openfig([cd, '\programFiles\bt_gui.fig'],'reuse');   % the same window is used each time within 1 matlab session.   
%                                          % this allows path information to be retained from subject to subject
% set(guifig, 'visible', 'on'); 

set(gcf,'CloseRequestFcn','closeGUI')   %calls 'closeGUI.m' when window is closed.
                                            % when the GUI window is closed,
                                            % the user is returned to the
                                            % start directory.
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = bt_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------end auto-generated code-------------------

%% set defaults ------------------------------------------------------
temp = get(handles.MarkerFileYN, 'Value');
if temp == 0
  set(handles.MarkerFileYN, 'Value', 0);
  set(handles.MarkerFile, 'Enable', 'off');
else
  set(handles.MarkerFile, 'Enable', 'on');
end


temp = get(handles.MarkerFileYN, 'Enable');
if strcmp(temp, 'off')
    set(handles.MarkerFile, 'Enable', 'off');
end

temp = get(handles.NoiseFileYN, 'Value');
if temp == 0
    set(handles.NoiseFileYN, 'Value', 0);
    set(handles.NoiseFile, 'Enable', 'off');
else
     set(handles.NoiseFile, 'Enable', 'on');
end

temp = get(handles.NoiseFile, 'Enable');
if strcmp(temp, 'off')
    set(handles.NoiseFile, 'Enable', 'off');
    set(handles.PlotNoiseFileYN,'Enable', 'off');
end

temp = get(handles.PlotNoiseFileYN, 'Value');
if temp == 0;
    set(handles.PlotNoiseFileYN, 'Value', 0);
end

temp = get(handles.StimFileYN, 'Value');
if temp == 0;
    set(handles.StimFileYN, 'Value', 0);
    set(handles.StimFile, 'Enable', 'off');
else
    set(handles.StimFile, 'Enable', 'on');
end

temp = get(handles.StimFile, 'Enable');
if strcmp(temp, 'off')
    set(handles.StimFile, 'Enable', 'off');
end

temp = get(handles.XL, 'String');
if isempty(char(temp));
     
 set(handles.XL, 'Value', 1);
end


temp = get(handles.StartTime, 'String');

if isempty(char(temp));
    set(handles.StartTime, 'String', 10);
end

temp = get(handles.StopTime, 'String');

if isempty(char(temp));
    set(handles.StopTime, 'String', 50);
end

temp = get(handles.FreqStartTime, 'String');
if isempty(char(temp));
   set(handles.FreqStartTime, 'String', 10);
end

temp = get(handles.FreqStopTime, 'String');
if isempty(char(temp));
   set(handles.FreqStopTime, 'String', 50);
end


temp = get(handles.RRStartTime, 'String');
if isempty(char(temp));
   set(handles.RRStartTime, 'String', 10);
end
temp = get(handles.RRStopTime, 'String');
if isempty(char(temp));
   set(handles.RRStopTime, 'String', 50);
end


temp = get(handles.SRStartTime, 'String');
if isempty(char(temp));
   set(handles.SRStartTime, 'String', 0);
end
temp = get(handles.SRStopTime, 'String');
if isempty(char(temp));
   set(handles.SRStopTime, 'String', 40);
end

temp = get(handles.Freq1Start, 'String');
if isempty(char(temp));
   set(handles.Freq1Start, 'String', 80);
end


temp = get(handles.Freq1Stop, 'String');
if isempty(char(temp));
    set(handles.Freq1Stop, 'String', 120);
end

temp = get(handles.Freq2Start, 'String');
if isempty(char(temp));
   set(handles.Freq2Start, 'String', 200);
end


temp = get(handles.Freq2Stop, 'String');
if isempty(char(temp));
    set(handles.Freq2Stop, 'String', 500);
end


temp = get(handles.Freq3Start, 'String');
if isempty(char(temp));
   set(handles.Freq3Start, 'String', 700);
end


temp = get(handles.Freq3Stop, 'String');
if isempty(char(temp));
    set(handles.Freq3Stop, 'String', 1200);
end


temp = get(handles.FFTScaleYN, 'Value');
if temp == 0;
   set(handles.FFTScaleYN, 'Value', 0);  % scale spectrum to peak uV?
end


% Defaults for marked peaks: names and polarities
temp = get(handles.Mark1, 'String');
if isempty(char(temp));
   set(handles.Mark1, 'String', 1);  
end

temp = get(handles.Pos1, 'Value');
if temp == 0;
   set(handles.Pos1, 'Value', 0);  
end

temp = get(handles.Mark2, 'String');
if isempty(char(temp));
   set(handles.Mark2,  'String', 2);  
end

temp = get(handles.Pos2, 'Value');
if temp == 0;
   set(handles.Pos2, 'Value', 0);  
end

temp = get(handles.Mark3, 'String');
if isempty(char(temp));
   set(handles.Mark3,  'String', 3); 
end

temp = get(handles.Pos3, 'Value');
if temp == 0;
   set(handles.Pos3, 'Value', 0);  
end


temp = get(handles.Mark4, 'String');
if isempty(char(temp));
   set(handles.Mark4,  'String', 4);  
end

temp = get(handles.Pos4, 'Value');
if temp == 0;
   set(handles.Pos4, 'Value', 0); 
end


temp = get(handles.Mark5, 'String');
if isempty(char(temp));
   set(handles.Mark5,  'String', 5);  
end

temp = get(handles.Pos5, 'Value');
if temp == 0;
   set(handles.Pos5, 'Value', 0);  
end

temp = get(handles.Mark6, 'String');
if isempty(char(temp));
   set(handles.Mark6,  'String', 6); 
end

temp = get(handles.Pos6, 'Value');
if temp == 0;
   set(handles.Pos6, 'Value', 0);  
end



temp = get(handles.Mark7, 'String');
if isempty(char(temp));
   set(handles.Mark7,  'String', 7);  % scale spectrum to peak uV?
end

temp = get(handles.Pos7, 'Value');
if temp == 0;
   set(handles.Pos7, 'Value', 0);  % scale spectrum to peak uV?
end


temp = get(handles.Mark8, 'String');
if isempty(char(temp));
   set(handles.Mark8,  'String',8);  % scale spectrum to peak uV?
end

temp = get(handles.Pos8, 'Value');
if temp == 0;
   set(handles.Pos8, 'Value', 0);  % scale spectrum to peak uV?
end


temp = get(handles.Mark9, 'String');
if isempty(char(temp));
   set(handles.Mark9,  'String', 9);  % scale spectrum to peak uV?
end

temp = get(handles.Pos9, 'Value');
if temp == 0;
   set(handles.Pos9, 'Value', 0);  % scale spectrum to peak uV?
end

temp = get(handles.Mark10, 'String');
if isempty(char(temp));
   set(handles.Mark10,  'String', 10);  % scale spectrum to peak uV?
end

temp = get(handles.Pos10, 'Value');
if temp == 0;
   set(handles.Pos10, 'Value', 0);  % scale spectrum to peak uV?
end

%more marker defaults to avoid complicated conditional later. Do not alter

temp = get(handles.Lat1, 'String');
if isempty(char(temp));
  set(handles.Lat1, 'String', '0');
end

temp = get(handles.Lat2, 'String');
if isempty(char(temp));
  set(handles.Lat2, 'String', '0');
end

temp = get(handles.Lat3, 'String');
if isempty(char(temp));
  set(handles.Lat3, 'String', '0');
end

temp = get(handles.Lat4, 'String');
if isempty(char(temp));
  set(handles.Lat4, 'String', '0');
end

temp = get(handles.Lat5, 'String');
if isempty(char(temp));
  set(handles.Lat5, 'String', '0');
end

temp = get(handles.Lat6, 'String');
if isempty(char(temp));
  set(handles.Lat6, 'String', '0');
end

temp = get(handles.Lat7, 'String');
if isempty(char(temp));
  set(handles.Lat7, 'String', '0');
end

temp = get(handles.Lat8, 'String');
if isempty(char(temp));
  set(handles.Lat8, 'String', '0');
end

temp = get(handles.Lat9, 'String');
if isempty(char(temp));
  set(handles.Lat9, 'String', '0');
end


temp = get(handles.Lat10, 'String');
if isempty(char(temp));
  set(handles.Lat10, 'String', '0');
end





%% Some conditional enabling/disabling of buttons ------------------------
function MarkerFileYN_Callback(hObject, eventdata, handles)
MarkerFileYN = get(handles.MarkerFileYN, 'Value');
if MarkerFileYN == 1;
    set(handles.MarkerFile, 'Enable', 'on');
else
    set(handles.MarkerFile, 'Enable', 'off');
end

function NoiseFileYN_Callback(hObject, eventdata, handles)
NoiseFileYN = get(handles.NoiseFileYN, 'Value');
if NoiseFileYN == 1;
    set(handles.NoiseFile, 'Enable', 'on');
    set(handles.PlotNoiseFileYN, 'Enable', 'on');
else
    set(handles.NoiseFile, 'Enable', 'off');
    set(handles.PlotNoiseFileYN, 'Enable', 'off');
end

function StimFileYN_Callback(hObject, eventdata, handles)
StimFileYN = get(handles.StimFileYN, 'Value');
if StimFileYN == 1;
    set(handles.StimFile, 'Enable', 'on');
else
    set(handles.StimFile, 'Enable', 'off');
end

%% get the browse-button info. -------------------------------------------
function QuietFile_Callback(hObject, eventdata, handles)
    global QuietFile QuietFilePath
    [QuietFile QuietFilePath] = uigetfile('*.avg','Select brainstem response file');
    set(handles.QuietFileText, 'String', [ QuietFile]); % this displays

    % open avg file and determine how many channels there are. If more than 1, prompt the user to select the file to process.
    global chan
    if isempty(chan);
        f = openavg([QuietFilePath QuietFile]);
        if size(f.chan_names,1) >1
            [chan v]=  listdlg('PromptString', 'Select Channel to Process' ,'ListString', f.chan_names, 'SelectionMode', 'single') ;
        else
            chan = 1;

        end
        global chan
    end          % the file selected to the right of the browse button.

function MarkerFile_Callback(hObject, eventdata, handles)
    global MarkerFile MarkerFilePath
    [MarkerFile MarkerFilePath] = uigetfile('*.txt;*.xls','Select marker file');
    set(handles.MarkerFileText, 'String', [ MarkerFile]); % this displays 
               % the file selected to the right of the browse button.
               
function NoiseFile_Callback(hObject, eventdata, handles)
    global NoiseFile NoiseFilePath
    [NoiseFile NoiseFilePath] = uigetfile('*.avg','Select comparison brainstem response file');
    set(handles.NoiseFileText, 'String', [ NoiseFile]); % this displays 
    global chanN
   
    if isempty(chanN);
    f = openavg([NoiseFilePath NoiseFile]);
    if size(f.chan_names,1) >1
        [chanN v]=  listdlg('PromptString', 'Select Channel to Process' ,'ListString', f.chan_names, 'SelectionMode', 'single') ;
    else
     chanN = 1;
    end
    global chanN
    end
    
     global answerRR
    
    if isempty(answerRR)
        
    
    prompt={'Lag Range Start (ms):                       ','Lag Range Stop (ms):                      '};
    name='Inter-response Correlation Settings';
    numlines = 1;
    defaultanswer={'0','2'};
    
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';


  
        answerRR=inputdlg(prompt,name,numlines,defaultanswer, options);
        answerRR = str2num(char(answerRR));
        global answerRR
    end
               % the file selected to the right of the browse button.
    
function StimFile_Callback(hObject, eventdata, handles)
    global StimFile StimFilePath
    [StimFile StimFilePath] = uigetfile('*.avg','Select filtered stimulus file');
    set(handles.StimFileText, 'String', [StimFile]); % this displays 
    
    global chanStim
    if isempty(chanStim);
    f = openavg([StimFilePath StimFile]);
    if size(f.chan_names,1) >1
        [chanStim v]=  listdlg('PromptString', 'Select Channel to Process' ,'ListString', f.chan_names, 'SelectionMode', 'single') ;
    else
     chanStim = 1;
    end
    global chanStim
    end

    global answerSR
    
    if isempty(answerSR)
        
    
    prompt={'Lag Range Start (ms):','Lag Range Stop (ms):'};
    name='Stimulus-to-Response Correlation Settings';
    numlines = 1;
    defaultanswer={'6.9','9.6'};
  
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';


  
        answerSR=inputdlg(prompt,name,numlines,defaultanswer, options);
        answerSR = str2num(char(answerSR));
        global answerSR
    end
               % the file selected to the right of the browse button.


%% Execute (upon clicking of OK button) ----------------------------------
function okButton_Callback(hObject, eventdata, handles)


global guifig  answerSR answerRR
global Identifier MarkerFileYN MarkerFile NoiseFileYN NoiseFile
global StimFileYN StimFile StartTime StopTime Freq1Start Freq1Stop
global Freq2Start Freq2Stop Freq3Start Freq3Stop FreqStartTime 
global FreqStopTime RRStartTime RRStopTime SRStartTime SRStopTime
global QuietFile QuietFilePath MarkerFilePath XL
global NoiseFilePath StimFilePath FFTScaleYN PlotNoiseFileYN


%% Get all the inputed values. -------------------------------------------
Identifier = get(handles.Identifier, 'String');
MarkerFileYN = get(handles.MarkerFileYN, 'Value');
NoiseFileYN = get(handles.NoiseFileYN, 'Value');
PlotNoiseFileYN = get(handles.PlotNoiseFileYN, 'Value');
StimFileYN = get(handles.StimFileYN, 'Value');
XL = get(handles.XL, 'Value');

StartTime = get(handles.StartTime, 'String'); % RMS
StopTime = get(handles.StopTime, 'String');

FreqStartTime = get(handles.FreqStartTime, 'String'); % FFT time ranges
FreqStopTime = get(handles.FreqStopTime, 'String');


Freq1Start = get(handles.Freq1Start, 'String'); % FFT
Freq1Stop = get(handles.Freq1Stop, 'String');
Freq2Start = get(handles.Freq2Start, 'String');
Freq2Stop = get(handles.Freq2Stop, 'String');
Freq3Start = get(handles.Freq3Start, 'String');
Freq3Stop = get(handles.Freq3Stop, 'String');
FFTScaleYN = get(handles.FFTScaleYN, 'Value');
FreqStartTime = get(handles.FreqStartTime, 'String');
FreqStopTime = get(handles.FreqStopTime, 'String');

RRStartTime = get(handles.RRStartTime, 'String'); % R-R corrs
RRStopTime = get(handles.RRStopTime, 'String');

SRStartTime = get(handles.SRStartTime, 'String'); % S-R corrs
SRStopTime = get(handles.SRStopTime, 'String');

Mark1 = get(handles.Mark1, 'String');
Lat1 = get(handles.Lat1, 'String');
Pos1 = get(handles.Pos1, 'Value');
Mark2 = get(handles.Mark2, 'String');
Lat2 = get(handles.Lat2, 'String');
Pos2 = get(handles.Pos2, 'Value');
Mark3 = get(handles.Mark3, 'String');
Lat3 = get(handles.Lat3, 'String');
Pos3 = get(handles.Pos3, 'Value');
Mark4 = get(handles.Mark4, 'String');
Lat4 = get(handles.Lat4, 'String');
Pos4 = get(handles.Pos4, 'Value');
Mark5 = get(handles.Mark5, 'String');
Lat5 = get(handles.Lat5, 'String');
Pos5 = get(handles.Pos5, 'Value');
Mark6 = get(handles.Mark6, 'String');
Lat6 = get(handles.Lat6, 'String');
Pos6 = get(handles.Pos6, 'Value');
Mark7 = get(handles.Mark7, 'String');
Lat7 = get(handles.Lat7, 'String');
Pos7 = get(handles.Pos7, 'Value');
Mark8 = get(handles.Mark8, 'String');
Lat8 = get(handles.Lat8, 'String');
Pos8 = get(handles.Pos8, 'Value');
Mark9 = get(handles.Mark9, 'String');
Lat9 = get(handles.Lat9, 'String');
Pos9 = get(handles.Pos9, 'Value');
Mark10 = get(handles.Mark10, 'String');
Lat10 = get(handles.Lat10, 'String');
Pos10 = get(handles.Pos10, 'Value');

set(guifig, 'visible', 'off'); 

%% Run BT functions ------------------------------------------------------




global chan
global chanN
global chanStim


%doublecheck that the time ranges are valid.
 f = openavg([QuietFilePath QuietFile]);
 if str2double(StartTime)<f.xmin || str2double(StopTime)>f.xmax
     errordlg('Invalid RMS time tange');
 end
 

  if str2double(FreqStartTime)<f.xmin  ||  str2double(FreqStopTime)>f.xmax
    errordlg('Invalid FFT time range');
 end

 
 if str2double(RRStartTime)<f.xmin|| str2double(RRStopTime)>f.xmax
     errordlg('Invalid inter-response correlation time range');
 end
 

 
 if str2double(SRStartTime)<f.xmin||str2double(SRStopTime)>f.xmax
     errordlg('Invalid stimulus-response correlation time range');
 end
 
 

% 1) RMSes and SNR

 
   [RMS RMSprestim SNR] = bt_rms([QuietFilePath QuietFile], ...
       str2double(StartTime),str2double(StopTime), chan);

% 2) Frequency domain analyses
 

    %   next decide whether to scale and execute:
    if FFTScaleYN == 1
    [Freq1 Freq2 Freq3 FFTwave] = bt_fftsc([QuietFilePath QuietFile],str2double(FreqStartTime),...
        str2double(FreqStopTime),str2double(Freq1Start),str2double(Freq1Stop),...
        str2double(Freq2Start),str2double(Freq2Stop),str2double(Freq3Start),str2double(Freq3Stop), chan);
    else
    [Freq1 Freq2 Freq3 FFTwave] = bt_fft2([QuietFilePath QuietFile],str2double(FreqStartTime),...
        str2double(FreqStopTime),str2double(Freq1Start),str2double(Freq1Stop),...
        str2double(Freq2Start),str2double(Freq2Stop),str2double(Freq3Start),str2double(Freq3Stop), chan);
    end

% 3) inter-response correlations
if NoiseFileYN == 1; % do it or not?
    
    
     % next execute
    % this is a little kludgy:  Using bt_qncorr to compute straight r,
    % (eliminating the 2nd two args) then xcorrelation for cross-r and cross-lag
    RRStraightCorr = bt_qncorr2([QuietFilePath QuietFile], ...
        [NoiseFilePath NoiseFile], str2double(RRStartTime), str2double(RRStopTime),  answerRR(1), answerRR(2), chanN, chan);
    [RRLag RRCorr RRCgram RRlagaxis] = bt_xcorrelation2([NoiseFilePath NoiseFile], ...
        [QuietFilePath QuietFile], str2double(RRStartTime), ...
        str2double(RRStopTime), answerRR(1), answerRR(2), 'POSITIVE',  chanN, chan);
    
    if FFTScaleYN == 1
    [Freq1N Freq2N Freq3N FFTwaveN] = bt_fftsc([NoiseFilePath NoiseFile],str2double(FreqStartTime),...
        str2double(FreqStopTime),str2double(Freq1Start),str2double(Freq1Stop),...
        str2double(Freq2Start),str2double(Freq2Stop),str2double(Freq3Start),str2double(Freq3Stop), chan);
    else
    [Freq1N Freq2N Freq3N FFTwaveN] = bt_fft2([NoiseFilePath NoiseFile],str2double(FreqStartTime),...
        str2double(FreqStopTime),str2double(Freq1Start),str2double(Freq1Stop),...
        str2double(Freq2Start),str2double(Freq2Stop),str2double(Freq3Start),str2double(Freq3Stop), chan);
    end

    
    
end

% 4) stimulus-to-response correlation
if StimFileYN == 1; % do it or not?
    %   Default lag is between 6.6 and 9.6 ms and it seeks a
    %   positive correlation.  
    
    
    
    [SRLag SRCorr SRCgram SRlagaxis] = bt_xcorrelation2([QuietFilePath QuietFile], ...
        [StimFilePath StimFile], str2double(SRStartTime), ...
        str2double(SRStopTime), answerSR(1), answerSR(2), 'POSITIVE', chan, chanStim);
end

% 5) Markers
% Gather markers and assemble into "C" and or "LatsAndPols"
if MarkerFileYN % if marker file
    % if text file
    if strcmp(MarkerFile(end-2:end),'txt')
        % generate cell array "C" from text file
        fid = fopen([MarkerFilePath MarkerFile]);
        C = textscan(fid,'%s%n%s');
        fclose(fid);
        clear fid
    else % if Excel file
        % generate cell array "C" (matching format of above) from xls file
        [Num Txt] = xlsread([MarkerFilePath MarkerFile]);
        C{1} = Txt(:,1);
        C{2} = Num(:,1);
        C{3} = cellstr(num2str(Num(:,2))); % yech
        clear Num
    end
    % generate "LatsAndPols" from C
    LatsAndPols = zeros(10,2); % pad with zeros for when < 10 marks
    LatsAndPols(1:length(C{1}),1) = C{2}; % latencies
    LatsAndPols(1:length(C{1}),2) = str2double(C{3}); % polarities
    % To make code below less complicated, making C{1} always 10x1, with ''
    % padding out the remainder
    for x = size(C{1},1)+1:10
        C{1}{x} = '';
    end
else % if hand-entered
    % generate "LatsAndPols" from entries
    LatsAndPols = [str2double(Lat1), Pos1; str2double(Lat2), Pos2; ...
        str2double(Lat3), Pos3; str2double(Lat4), Pos4; ...
        str2double(Lat5), Pos5; str2double(Lat6), Pos6; ...
        str2double(Lat7), Pos7; str2double(Lat8), Pos8; ...
        str2double(Lat9), Pos9; str2double(Lat10), Pos10];
end

[Picked Auto] = bt_peaks3([QuietFilePath QuietFile], LatsAndPols, chan);



%% Populate bt_gui_out


% open output figure and prepare for being populated
% OutputFig = openfig([cd, '\programFiles\', 'bt_GUI_out.fig']);
OutputFig = openfig(fullfile(cd, 'programFiles', 'bt_GUI_out.fig'));
orient landscape
handles = guihandles(OutputFig);
guidata(OutputFig, handles);
    
% plot waveform(s)
% plot quiet
y = openavg([QuietFilePath QuietFile]);
xaxis = linspace(y.xmin,y.xmax,y.pnts);
% prestim baseline
Qsignal = y.signal(:,chan);
Qprestim = mean(Qsignal(1:ms2row(y,0)));
Qsignal = Qsignal-Qprestim;
% plot
h = findobj(OutputFig, 'tag', 'WavePlot');
if size(xaxis)~=size(Qsignal);
       xaxis = xaxis';
end


plot(h, xax is, Qsignal, 'k'), ylabel(h,'�V'), xlabel(h,'ms') % plot quiet
MaxAmp = max(max(y.signal(:,1)),abs(min(y.signal(:,1)))); % find max voltage
MaxAmp = ceil(MaxAmp.*10)/10; % round up nearest .1 �V
ylim(h, [-MaxAmp MaxAmp])
xlim(h, [xaxis(1) xaxis(end)])
% plot noise
if PlotNoiseFileYN == 1;
    z = openavg([NoiseFilePath NoiseFile]);
    % prestim baseline
    Nsignal = z.signal(:,chanN);
    Nprestim = mean(Nsignal(1:ms2row(z,0)));
    Nsignal = Nsignal-Nprestim;
    % plot
    hold(h,'on'), plot(h, xaxis, Nsignal,'r')
    hold(h,'off') % and back off to prepare for next time
    % legend (only if 2 waveforms are plotted)
    q = legend(h, QuietFile, NoiseFile);
    set(q, 'Interpreter', 'none');
    xlim(h, [xaxis(1) xaxis(end)])
end

% plot marked peaks, black for user-picked; red for auto-picked
% Relying on black to overlay red (based on plot order)
if ~isempty(Picked) % perform if not empty
    if (sum(sum(cell2mat(Picked))))~=0;
    hold(h,'on')
    plot(h, (Auto{1}(Picked{1}~=0)),(Auto{2}(Picked{1}~=0)),'ro')
    plot(h, (Picked{1}(Picked{1}~=0)),(Picked{2}(Picked{1}~=0)),'bo')
    hold(h,'off')
    end
end


% populate identifying info
set(handles.Identifier, 'String', Identifier);
set(handles.QuietFile, 'String', QuietFile);
set(handles.NoiseFile, 'String', NoiseFile);
set(handles.StimFile, 'String', StimFile);
set(handles.MarkerFile, 'String', MarkerFile);

% populate marked peaks table
if sum(sum(cell2mat(Picked))) % only if not all zeros
    % Will reset marker names at end, if marker file.
    % wave 1
    set(handles.Mark1, 'String', Mark1);
    set(handles.P_1lat, 'String', num2str(Picked{1}(1,1),'%3.2f'));
    set(handles.P_1amp, 'String', num2str(Picked{2}(1,1),'%3.3f'));
    set(handles.A_1lat, 'String', num2str(Auto{1}(1,1),'%3.2f'));
    set(handles.A_1amp, 'String', num2str(Auto{2}(1,1),'%3.3f'));
    % wave 2
    set(handles.Mark2, 'String', Mark2);    
    set(handles.P_2lat, 'String', num2str(Picked{1}(2,1),'%3.2f'));
    set(handles.P_2amp, 'String', num2str(Picked{2}(2,1),'%3.3f'));
    set(handles.A_2lat, 'String', num2str(Auto{1}(2,1),'%3.2f'));
    set(handles.A_2amp, 'String', num2str(Auto{2}(2,1),'%3.3f'));
    % wave 3
    set(handles.Mark3, 'String', Mark3);    
    set(handles.P_3lat, 'String', num2str(Picked{1}(3,1),'%3.2f'));
    set(handles.P_3amp, 'String', num2str(Picked{2}(3,1),'%3.3f'));
    set(handles.A_3lat, 'String', num2str(Auto{1}(3,1),'%3.2f'));
    set(handles.A_3amp, 'String', num2str(Auto{2}(3,1),'%3.3f'));
    % wave 4
    set(handles.Mark4, 'String', Mark4);    
    set(handles.P_4lat, 'String', num2str(Picked{1}(4,1),'%3.2f'));
    set(handles.P_4amp, 'String', num2str(Picked{2}(4,1),'%3.3f'));
    set(handles.A_4lat, 'String', num2str(Auto{1}(4,1),'%3.2f'));
    set(handles.A_4amp, 'String', num2str(Auto{2}(4,1),'%3.3f'));
    % wave 5
    set(handles.Mark5, 'String', Mark5);    
    set(handles.P_5lat, 'String', num2str(Picked{1}(5,1),'%3.2f'));
    set(handles.P_5amp, 'String', num2str(Picked{2}(5,1),'%3.3f'));
    set(handles.A_5lat, 'String', num2str(Auto{1}(5,1),'%3.2f'));
    set(handles.A_5amp, 'String', num2str(Auto{2}(5,1),'%3.3f'));
    % wave 6
    set(handles.Mark6, 'String', Mark6);    
    set(handles.P_6lat, 'String', num2str(Picked{1}(6,1),'%3.2f'));
    set(handles.P_6amp, 'String', num2str(Picked{2}(6,1),'%3.3f'));
    set(handles.A_6lat, 'String', num2str(Auto{1}(6,1),'%3.2f'));
    set(handles.A_6amp, 'String', num2str(Auto{2}(6,1),'%3.3f'));
    % wave 7
    set(handles.Mark7, 'String', Mark7);    
    set(handles.P_7lat, 'String', num2str(Picked{1}(7,1),'%3.2f'));
    set(handles.P_7amp, 'String', num2str(Picked{2}(7,1),'%3.3f'));
    set(handles.A_7lat, 'String', num2str(Auto{1}(7,1),'%3.2f'));
    set(handles.A_7amp, 'String', num2str(Auto{2}(7,1),'%3.3f'));
    % wave 8
    set(handles.Mark8, 'String', Mark8);    
    set(handles.P_8lat, 'String', num2str(Picked{1}(8,1),'%3.2f'));
    set(handles.P_8amp, 'String', num2str(Picked{2}(8,1),'%3.3f'));
    set(handles.A_8lat, 'String', num2str(Auto{1}(8,1),'%3.2f'));
    set(handles.A_8amp, 'String', num2str(Auto{2}(8,1),'%3.3f'));
    % wave 9
    set(handles.Mark9, 'String', Mark9);    
    set(handles.P_9lat, 'String', num2str(Picked{1}(9,1),'%3.2f'));
    set(handles.P_9amp, 'String', num2str(Picked{2}(9,1),'%3.3f'));
    set(handles.A_9lat, 'String', num2str(Auto{1}(9,1),'%3.2f'));
    set(handles.A_9amp, 'String', num2str(Auto{2}(9,1),'%3.3f'));
    % wave 10
    set(handles.Mark10, 'String', Mark10);    
    set(handles.P_10lat, 'String', num2str(Picked{1}(10,1),'%3.2f'));
    set(handles.P_10amp, 'String', num2str(Picked{2}(10,1),'%3.3f'));
    set(handles.A_10lat, 'String', num2str(Auto{1}(10,1),'%3.2f'));
    set(handles.A_10amp, 'String', num2str(Auto{2}(10,1),'%3.3f'));
    % mark names if marker file
    if MarkerFileYN == 1
            set(handles.Mark1, 'String', C{1}(1));
            set(handles.Mark2, 'String', C{1}(2));
            set(handles.Mark3, 'String', C{1}(3));
            set(handles.Mark4, 'String', C{1}(4));
            set(handles.Mark5, 'String', C{1}(5));
            set(handles.Mark6, 'String', C{1}(6));
            set(handles.Mark7, 'String', C{1}(7));
            set(handles.Mark8, 'String', C{1}(8));
            set(handles.Mark9, 'String', C{1}(9));
            set(handles.Mark10, 'String', C{1}(10));
    end
end

% populate RMS section
    set(handles.SNR, 'String', num2str(SNR,'%2.2f'));
    set(handles.RMS, 'String', num2str(RMS,'%2.3f'));
    set(handles.RMSprestim, 'String', num2str(RMSprestim,'%2.3f'));
    set(handles.RMStime, 'String', [StartTime '-' StopTime]);

% populate frequency analysis section
    set(handles.FreqTime, 'String', [FreqStartTime '-' FreqStopTime]);
    set(handles.Freq1Range, 'String', [Freq1Start '-' Freq1Stop]);
    set(handles.Freq2Range, 'String', [Freq2Start '-' Freq2Stop]);
    set(handles.Freq3Range, 'String', [Freq3Start '-' Freq3Stop]);    
    set(handles.Freq1, 'String', num2str(Freq1,'%2.3f'));
    set(handles.Freq2, 'String', num2str(Freq2,'%2.3f'));
    set(handles.Freq3, 'String', num2str(Freq3,'%2.3f'));    

% populate correlational analyses section
if StimFileYN == 1;
    set(handles.SRTime, 'String', [SRStartTime '-' SRStopTime]);
    set(handles.SRCorr, 'String', num2str(SRCorr,'%2.2f'));
    set(handles.SRLag, 'String', num2str(SRLag,'%2.2f'));
end

if NoiseFileYN == 1;
    set(handles.RRTime, 'String', [RRStartTime '-' RRStopTime]);
    set(handles.RRStraightCorr, 'String', num2str(RRStraightCorr,'%2.2f'));
    set(handles.RRCorr, 'String', num2str(RRCorr,'%2.2f'));
    set(handles.RRLag, 'String', num2str(RRLag,'%2.2f'));
end

% Other plots: FFTs and Cgrams
h = findobj(OutputFig, 'tag', 'FreqPlot');
xaxis = (0:1:length(FFTwave)-1)';
plot(h, xaxis, FFTwave, 'k'), xlabel(h,'Hz'), xlim(h,[0 1500])

if PlotNoiseFileYN == 1;
   hold(h,'on')
    plot(h, xaxis, FFTwaveN, 'r'), xlabel(h,'Hz'), xlim(h,[0 1500]);
end

ylabel(h, 'peak �V')
h = findobj(OutputFig, 'tag', 'FreqPlotZoom');
xaxis = (0:1:length(FFTwave)-1)';
plot(h, xaxis, FFTwave, 'k'), xlabel(h,'Hz'), xlim(h,[0 400])

if PlotNoiseFileYN == 1;
    hold(h,'on')
    plot(h, xaxis, FFTwaveN, 'r'), xlabel(h,'Hz'), xlim(h,[0 400])
end

if StimFileYN == 1;
    h = findobj(OutputFig, 'tag', 'SRCgram');
    plot(h, SRlagaxis, SRCgram, 'k')
    hold(h,'on')
    plot(h, SRLag,SRCorr,'r*')
    hold(h,'off')
    xlim(h, [SRlagaxis(1) SRlagaxis(end)])
    ylim([-1.01 1.01])
end

if NoiseFileYN == 1;
    h = findobj(OutputFig, 'tag', 'RRCgram');
    plot(h, RRlagaxis, RRCgram, 'k'), xlabel(h,'lag (ms)')
    hold(h,'on')
    plot(h, RRLag,RRCorr,'r*')
    hold(h,'off')
    xlim(h, [RRlagaxis(1) RRlagaxis(end)])
    ylim([-1.01 1.01])
end

%% Save to Excel

if XL == 1
    % align output: identifiying info and labels in cell array, data in num array.
    ID{1} = Identifier;
    ID{2} = ['Primary file: ' QuietFile];
    ID{3} = ['Comparison file: ' NoiseFile];
    ID{4} = ['Stimulus File: ' StimFile];
    ID{5} = ['Marker File: ' MarkerFile];
    Labels = {...
        'FFRTimeStart','FFRTimeStop','ResponseRMS','PrestimRMS','SNR',... %1-5
        'Freq1Start','Freq1Stop','FreqAmp1', ...
        'Freq2Start','Freq2Stop','FreqAmp2', ...
        'Freq3Start','Freq3Stop','FreqAmp3', ... %6-14
        'SRTimeStart','SRTimeStop','SRr','SRlag', ... 
        'RRTimeStart','RRTimeStop','RRunshiftedr','RRr','RRlag', ... %15-23
        'Mark1LatPicked','Mark1AmpPicked', ...
          'Mark1LatAuto','Mark1AmpAuto', ... %24-27
        'Mark2LatPicked','Mark2AmpPicked', ...
          'Mark2LatAuto','Mark2AmpAuto', ... %28-31
        'Mark3LatPicked','Mark3AmpPicked', ...
          'Mark3LatAuto','Mark3AmpAuto', ... %32-35
        'Mark4LatPicked','Mark4AmpPicked', ...
          'Mark4LatAuto','Mark4AmpAuto', ... %36-39
        'Mark5LatPicked','Mark5AmpPicked', ...
          'Mark5LatAuto','Mark5AmpAuto', ... %40-43
        'Mark6LatPicked','Mark6AmpPicked', ...
          'Mark6LatAuto','Mark6AmpAuto', ... %44-47
        'Mark7LatPicked','Mark7AmpPicked', ...
          'Mark7LatAuto','Mark7AmpAuto', ... %48-51
        'Mark8LatPicked','Mark8AmpPicked', ...
          'Mark8LatAuto','Mark8AmpAuto', ... %52-55
        'Mark9LatPicked','Mark9AmpPicked', ...
          'Mark9LatAuto','Mark9AmpAuto', ... %56-59
        'Mark10LatPicked','Mark10AmpPicked', ...
          'Mark10LatAuto','Mark10AmpAuto'};  %60-63

    % Exported data will contain -999 for N/A cells.
    Data = ones(1,63).*-999;

    % Populate "Data" with analysis results
    % First the marks, loopable, but vaguely tricky.
    if sum(Picked{1}~=0) % if Picked latency contains anything but zeros
        for x = 1:sum(Picked{1}~=0) % # of non-zeros
            Data((x-1).*4+24) = Picked{1}(x);
            Data((x-1).*4+25) = Picked{2}(x);            
            Data((x-1).*4+26) = Auto{1}(x);            
            Data((x-1).*4+27) = Auto{2}(x);                        
        end
    end

    % Now march thru the rest, tedious but straighforward.
    Data(1) = str2double(StartTime);
    Data(2) = str2double(StopTime);
    Data(3) = RMS;
    Data(4) = RMSprestim;
    Data(5) = SNR;
    Data(6) = str2double(Freq1Start);
    Data(7) = str2double(Freq1Stop);
    Data(8) = Freq1;
    Data(9) = str2double(Freq2Start);
    Data(10) = str2double(Freq2Stop);
    Data(11) = Freq2;
    Data(12) = str2double(Freq3Start);
    Data(13) = str2double(Freq3Stop);
    Data(14) = Freq3;
    if StimFileYN == 1;
        Data(15) = str2double(SRStartTime);
        Data(16) = str2double(SRStopTime);
        Data(17) = SRCorr;
        Data(18) = SRLag;
    end
    if NoiseFileYN == 1;
        Data(19) = str2double(RRStartTime);
        Data(20) = str2double(RRStopTime);
        Data(21) = RRStraightCorr;
        Data(22) = RRCorr;
        Data(23) = RRLag;
    end

    if Identifier
        bt_xlswrite(Data,ID,Labels,fullfile(cd, 'outputFiles', Identifier),'Sheet1');
    else
        bt_xlswrite(Data,ID,Labels,fullfile(cd, 'outputFiles', QuietFile(1:end-4)),'Sheet1');
    end
end



