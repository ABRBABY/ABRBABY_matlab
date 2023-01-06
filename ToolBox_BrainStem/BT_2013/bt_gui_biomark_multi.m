function varargout = bt_gui_BioMARK(varargin)
% The brainstem toolbox  is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation.
%
%~ The brainstem toolbox is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%~ See the GNU General Public License for more details <http://www.gnu.org/licenses/>/
%	Copyright 2007-2008, Trent Nicol and Erika Skoe.
% Updated, Jan 2010 by Erika Skoe to fix MATLAB2010-related bug.
% 
% Dependencies:
% brainstem toolbox m-files: bt_ascii2avg, bt_fft2, bt_fftsc,
%          bt_peaks2, bt_qncorr, bt_rms, bt_xlswrite
% other m-files: closestrc, ms2row, nancorrcoef, openavg, xcorrelation (at
%          least)
% fig-files: bt_gui_BioMARK, bt_gui_BioMARK_out
% avg-files: BioMARK-Template, Filtered_BioMARK_stimulus.avg
%

addpath([cd, '\programFiles'])
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @bt_gui_BioMARK_OpeningFcn, ...
    'gui_OutputFcn',  @bt_gui_BioMARK_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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

% --- Executes just before bt_gui_BioMARK is made visible.
function bt_gui_BioMARK_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bt_gui_BioMARK (see VARARGIN)
% guifig = openfig([cd, '\programFiles\bt_gui_BioMARK.fig'],'reuse');   % the same window is used each time within 1 matlab session.
%                                          % this allows path information to be retained from subject to subject
% set(guifig, 'visible', 'on');
% Choose default command line output for bt_gui_BioMARK
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = bt_gui_BioMARK_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------END OF AUTO-GENERATED CODE-------------------
%% set some defaults

set(handles.NoiseFileYN, 'Value', 0);
set(handles.NoiseFile, 'Enable', 'off');
set(handles.PlotNoiseFileYN, 'Value', 0);
set(handles.StimFileYN, 'Value', 0);
set(handles.StimFile, 'Enable', 'off');
set(handles.StimFileText, 'string', 'filtered_BioMARK_stimulus.avg');
global StimFile StimFilePath
StimFilePath = [cd, '\ProgramFiles\'];
StimFile =  'filtered_BioMARK_stimulus.avg';

set(handles.StartTime, 'String', 11.4);
set(handles.StopTime, 'String', 40.6);

set(handles.Freq1Start, 'String', 103);
set(handles.Freq1Stop, 'String', 121);
set(handles.Freq2Start, 'String', 454);
set(handles.Freq2Stop, 'String', 719);
set(handles.Freq3Start, 'String', 721);
set(handles.Freq3Stop, 'String', 1155);
set(handles.FFTScaleYN, 'Value', 0);

set(handles.FreqStartTime, 'String', 11.4);
set(handles.FreqStopTime, 'String', 40.6);
set(handles.RRStartTime, 'String', 11.4);
set(handles.RRStopTime, 'String',40.6);
set(handles.SRStartTime, 'String', 10);
set(handles.SRStopTime, 'String', 40);
set(handles.ExcelFileYN, 'Value', 1);

%% Get data (checkboxes, browse buttons)

% Checkboxes
% Noise file
function NoiseFileYN_Callback(hObject, eventdata, handles)
NoiseFileYN = get(handles.NoiseFileYN, 'Value');
% conditional enabling of corresponding browse button
if NoiseFileYN == 1;
    set(handles.NoiseFile, 'Enable', 'on');
else
    set(handles.NoiseFile, 'Enable', 'off');
end

% Stimulus file
function StimFileYN_Callback(hObject, eventdata, handles)
StimFileYN = get(handles.StimFileYN, 'Value');
% conditional enabling of corresponding browse button
if StimFileYN == 1;
    set(handles.StimFile, 'Enable', 'on');
else
    set(handles.StimFile, 'Enable', 'off');
end

function QuietFile_Callback(hObject, eventdata, handles)

global QuietFileTxt FilePath
[QuietFileTxt FilePath] = uigetfile('*.txt','Select quiet BioMAP ASCII file(s)', 'MultiSelect', 'on');  %updated to be MATLAB2010 compatible



if iscell(QuietFileTxt) == 0  %Handles the case when only 1 file is selected.
    temp = QuietFileTxt;
    clear QuietFileTxt;
    QuietFileTxt{1} = temp;
end

set(handles.QuietFileText, 'String', [num2str(length(QuietFileTxt)) ' file(s) selected']); % this displays

% Noise file browse
function NoiseFile_Callback(hObject, eventdata, handles)
global NoiseFileTxt NoiseFilePath
[NoiseFileTxt NoiseFilePath] = uigetfile('*.txt','Select noise BioMAP ASCII file(s)', 'MultiSelect', 'on'); %updated to be MATLAB2010 compatible



if iscell(NoiseFileTxt) == 0  %Handles the case when only 1 file is selected.
    temp = NoiseFileTxt;
    clear NoiseFileTxt;
    NoiseFileTxt{1} = temp;
end

set(handles.NoiseFileText, 'String', [num2str(length(NoiseFileTxt)) ' file(s) selected']); % this displays

% Stimulus file browse
function StimFile_Callback(hObject, eventdata, handles)
global StimFile StimFilePath
[StimFile StimFilePath] = uigetfile('*.avg','Select filtered stimulus file');
set(handles.StimFileText, 'String', [ StimFile]); % this displays
% the file selected to the right of the browse button.



%% Execute (upon clicking of OK button) ----------------------------------
function okButton_Callback(hObject, eventdata, handles)

global guifig
global MarkerFile NoiseFileYN NoiseFile NoiseFileTxt
global StimFileYN StimFile StartTime StopTime Freq1Start
global Freq1Stop Freq2Start Freq2Stop Freq3Start Freq3Stop FreqStartTime
global FreqStopTime RRStartTime RRStopTime SRStartTime SRStopTime
global QuietFile QuietFileTxt FilePath MarkerFilePath ExcelFileYN
global NoiseFilePath StimFilePath FFTScaleYN PlotNoiseFileYN OutputFig


if iscell(QuietFileTxt) == 0  %Handles the case when only 1 file is selected.
    temp = QuietFileTxt;
    clear QuietFileTxt;
    QuietFileTxt{1} = temp;
end


if iscell(NoiseFileTxt) == 0  %Handles the case when only 1 file is selected.
    temp = NoiseFileTxt;
    clear NoiseFileTxt;
    NoiseFileTxt{1} = temp;
end

%% Get all the inputed values. -------------------------------------------
for X = 1:length(QuietFileTxt)
    if X == 1;

        NoiseFileYN = get(handles.NoiseFileYN, 'Value');
        PlotNoiseFileYN = get(handles.PlotNoiseFileYN, 'Value');
        StimFileYN = get(handles.StimFileYN, 'Value');
        ExcelFileYN = get(handles.ExcelFileYN, 'Value');

        StartTime = get(handles.StartTime, 'String'); %RMS plus default for others
        StopTime = get(handles.StopTime, 'String'); %RMS plus default for others

        Freq1Start = get(handles.Freq1Start, 'String');
        Freq1Stop = get(handles.Freq1Stop, 'String');
        Freq2Start = get(handles.Freq2Start, 'String');
        Freq2Stop = get(handles.Freq2Stop, 'String');
        Freq3Start = get(handles.Freq3Start, 'String');
        Freq3Stop = get(handles.Freq3Stop, 'String');
        FFTScaleYN = get(handles.FFTScaleYN, 'Value');

        FreqStartTime = get(handles.FreqStartTime, 'String');
        FreqStopTime = get(handles.FreqStopTime, 'String');

        RRStartTime = get(handles.RRStartTime, 'String');
        RRStopTime = get(handles.RRStopTime, 'String');

        SRStartTime = get(handles.SRStartTime, 'String');
        SRStopTime = get(handles.SRStopTime, 'String');

        set(guifig, 'visible', 'off');
    end
    %% Run BT functions ------------------------------------------------------


    % 0) Convert .txt file to corresponding .avg and .mrk files, and assign
    %    them variable names
    bt_biologic2avg([FilePath QuietFileTxt{X}])
    QuietFile = [QuietFileTxt{X}(1:end-4) '.avg'];
    MarkerFile = [QuietFileTxt{X}(1:end-4) '.mrk'];


    movefile([FilePath QuietFile], [cd, '\outputFiles\', QuietFile]);
    movefile([FilePath MarkerFile], [cd, '\outputFiles\', MarkerFile]);

    QuietFile = [cd, '\outputFiles\', QuietFile];
    MarkerFile = [cd, '\outputFiles\', MarkerFile];

    [p QuietFileName ext] = fileparts(QuietFile) ;
    [p MarkerFileName ext] = fileparts(MarkerFile) ;

    chan = 1;
    if NoiseFileYN == 1; % do it or not?
        bt_biologic2avg([NoiseFilePath NoiseFileTxt{X}])
        NoiseFile = [NoiseFileTxt{X}(1:end-4) '.avg'];
        movefile([NoiseFilePath NoiseFile], [cd, '\outputFiles\', NoiseFile]);
        NoiseFile = [cd, '\outputFiles\', NoiseFile];
        [p NoiseFileName ext] = fileparts(NoiseFile) ;
    else
        NoiseFileName = ''; % to avoid ghost filename appearing on output
    end




    % 1) RMSes and SNR
    [RMS RMSprestim SNR] = bt_rms([QuietFile], ...
        bt_str2double(StartTime),bt_str2double(StopTime), chan);

    % 2) Frequency domain analyses
    %   First, determine time range:


    %   next decide whether to scale and execute:
    if X == 1
        FreqStopTime = bt_str2double(FreqStopTime);
        FreqStartTime = bt_str2double(FreqStartTime);

        if (FreqStopTime-FreqStartTime)<4
            FreqStopTime=FreqStartTime +4;
            display(['FFT time range too small: user inputted FFT stop time reset to ' num2str(FreqStopTime)]);
        end
    end

    if FFTScaleYN == 1
        [Freq1 Freq2 Freq3 FFTwave xaxisHz] = bt_fftsc([QuietFile],FreqStartTime,...
            FreqStopTime,bt_str2double(Freq1Start),bt_str2double(Freq1Stop),...
            bt_str2double(Freq2Start),bt_str2double(Freq2Stop),bt_str2double(Freq3Start),bt_str2double(Freq3Stop), chan);
    else
        [Freq1 Freq2 Freq3 FFTwave xaxisHz] = bt_fftBioMARK([QuietFile],FreqStartTime,...
            FreqStopTime,bt_str2double(Freq1Start),bt_str2double(Freq1Stop),...
            bt_str2double(Freq2Start),bt_str2double(Freq2Stop),bt_str2double(Freq3Start),bt_str2double(Freq3Stop), 0, chan);
    end

    % 3) quiet-to-noise correlations
    if NoiseFileYN == 1; % do it or not?

        % next execute
        % this is a little kludgy:  Using bt_qncorr to compute straight r,
        % eliminating the 2nd two args) then xcorrelation for cross-r and cross-lag
        RRStraightCorr = bt_qncorr([QuietFile], ...
            [ NoiseFile], bt_str2double(RRStartTime), bt_str2double(RRStopTime), chan);
        [RRLag RRCorr RRCgram RRlagaxis] = bt_xcorrelation([NoiseFile], ...
            [QuietFile], bt_str2double(RRStartTime), ...
            bt_str2double(RRStopTime), 0, 2, 'POSITIVE', chan);


        if FFTScaleYN == 1
            s = 1;
            [Freq1N Freq2N Freq3N FFTwaveN xaxisHz] = bt_fftsc([NoiseFile],FreqStartTime,...
                FreqStopTime,bt_str2double(Freq1Start),bt_str2double(Freq1Stop),...
                bt_str2double(Freq2Start),bt_str2double(Freq2Stop),bt_str2double(Freq3Start),bt_str2double(Freq3Stop), chan);
        else
            s = 0;
            [Freq1N Freq2N Freq3N FFTwaveN xaxisHz] =  bt_fftBioMARK([NoiseFile],FreqStartTime,...
                FreqStopTime,bt_str2double(Freq1Start),bt_str2double(Freq1Stop),...
                bt_str2double(Freq2Start),bt_str2double(Freq2Stop),bt_str2double(Freq3Start),bt_str2double(Freq3Stop), 0, chan);




        end
    end

    % 4) stimulus-to-response correlation
    % next execute.  Default lag is between 6.6 and 9.6 ms and it
    % seeks a positive correlation.  Lag range and sign can be changed by user, below.
    [SRLag SRCorr SRCgram SRlagaxis] = bt_xcorrelation([QuietFile], ...
        [StimFilePath StimFile], bt_str2double(SRStartTime), ...
        bt_str2double(SRStopTime), 6.6, 9.6, 'POSITIVE', chan);


    % 5) Marker report
    [Picked Auto] = bt_peaks2([QuietFile],[MarkerFile]);


    %% Populate bt_gui_BioMARK_out
    % open output figure and prepare for being populated
    OutputFig = openfig([cd, '\programFiles\', 'bt_gui_BioMARK_out.fig']);
    handles = guihandles(OutputFig);
    orient landscape
    guidata(OutputFig, handles);

    % plot waveform(s)
    % plot quiet
    y = openavg([QuietFile]);
    xaxis = linspace(y.xmin,y.xmax,y.pnts);

    % prestim baseline
    Qsignal = y.signal;
    Qprestim = mean(Qsignal(1:ms2row(y,0)));
    Qsignal = Qsignal-Qprestim;

    if size(xaxis,1)~=size(Qsignal,1)
        xaxis = xaxis';
    end
    % plot
    h = findobj(OutputFig, 'tag', 'WavePlot');

    plot(h, xaxis, Qsignal, 'color', 'k'), ylabel(h,'µV'), xlabel(h,'ms') % plot quiet
    set(h, 'xlim' ,[xaxis(1) xaxis(end)]);

    MaxAmp = max(max(y.signal),abs(min(y.signal))); % find max voltage
    MaxAmp = ceil(MaxAmp.*10)/10; % round up nearest .1 µV
    ylim(h, [-MaxAmp MaxAmp])
    % legend(h, QuietFile)
    % plot noise
    if PlotNoiseFileYN == 1;
        z = openavg([ NoiseFile]);
        % prestim baseline
        Nsignal = z.signal;
        Nprestim = mean(Nsignal(1:ms2row(z,0)));
        Nsignal = Nsignal-Nprestim;
        % plot
        hold(h,'on'), plot(h, xaxis, Nsignal,'r')
        hold(h,'off') % and back off to prepare for next time

        l = legend(h, QuietFileName, NoiseFileName);
        set(l, 'Interpreter', 'None');
    end

    % populate identifying info

    set(handles.QuietFile, 'String', QuietFileName);
    set(handles.NoiseFile, 'String', NoiseFileName);
    set(handles.StimFile, 'String', StimFile);
    set(handles.MarkerFile, 'String', MarkerFileName);

    % populate marked peaks table (and plot marks)
    if exist('Picked','var')

        % Set Autos to NaN if not different from Picked
        Auto{1}(Auto{1}==Picked{1}) = NaN;
        Auto{2}(Auto{2}==Picked{2}) = NaN;
        Auto{3}(Auto{3}==Picked{3}) = NaN;

        if size(Picked,1) == 0
            NumMarks = 0;
        else
            NumMarks = size(Picked{1},1);
        end

        if NumMarks == 7
            % wave O
            set(handles.P_Olat, 'String', num2str(Picked{1}(7,1),'%3.2f'));
            set(handles.P_Oamp, 'String', num2str(Picked{2}(7,1),'%3.3f'));
            set(handles.A_Olat, 'String', num2str(Auto{1}(7,1),'%3.2f'));
            set(handles.A_Oamp, 'String', num2str(Auto{2}(7,1),'%3.3f'));
        end
        if NumMarks > 5
            % wave F
            set(handles.P_Flat, 'String', num2str(Picked{1}(6,1),'%3.2f'));
            set(handles.P_Famp, 'String', num2str(Picked{2}(6,1),'%3.3f'));
            set(handles.A_Flat, 'String', num2str(Auto{1}(6,1),'%3.2f'));
            set(handles.A_Famp, 'String', num2str(Auto{2}(6,1),'%3.3f'));
        end
        if NumMarks > 4
            % wave E
            set(handles.P_Elat, 'String', num2str(Picked{1}(5,1),'%3.2f'));
            set(handles.P_Eamp, 'String', num2str(Picked{2}(5,1),'%3.3f'));
            set(handles.A_Elat, 'String', num2str(Auto{1}(5,1),'%3.2f'));
            set(handles.A_Eamp, 'String', num2str(Auto{2}(5,1),'%3.3f'));
        end
        if NumMarks > 3
            % wave D
            set(handles.P_Dlat, 'String', num2str(Picked{1}(4,1),'%3.2f'));
            set(handles.P_Damp, 'String', num2str(Picked{2}(4,1),'%3.3f'));
            set(handles.A_Dlat, 'String', num2str(Auto{1}(4,1),'%3.2f'));
            set(handles.A_Damp, 'String', num2str(Auto{2}(4,1),'%3.3f'));
        end
        if NumMarks > 2
            % wave C
            set(handles.P_Clat, 'String', num2str(Picked{1}(3,1),'%3.2f'));
            set(handles.P_Camp, 'String', num2str(Picked{2}(3,1),'%3.3f'));
            set(handles.A_Clat, 'String', num2str(Auto{1}(3,1),'%3.2f'));
            set(handles.A_Camp, 'String', num2str(Auto{2}(3,1),'%3.3f'));
        end
        if NumMarks > 1
            % wave A
            set(handles.P_Alat, 'String', num2str(Picked{1}(2,1),'%3.2f'));
            set(handles.P_Aamp, 'String', num2str(Picked{2}(2,1),'%3.3f'));
            set(handles.A_Alat, 'String', num2str(Auto{1}(2,1),'%3.2f'));
            set(handles.A_Aamp, 'String', num2str(Auto{2}(2,1),'%3.3f'));
            % V/A interpeak lat
            set(handles.P_VAlat, 'String', num2str(Picked{3}(1,1),'%3.2f'));
            set(handles.A_VAlat, 'String', num2str(Auto{3}(1,1),'%3.2f'));
            % V/A interpeak amp
            set(handles.P_VAamp, 'String', num2str(Picked{3}(2,1),'%3.3f'));
            set(handles.A_VAamp, 'String', num2str(Auto{3}(2,1),'%3.3f'));
            % V/A slope
            set(handles.P_slope, 'String', num2str(Picked{3}(3,1),'%3.2f'));
            set(handles.A_slope, 'String', num2str(Auto{3}(3,1),'%3.2f'));
            % V/A area
            set(handles.P_area, 'String', num2str(Picked{3}(4,1),'%3.3f'));
            set(handles.A_area, 'String', num2str(Auto{3}(4,1),'%3.3f'));
        end
        if NumMarks > 0
            % wave V
            set(handles.P_Vlat, 'String', num2str(Picked{1}(1,1),'%3.2f'));
            set(handles.P_Vamp, 'String', num2str(Picked{2}(1,1),'%3.3f'));
            set(handles.A_Vlat, 'String', num2str(Auto{1}(1,1),'%3.2f'));
            set(handles.A_Vamp, 'String', num2str(Auto{2}(1,1),'%3.3f'));
            hold(h,'on')
            plot(h, cell2mat(Auto(1)),cell2mat(Auto(2)),'ro')
            plot(h, cell2mat(Picked(1)),cell2mat(Picked(2)),'bo')
            hold(h,'off')
        end
    end

    % populate RMS section
    set(handles.SNR, 'String', num2str(SNR,'%2.2f'));
    set(handles.RMS, 'String', num2str(RMS,'%2.3f'));
    set(handles.RMSprestim, 'String', num2str(RMSprestim,'%2.3f'));
    set(handles.RMStime, 'String', [StartTime '-' StopTime]);

    % populate frequency analysis section
    set(handles.FreqTime, 'String', [num2str(FreqStartTime) '-' num2str(FreqStopTime)]);
    set(handles.Freq1Range, 'String', [Freq1Start '-' Freq1Stop]);
    set(handles.Freq2Range, 'String', [Freq2Start '-' Freq2Stop]);
    set(handles.Freq3Range, 'String', [Freq3Start '-' Freq3Stop]);
    set(handles.Freq1, 'String', num2str(Freq1,'%2.3f'));
    set(handles.Freq2, 'String', num2str(Freq2,'%2.3f'));
    set(handles.Freq3, 'String', num2str(Freq3,'%2.3f'));

    % populate correlational analyses section
    set(handles.SRTime, 'String', [SRStartTime '-' SRStopTime]);
    set(handles.SRCorr, 'String', num2str(SRCorr,'%2.2f'));
    set(handles.SRLag, 'String', num2str(SRLag,'%2.2f'));

    if NoiseFileYN == 1;
        set(handles.RRTime, 'String', [RRStartTime '-' RRStopTime]);
        set(handles.RRStraightCorr, 'String', num2str(RRStraightCorr,'%2.2f'));
        set(handles.RRCorr, 'String', num2str(RRCorr,'%2.2f'));
        set(handles.RRLag, 'String', num2str(RRLag,'%2.2f'));
    end

    % Other plots: FFTs and Cgrams
    h = findobj(OutputFig, 'tag', 'FreqPlot');

    plot(h, xaxisHz, FFTwave, 'k'), xlabel(h,'Hz'), xlim(h,[0 1500])

    if PlotNoiseFileYN == 1;

        hold(h,'on');
        plot(h, xaxisHz, FFTwaveN, 'r'), xlabel(h,'Hz'), xlim(h,[0 1500]);
    end


    ylabel(h, 'peak µV')
    h = findobj(OutputFig, 'tag', 'FreqPlotZoom');

    plot(h, xaxisHz, FFTwave, 'k'), xlabel(h,'Hz'), xlim(h,[0 400])

    if PlotNoiseFileYN == 1;
        hold(h,'on');
        plot(h, xaxisHz, FFTwaveN, 'r'), xlabel(h,'Hz'), xlim(h,[0 400])
    end


    h = findobj(OutputFig, 'tag', 'SRCgram');
    plot(h, SRlagaxis, SRCgram, 'k')
    hold(h,'on')
    plot(h, SRLag,SRCorr,'r*')
    hold(h,'off')
    set(h, 'xlim', [SRlagaxis(1) SRlagaxis(end)]);

    h = findobj(OutputFig, 'tag', 'RRCgram');
    xlabel(h,'lag (ms)')
    if NoiseFileYN == 1;
        plot(h, RRlagaxis, RRCgram, 'k')
        hold(h,'on')
        plot(h, RRLag,RRCorr,'r*')
        hold(h,'off')
        set(h, 'xlim', [RRlagaxis(1) RRlagaxis(end)]);
    end

    %% Save to Excel

    if ExcelFileYN == 1
        % align output: identifiying info and labels in cell array, data in num array.

        ID{1} = ['QuietFile: ' QuietFileName];
        ID{2} = ['NoiseFile: ' NoiseFileName];
        ID{3} = ['StimFile: ' StimFile];

        Labels = {...
            'FileName', 'VlatPicked','AlatPicked','ClatPicked','DlatPicked','ElatPicked','FlatPicked',...
            'OlatPicked',...  %1-7
            'VampPicked','AampPicked','CampPicked','DampPicked','EampPicked',...
            'FampPicked','OampPicked',...  %8-14
            'VAlatPicked','VAampPicked','VAslopePicked','VAareaPicked',... %15-18
            'VlatAuto','AlatAuto','ClatAuto','DlatAuto','ElatAuto','FlatAuto',...
            'OlatAuto',... %19-25
            'VampAuto','AampAuto','CampAuto','DampAuto','EampAuto',...
            'FampAuto','OampAuto',... %26-32
            'VAlatAuto','VAampAuto','VAslopeAuto','VAareaAuto',... %33-36
            'FFRTimeStart','FFRTimeStop','ResponseRMS','PrestimRMS','SNR',... %37-41
            'FreqTimeStart','FreqTimeStop','Freq1Start','Freq1Stop','FreqAmp1',...
            'Freq2Start','Freq2Stop','FreqAmp2','Freq3Start','Freq3Stop','FreqAmp3',... %42-52
            'SRTimeStart','SRTimeStop','SRr','SRlag',... %53-56
            'RRTimeStart','RRTimeStop','RRunshiftedr','RRr','RRlag'}; %57-61

        % Exported data will contain -999 for N/A cells.
        Data = ones(1,61).*-999;

        % Populate "Data" with analysis results
        % First the marks, loopable, but vaguely tricky.
        if NumMarks
            for x = 1:NumMarks
                Data(x) = Picked{1}(x);
                Data(x+7) = Picked{2}(x);
                Data(x+18) = Auto{1}(x);
                Data(x+25) = Auto{2}(x);
            end
        end
        if NumMarks > 1
            for x = 1:4
                Data(x+14) = Picked{3}(x);
                Data(x+32) = Auto{3}(x);
            end
        end

        % Now march thru the rest, tedious but straighforward.
        Data(37) = bt_str2double(StartTime);
        Data(38) = bt_str2double(StopTime);
        Data(39) = RMS;
        Data(40) = RMSprestim;
        Data(41) = SNR;
        Data(42) = FreqStartTime;
        Data(43) = FreqStopTime;
        Data(44) = bt_str2double(Freq1Start);
        Data(45) = bt_str2double(Freq1Stop);
        Data(46) = Freq1;
        Data(47) = bt_str2double(Freq2Start);
        Data(48) = bt_str2double(Freq2Stop);
        Data(49) = Freq2;
        Data(50) = bt_str2double(Freq3Start);
        Data(51) = bt_str2double(Freq3Stop);
        Data(52) = Freq3;
        Data(53) = bt_str2double(SRStartTime);
        Data(54) = bt_str2double(SRStopTime);
        Data(55) = SRCorr;
        Data(56) = SRLag;
        if NoiseFileYN == 1;
            Data(57) = bt_str2double(RRStartTime);
            Data(58) = bt_str2double(RRStopTime);
            Data(59) = RRStraightCorr;
            Data(60) = RRCorr;
            Data(61) = RRLag;
        end


        biomarkXL = 'bt_gui_BioMARK.xls';

        if ~exist([cd, '\outputFiles\', biomarkXL])
            bt_xlswrite('BioMARK data', '', '', [cd, '\outputFiles\', biomarkXL]);
        end

        [nextline nextcolumn] = findNextXLLine([cd, '\outputFiles\', biomarkXL], 'Sheet1');


        if X == 1;
            xlswrite2([cd, '\outputFiles\', biomarkXL], Labels, 'Sheet1', 'A1');
        end

        xlswrite2([cd, '\outputFiles\', biomarkXL], {QuietFileName}, 'Sheet1', ['A', num2str(nextline)]);
        xlswrite2([cd, '\outputFiles\', biomarkXL], Data, 'Sheet1', ['B', num2str(nextline)]);


    end
end
