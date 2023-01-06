function varargout = bt_ptgui(varargin)
% The brainstem toolbox  is free software: you can redistribute it and/or modify
    %~ it under the terms of the GNU General Public License as published by
    %~ the Free Software Foundation.
%
    %~ The brainstem toolbox is distributed in the hope that it will be useful,
    %~ but WITHOUT ANY WARRANTY; without even the implied warranty of
    %~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    %~ See the GNU General Public License for more details <http://www.gnu.org/licenses/>/
%	Copyright 2007-2008, Trent Nicol and Erika Skoe.  
%
% Dependencies:
% brainstem toolbox m-files: bt_fft2, bt_fftsc, 
%          bt_peaks2, bt_qncorr, bt_rms, bt_xlswrite
% other m-files: closestrc, ms2row, nancorrcoef, openavg, xcorrelation (at
%          least)
% fig-files: bt_ptgui, bt_ptgui_out
% avg-files: 
%
% Sept-Dec, 2007, Trent Nicol.  Based on "Brainstem Toolbox" by Erika Skoe and
% Trent Nicol, and much underlying code by Erika Skoe.
addpath([cd, '\programFiles'])

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bt_ptgui_OpeningFcn, ...
                   'gui_OutputFcn',  @bt_ptgui_OutputFcn, ...
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
function bt_ptgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bt_bui (see VARARGIN)
% guifig = openfig([cd, '\programFiles\bt_ptgui.fig'],'reuse');   % the same window is used each time within 1 matlab session.   
%                                          % this allows path information to be retained from subject to subject
%                                          set(guifig, 'visible', 'on'); 
% Choose default command line output for bt_bui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = bt_ptgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ----------------------------END OF AUTO-GENERATED CODE-------------------

%% set defaults ------------------------------------------------------
warning 'off' 'all'


temp = get(handles.block, 'String');
if isempty(char(temp));
  set(handles.block, 'String', 40);
end

temp = get(handles.step, 'String');
if isempty(char(temp));
 set(handles.step, 'String', 1);
end

temp = get(handles.startSTIM, 'String');
if isempty(char(temp));
set(handles.startSTIM, 'String', 0);
end

temp = get(handles.endSTIM, 'String');
if isempty(char(temp));
set(handles.endSTIM, 'String', 150);
end

temp = get(handles.neural_lag, 'String');
if isempty(char(temp));
set(handles.neural_lag, 'String', 10);
end

temp = get(handles.channel, 'String');
if isempty(char(temp));

set(handles.channel, 'String', 2);
end

temp = get(handles.minFrequency, 'String');
if isempty(char(temp));
set(handles.minFrequency, 'String', 0);
end

temp = get(handles.maxFrequency, 'String');
if isempty(char(temp));

set(handles.maxFrequency, 'String', 120);
end

temp = get(handles.minFrequency, 'String');
if isempty(char(temp));


set(handles.minFrequency_stim, 'String', 80);
end


temp = get(handles.minFrequency_stim, 'String');
if isempty(char(temp));

set(handles.minFrequency_stim, 'String', 80);
end

temp = get(handles.maxFrequency_stim, 'String');
if isempty(char(temp));

set(handles.maxFrequency_stim, 'String', 120);
end



temp = get(handles.ExcelYN, 'String');
if isempty(char(temp));


set(handles.ExcelYN, 'Value', 0);
end




 
%% get the browse-button info. -------------------------------------------
function response_name_Callback(hObject, eventdata, handles)
    global response_name response_name_path
    [response_name response_name_path] = uigetfile('*.avg','Select response to be analyzed');
    % display the file selected to the right of the browse button.
    set(handles.response_name_text, 'String', [ response_name]); 

function stimulus_name_Callback(hObject, eventdata, handles)
    global stimulus_name stimulus_name_path
    [stimulus_name stimulus_name_path] = uigetfile('*.avg','Select stimulus to be analyzed');
    % display the file selected to the right of the browse button.
    set(handles.stimulus_name_text, 'String', [ stimulus_name]); 
    
    
%% Execute (upon clicking of OK button) ----------------------------------
function okButton_Callback(hObject, eventdata, handles)

hwaitbar = msgbox('Please wait...');


global guifig 
global Identifier ExcelYN response_name response_name_path
global stimulus_name stimulus_name_path


%% Get inputed values. -------------------------------------------
block = str2double(get(handles.block, 'String'));
step = str2double(get(handles.step, 'String'));
startSTIM = str2double(get(handles.startSTIM, 'String'));
endSTIM = str2double(get(handles.endSTIM, 'String'));
neural_lag = str2double(get(handles.neural_lag, 'String'));
channel = str2double(get(handles.channel, 'String'));
minFrequency = str2double(get(handles.minFrequency, 'String'));
maxFrequency = str2double(get(handles.maxFrequency, 'String'));
minFrequency_stim = str2double(get(handles.minFrequency_stim, 'String'));
maxFrequency_stim = str2double(get(handles.maxFrequency_stim, 'String'));
ExcelYN = get(handles.ExcelYN, 'Value');
Identifier = get(handles.Identifier, 'String');

% set some more defaults and compute some more values from inputs
stim_channel = 1; % i.e. if stereo file, will only work on 1st channel
startRESP = startSTIM + neural_lag;

                                 % this allows path information to be retained from subject to subject
set(guifig, 'visible', 'off'); 
%% Run pitchtrack functions ------------------------------------------------------

% pitchtrack response
[time autocorr lag FFT_resp freqaxis prestimFFT totalblocks]= pitchtrack([response_name_path response_name], block, step, startRESP, channel, 0);

% pitchtrack stimulus (make conditional)
[time_stim autocorr_stim lag_stim FFT_stim null freqaxis_stim totalblocks_stim]=pitchtrack([stimulus_name_path stimulus_name], block, step, startSTIM, stim_channel,0);

time_stim = time_stim + neural_lag;  %for plotting purpose shift stimulus forward in time.

[v stopPT]=closestrc(time_stim, endSTIM+(block./2));

autocorr_stim(:, stopPT+1:end)=[];
time_stim(stopPT+1:end)=[];
FFT_stim(:, stopPT+1:end)=[];
totalblocks_stim = stopPT;



%Extract F0 from stimulus using Autocorrelation Method 
freqAC_vector = 1000./lag_stim;
[s LagStart_stim]=closestrc(freqAC_vector,maxFrequency_stim);  %s is dummy variable
[s LagStop_stim] =closestrc(freqAC_vector,minFrequency_stim); 
[R_stim index] = max(autocorr_stim(LagStart_stim:LagStop_stim, :));
FreqAC_stim=freqAC_vector(LagStart_stim+index-1);
clear index s

%Extract F0 from response using Autocorrelation Method 
freqAC_vector = 1000./lag;
[s LagStart]=closestrc(freqAC_vector,maxFrequency);  %s is dummy variable
[s LagStop] =closestrc(freqAC_vector,minFrequency); 
[R index] = max(autocorr(LagStart:LagStop, :));
FreqAC=freqAC_vector(LagStart+index-1);
LAG = 1000./FreqAC;

%Extract F0 from stimulus using FFT Method 
[FFTAMP_stim index]= max(FFT_stim(minFrequency_stim+1:maxFrequency_stim+1, :));
FreqFFT_stim = minFrequency_stim+index-1;
FreqFFT_stim = FreqFFT_stim';

%Extract F0 from response using FFT Method 
[MaxFFTAMP index]= max(FFT_resp(minFrequency+1:maxFrequency+1, :));
FreqFFT = (minFrequency+index-1)'; %flip direction so that it matches FreqAC

%% Determine whether each extracted frequency was above the noise floor (NF)
% When the pitch-track is plotted a small gray dot will appear over the time
% ranges where the extracted frequency is below the noise floor The dot will be located above the plot and the location 
% is based on what was inputted for maxFrequency (maxFrequency+4).  The total (i.e. total_belowNF) that gets exported is based 
% only on the total number of blocks in the stimulus and not the total
% number of blocks in the response.  In other words, it only includes the blocks that were used to calculate pitch error.

for x = 1:totalblocks;
    F0_AMPusingACfreqs(x,1) = FFT_resp(round(FreqAC(x))-1,x);  %Find amplitude using frequencies extracted using AC method
end

mean_F0_AMP = mean(F0_AMPusingACfreqs(1:totalblocks_stim,1));
F0_AMP_prestim = prestimFFT(round(FreqAC)-1);
PITCH_SNR = F0_AMPusingACfreqs./F0_AMP_prestim ;

%Determine whether the extracted frequency was also the spectral maximum
for x = 1:totalblocks;
    peaks{x}=localmax(FFT_resp(:,x));  %Find amplitude using frequencies extracted using AC method
    [closest_spectralpeak(x) row(x) column(x)]=closestrc(peaks{x}, FreqAC(x));
    if  FreqFFT(x)== (closest_spectralpeak(x)-1)
       notspectralMax(x) = 0;
    else
       notspectralMax(x) = 1;
    end
end
   
% When the pitch-track is plotted a small gray dot will appear over the time
% ranges where the extracted frequency is not at the spectral max. The dot will be located above the plot and the location 
% is based on what was inputted for maxFrequency.  The total (i.e. total_notatSpectralMax) that gets exported is based 
% only on the total number of blocks in the stimulus and not the total number of blocks in the response.  In other words, \
% it only includes the blocks that were used to calculate pitch error.

total_notatSpectralMax = sum(notspectralMax(1:totalblocks_stim));
plot_notatspectralMax(notspectralMax==1)=maxFrequency+5; 
plot_notatspectralMax(notspectralMax==0)=NaN;

plot_belowNF(PITCH_SNR<1)=maxFrequency+4;
plot_belowNF(PITCH_SNR>=1)=NaN;

% Calculate Final Measures:
Method = get(handles.Method, 'Value');

switch Method
    case 1 %Autocorrelation
        PITCH_ERROR = mean(abs(FreqAC(1:totalblocks_stim)-FreqAC_stim));  %Measured in Hz.
        
    case 2 %FFT
        PITCH_ERROR = mean(abs(FreqFFT(1:totalblocks_stim)-FreqFFT_stim));  %Measured in Hz.
       
end

%on the off chance that one or more of the Rs is exactly 1, we must set
%these Rs to 0.999999 to get a valid number for fisher (i.e. not inf)If you are are reading this you are as annoyed as we are.
PITCH_STRENGTH = mean(R(1:totalblocks_stim));
Rtemp=R;
Rtemp(Rtemp==1)=0.999999;
PITCH_STRENGTH2 = fisherinv(mean(fisher(Rtemp(1:totalblocks_stim))));
CORR = corrcoef(FreqAC(1:totalblocks_stim), FreqAC_stim);  %correlation between stimulus and response f0 contour
PITCH_SRCORR =CORR(1,2); %the first number is always 1, need to take second 
total_belowNF =  sum(PITCH_SNR(1:totalblocks_stim)<1);
total_notatSpectralMax = sum(notspectralMax(1:totalblocks_stim));


close(hwaitbar) 
%% Populate output figure

% open and prep output figure 
OutputFig = openfig([cd, '\programFiles\', 'bt_ptgui_out.fig']);
orient landscape
handles = guihandles(OutputFig);
guidata(OutputFig, handles);

% fill in numbers
set(handles.Identifier, 'String', Identifier);
set(handles.StimFile, 'String', stimulus_name);
set(handles.ResponseFile, 'String', response_name);
set(handles.PitchError, 'String', num2str(PITCH_ERROR,'%3.2f'));
set(handles.PitchStrength, 'String', num2str(PITCH_STRENGTH2,'%3.2f'));
set(handles.Corr, 'String', num2str(PITCH_SRCORR,'%3.2f'));
set(handles.BelowNF, 'String', num2str(total_belowNF,'%3.0f'));
set(handles.NotSpectMax, 'String', num2str(total_notatSpectralMax ,'%3.0f'));

%plot pitch track
h = findobj(OutputFig, 'tag', 'TrackPlot');
colormap(h, 'hot')

switch Method
    case 1 %Autocorrelation
        plot(h, time, FreqAC, 's', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6);
        hold(h,'on')
        plot(h, time_stim, FreqAC_stim, 'k', 'LineWidth', 2);
        
    case 2 %FFT
        plot(h, time, FreqFFT, 's', 'color',  [1 0.7  0], 'MarkerFaceColor', 'y',  'MarkerSize', 6);
        hold(h,'on')
        plot(h, time_stim, FreqFFT_stim, 'k', 'LineWidth', 2);
       
end

xlabel(h, 'Time (ms)', 'FontSize',12);
ylabel(h, 'Frequency (Hz)', 'FontSize',12);
set(h, 'FontSize', 10);
ylim(h, [min(FreqAC)-20 max(FreqAC)+20])
% hold(h,'on')
plot(h, time, plot_belowNF, '*', 'color', 'b' , 'MarkerSize', 2 );
plot(h, time, plot_notatspectralMax, '*', 'color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 2);
xlim(h, [time(1) time(end)]);
title(h, 'Pitch Track', 'FontSize', 14);
hleg = legend('Response', 'Stimulus', 'below Noise Floor', 'Not Spectral Max', 'Location', 'SouthWest');
set(hleg, 'FontSize', 6);
set(h, 'Box', 'off');
xlim(h, [time(1) time(totalblocks_stim)]);

%plot autocorrelogram
h = findobj(OutputFig, 'tag', 'AutoPlot');
imagesc(time, lag, autocorr, 'parent', h);
ylabel(h, 'Lag (ms)', 'FontSize', 12);
set(h, 'FontSize', 10);
hold(h,'on')
caxis(h, [-1 1]);
switch Method
    case 1 %Autocorrelation
       plot(h, time, LAG, 'k', 'LineWidth', 2);
        
    case 2 %FFT
        plot(h, time, 1000./FreqFFT, 'k', 'LineWidth', 2);
       
end


xlim(h, [time(1) time(totalblocks_stim)]);
title(h, 'Running Autocorrelogram with Pitch Track', 'FontSize', 14);
c = findobj(OutputFig, 'tag', 'AutoColor');
colorbar('peer',c,[0.937 0.537 0.036 0.399])
caxis(c, [-1 1]);
set(c, 'FontSize', 8);
set(get(c,'ylabel'),'string','correlation coefficient (r)'); 
set(c, 'Box', 'off');
% must do this a second time to actually set axis.
colorbar('peer',c,[0.937 0.537 0.036 0.399])

%plot spectrogram
h = findobj(OutputFig, 'tag', 'FFTPlot');
FFT_cut = (FFT_resp(1:701,:));
% FFT_cut(FFT_cut<15)=NaN;  %only plot amplitudes above 15
imagesc(time(1:totalblocks_stim), freqaxis(1:701), FFT_cut(:, 1:totalblocks_stim), 'parent', h);
hold(h, 'on')
switch Method
    case 1 %Autocorrelation
       plot(h, time, FreqAC, 'k', 'LineWidth', 2);
        
    case 2 %FFT
        plot(h, time, FreqFFT, 'k', 'LineWidth', 2);
       
end

xlabel(h, 'Time (ms)', 'FontSize', 12);
ylabel(h, 'Frequency (Hz)', 'FontSize', 12);
set(h, 'FontSize', 10);
set(h, 'YDir', 'normal');
xlim(h, [time(1) time(totalblocks_stim)]);
title(h, 'Running FFT', 'FontSize', 14);
c = findobj(OutputFig, 'tag', 'FFTColor');
colorbar('peer',c,[0.937 0.07 0.036 0.399])
set(c, 'FontSize', 8);
set(get(c,'ylabel'),'string','microvolts'); 
set(c, 'Box', 'off');
% must do this a second time to actually set axis.
colorbar('peer',c,[0.937 0.07 0.036 0.399])

%% Save to Excel

if ExcelYN == 1
    % align output: identifiying info and labels in cell array, data in num array.
    ID{1} = Identifier;
    ID{2} = ['Response: ' response_name];
    ID{3} = ['Stimulus: ' stimulus_name];
    ID{4} = ['Number of time chunks: ' num2str(totalblocks_stim)];
    Labels = {...
        'Pitch Error','Pitch Strength','Correlation','# below Noise Floor','# Not Spectral Max'};
    Data(1) = PITCH_ERROR;
    Data(2) = PITCH_STRENGTH2;
    Data(3) = PITCH_SRCORR;
    Data(4) = total_belowNF;
    Data(5) = total_notatSpectralMax;
    
    if Identifier
        bt_xlswrite(Data,ID,Labels,[cd, '\outputFiles\', Identifier, '_pt'],'Sheet1');
    else
        bt_xlswrite(Data,ID,Labels,[cd, '\outputFiles\',response_name(1:end-4)], '_pt','Sheet1');
    end
end

%% print
