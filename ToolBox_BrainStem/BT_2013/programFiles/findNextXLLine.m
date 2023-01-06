function [nextline nextcolumn] = findNextXLLine(filename, sheetname)
% Description:  returns the next free line  and next free column in the specified Excel file.
%               When not specified, 'sheetname' defaults to 'Sheet1'. User will be
%               alerted when specified sheet does not exist.

% Author: eeskoe May 12, 2006
% Warning:  This file has not been exhaustively debugged. 


% call Excel
Excel = actxserver('Excel.Application');
% By setting visible to zero, the action will go on behind the scenes.
set(Excel,'Visible',0);

% This turns off Excel's alert windows, so that the user isn't asked
% whether or not he wants to overwrite existing files.  The reason for
% turning this off is that the user is prompted at the  beginning of
% ABR Batch 2 if XXX-FFT.xls exists, and is given a choice about continuing.
% If the user selects, "Yes, Continue", files will get overwritten without being prompted by Excel.
% By alterting the user early on this saves having to abort processing midway
% through.  
Excel.DisplayAlerts = 0;
Workbooks = Excel.Workbooks;             % Call up the Excel Workbook obect
% create new workbook.


ExcelWorkbook = Excel.workbooks.Open([filename]);
ActiveWorkbook = Excel.ActiveWorkbook;

if nargin<2
    sheetname ='Sheet1';
end

Sheets = Excel.ActiveWorkBook.Sheets;
try 
    ws = get(Sheets, 'Item', sheetname);
    
catch
    disp('Error:  Worksheet not found');
    release(ActiveWorkbook);
    invoke(Excel, 'Quit');
    delete(Excel)
    return
end

invoke(ws, 'activate');
DataRange = Excel.ActiveSheet.UsedRange;
rawData = DataRange.Value;
nextline = size(rawData, 1) +1;
nextcolumn = size(rawData,2)+1;

release(ActiveWorkbook);
invoke(Excel, 'Quit');
delete(Excel)

