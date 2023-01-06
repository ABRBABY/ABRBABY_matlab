function RowNum = ms2row2(AvgStruct,Ms)
% calculates row number of a given latency for a Neuroscan .avg file.
% example, assuming AVG is a structure array resulting from a eeg_load
% command,
%
% rownumber = ms2row2(AVG,23.5)
Ms = Ms/1000;

PtsPerMs = (AvgStruct.header.pnts-1)/(AvgStruct.header.xmax-AvgStruct.header.xmin);
Offset = PtsPerMs.*AvgStruct.header.xmin-1;
RowNum=round(Ms.*PtsPerMs-Offset);