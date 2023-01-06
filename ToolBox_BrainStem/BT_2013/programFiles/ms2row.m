function RowNum = ms2row(AvgStruct,Ms)
% calculates row number of a given latency for a Neuroscan .avg file.
% example, assuming AVG is a structure array resulting from a openavg
% command,
%
% rownumber = ms2row(AVG,23.5)
% returns the row number corresponding to 23.5 ms.
            
PtsPerMs = (AvgStruct.pnts-1)/(AvgStruct.xmax-AvgStruct.xmin);
Offset = PtsPerMs.*AvgStruct.xmin-1;
RowNum=round(Ms.*PtsPerMs-Offset);