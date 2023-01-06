function r = fisherinv(zprime)
% Description: performs inverse Fisher's transform
% Use this function to convert a z' score
%   to a correlation coefficient.
% Usage:  fisherinv(zprime)
r = (exp(zprime./0.5)-1)./(exp(zprime./0.5)+1);

