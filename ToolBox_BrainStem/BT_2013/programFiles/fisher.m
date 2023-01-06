function zprime = fisher(r)
% performs Fisher's transform
% Use this function to perform hypothesis testing on the correlation coefficient.
% Usage:  fisher(r)
zprime = log((1+r)./(1-r))./2;