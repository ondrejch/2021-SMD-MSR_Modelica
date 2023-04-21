function [y] = unscore(x, meanval,stdval)

%UNSCORE    Un-scale MCUV scaled data.
%   This function un-scales mean center, unit variance scaled data.
%   
%   X = UNSCORE(XS,XMean,XSTD) un-scales XS with the scaling means XMean 
%   and standard deviations XSTD.
%
%   The dimension of X and XS is NumObs x NumVars, where NumObs is the 
%   number of observations and NumVars is the number of variables.  The 
%   dimension of XMean and XSTD is 1 x NumVar.
%
%   See also ZSCORE1, SMODDATA, USMODDATA.

%   J. Wesley Hines
%   The University of Tennessee, Knoxville
%   Nuclear Engineering Department
%   Last Update:    6/10/1998
%
%   Copyright (c)

[nrows,ncols]=size(x);

y = x .* (ones(nrows,1)*stdval);  % Un-normalize from unit variance.

y = y + ones(nrows,1)*meanval; 	 % Add back on the mean.

return;

%--------------------------------------------------------------------------