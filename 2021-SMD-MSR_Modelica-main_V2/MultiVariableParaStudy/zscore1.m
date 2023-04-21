function [y,meanval,stdval] = zscore1(x, mvin,stdin)

%ZSCORE1    MCUV scale data.
%   This function mean center, unit variance (MCUV) scales data.
%
%   [XS XMean XSTD] = ZSCORE1(X) MCUV scales X, returning the scaled data 
%   XS, scaling mean XMean, and scaling standard deviation XSTD.
%
%   XS = ZSCORE1(X,XMean,XSTD) MCUV scales X with the scaling mean XMean 
%   and scaling standard deviation XSTD.  This implementation is used to 
%   scale new data according to traininig parameters. 
%
%   The dimension of X and XS is NumObs x NumVars, where NumObs is the 
%   number of observations and NumVars is the number of variables.  The 
%   dimension of XMean and XSTD is 1 x NumVar.
%
%   See also UNSCORE, SMODDATA, USMODDATA.

%   J. Wesley Hines
%   The University of Tennessee, Knoxville
%   Nuclear Engineering Department
%   Last Update:    11/5/2005
%
%   Copyright (c)

[nrows,ncols]=size(x);

if nargin == 1
   meanval = mean(x);  	% calculate mean values
else
   meanval = mvin;      % use previously calculated value
end

y = x - ones(nrows,1)*meanval; 	% subtract off mean 

if nargin == 1
   stdval = std(y);	% calculate the SD
else
   stdval = stdin;	% use previously calculated value
end

% normalize to unit variance
if stdval~=0;
    y = y ./ (ones(nrows,1)*stdval);  
end;

return;

%--------------------------------------------------------------------------