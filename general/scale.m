%
% Scales matrix to specified limits.
%
% USAGE:
%   scaled = scale(mat,loset,hiset,pct,limit);
%
%   mat: input matrix
%   loset, hiset: low and high limits for scaled output
%   pct: set output limits to input value at (pct) and (100-pct) percentiles
%        (default = 0 uses strict min and max of mat)
%   limit: cap any output values outside specified range (i.e. only applies
%        when using pct)

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 22-Apr-2010 11:12:54 
%---------------------------- 

function scaled = scale(mat,loset,hiset,pct,limit);

if~exist('limit') || isempty(limit) limit = 0; end
if ~exist('pct') || isempty(pct) pct = 0; end
if numel(pct) == 1 pct = [pct,pct]; end

if hiset < loset
    temp = loset;
    loset = hiset;
    hiset = temp;
end
rangeset = hiset - loset;

if pct == 0
    lo = min(reshape(mat,[],1));
    hi = max(reshape(mat,[],1));
else
    lo = prctile(reshape(mat,[],1),pct(1));
    hi = prctile(reshape(mat,[],1),100-pct(2));
end
range = hi - lo;

scaled = ((mat - lo) / range * rangeset) + loset;

if limit
    scaled = max(loset,min(hiset,scaled));
end

