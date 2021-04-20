%
% Add orthogonal vectors, expanding each.
%
% USAGE:
%   output = addortho(vr, vc)
%
%   vr: row vector          e.g. [1 2 3 4]
%   vc: column vector       e.g. [2; 5]
%   output: summed matrix   e.g. [3 4 5 6;
%                                 6 7 8 9];

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 25-Sep-2009 11:02:56 
%---------------------------- 

function output = addortho(vr, vc)

sizes = [size(vr), size(vc)];

if size(vr,1) == 0 
    output = vc; 
elseif size(vc,1) == 0
    output = vr;
elseif sizes(1) == 1 && sizes(4) == 1
    v = vr; vr = vc; vc = v;
    output = repmat(vc,size(vr,1),1) + repmat(vr,1,size(vc,2));
elseif sizes(2) == 1 | sizes(3) == 1
    output = repmat(vc,size(vr,1),1) + repmat(vr,1,size(vc,2));
else
    output = [];
    error('need one column and one row vector');
end
