%
% Displays a progress note in the MATLAB window that overwrites itself.
%
% USAGE:
%
%  <before loop>
%       i = -1;
%
%  <in loop>
%       i = showprog([ #1 #2 ], i);
%
%  <output>
%       #1 of #2

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 27-Apr-2009 11:28:20 
%---------------------------- 

function chnum = showprog(val,bschars)
if nargin < 2 bschars = 0; end

if isstr(val)
    str = val;
else
    if length(val)==1
        str = num2str(val);
    else
        str = [num2str(val(1)),' of ',num2str(val(2))];
    end
end
        
chnum = length(str);
bs = repmat('\b',1,bschars+2);

fprintf(1,[bs,'\r',str,'\r']);

