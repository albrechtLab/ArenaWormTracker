%
% Get user input, use for both numerical and string inputs.
%
% USAGE:
%   response = dagetnum(prompts, defaults, title)

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 29-Aug-2008 17:01:47 
%---------------------------- 

function response = dagetnum(prompts, defaults, title)

if nargin < 3 title = ''; end
if nargin < 2 defaults = zeros(size(prompts)); end
if ~iscell(prompts) prompts = {prompts}; end
    
if size(prompts,2) ~= size(defaults,2)
    error('dimensions don''t match');
    return
end

if ~iscell(defaults) defaults = cellstr(num2str(defaults'))'; end

answercell = inputdlg(prompts,title,1,defaults);

response = [];
for i = 1:size(answercell,1)
    response(i).num = str2num(char(answercell(i)));
    response(i).char = char(answercell(i));
end
