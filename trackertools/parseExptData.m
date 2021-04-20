%
% Reads in text file containing experiment data.
% SAMPLE TEXT FILE:
%     ----------------------
%     Genotype,	N2
%     Stimulus,	1,4,IAA,1e-7
%         ,	7,7.5,IAA,1e-7
%         ,	8,8.5,IAA,1e-7
%         ,	9,9.5,IAA,1e-7
%         ,	11,12,IAA,1e-7
%         ,	13,14,IAA,1e-7
%     NumCycles,	4
%     CycleLen,	15
%     ----------------------
%   all times in minutes; Stimulus definition lists on time,off time,stimulus type,concentration
%
% USAGE:
%   data = parseExptData(filename)
%
%   filename: text filename
%   data: structure containing experiment data

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 30-Mar-2011 13:23:10 
%---------------------------- 

function data = parseExptData(filename)

try
    fid = fopen(filename);
    rawdata = textscan(fid,'%s%s%s%s%s%s%s%s%s%s','delimiter',',');
    fclose(fid);
catch
    disp(['Error reading file:',filename])
end

datacell = [];
for col = 1:length(rawdata)
    datacell = safecat(2,datacell,rawdata{col},{''});
end

[pathname,basename] = fileparts(filename);

timerange = 0;
data = struct('Name',filename,...
              'Path',pathname);
              
for row = 1:size(datacell,1)
    trimcell = [];
    for col = 2:size(datacell,2)
        value = datacell(row,col);
        if length(char(value)) > 0 
            trimcell = [trimcell, value];
        end
        if length(trimcell) == 0 trimcell = {}; end
    end
        
    varname = char(datacell(row,1));
    if length(varname) > 0 && ~strcmp(varname,'Stimulus')
        data = setfield(data,varname,trimcell);
    else
        timerange = timerange + 1;
        data.Stimulus(timerange).Time = str2num(char(trimcell(1:2)))';
        data.Stimulus(timerange).Odor = char(trimcell(3));
        data.Stimulus(timerange).Pattern = str2num(char(trimcell(4:end)))';
    end
end

