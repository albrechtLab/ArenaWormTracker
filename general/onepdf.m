%
% Create a single PDF document from multiple individual PDFs. Uses
% GPL-licensed PDF ToolKit (pdftk):
%       www.pdflabs.com
% 
% Requires pdftk.exe in any MATLAB path.
%       
% USAGE:
%   onepdf(pdfnames,outputname,deletefiles,Overwrite)
%
%   pdfnames: input filenames (can use wildcard *)
%   outputname: name of nultipage output PDF
%   deletefiles: removes original PDF files after combining
%   Overwrite: allow output file to be overwritten if it already exists

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 14-Apr-2011 10:22:49 
%---------------------------- 

function onepdf(pdfnames,outputname,deletefiles,Overwrite)

if ~exist('pdfnames') || isempty(pdfnames) 
    [f,p] = uigetfile('MultiSelect','on');
    pdfnames = []; 
    for i=1:length(f); 
        pdfnames = cat(2,pdfnames,['"',p,char(f(i)),'" ']); 
    end  
end
if ~exist('outputname') || isempty(outputname) outputname = [p,'combinedPDF.pdf']; end
if ~exist('deletefiles') || isempty(deletefiles) deletefiles = true; end
if ~exist('Overwrite') || isempty(Overwrite) Overwrite = true; end

if isempty(strfind(outputname,'.pdf')) outputname = [outputname,'.pdf']; end

D = dir(outputname); OutputExist = size(D,1)>0;
if OutputExist
    if ~Overwrite
        error(['Output file ',PDFname,' exists.']);
        return;
    end
end

pdfpath = which('pdftk.exe');
if isempty(pdfpath)
    error('Can''t find pdftk.exe to create multipage pdf.  Place this file (http://www.accesspdf.com/pdftk/) in the matlab path.');
    return;
end

%Prevent filename issues
% CleanPDFname = strrepl(PDFname,'Documents and Settings','DOCUME~1');
% CleanPDFname = strrepl(CleanPDFname,'All Users','ALLUSE~1');
% CleanPDFname = strrepl(CleanPDFname,' ','_');
% CleanPDFlist = strrepl(PDFlist,'Documents and Settings','DOCUME~1');
% CleanPDFlist = strrepl(CleanPDFlist,'All Users','ALLUSE~1');

% if (length(strfind(CleanPDFname, ' ')) > 0 | length(strfind(CleanPDFlist, ' ')) > 0)
%     disp(['Filename Error: Check for spaces in file path: ',CleanPDFname]);
% end

% com = sprintf('!"%s" %s cat output %s',pdfpath,CleanPDFlist,CleanPDFname);
com = sprintf('!"%s" %s cat output %s',pdfpath,pdfnames,outputname);

currentpath = pwd;
if currentpath(1) == '\'
    p = path; lastpath = p(max(find(p == ';'))+1:length(p));
    cd(lastpath);
end

if (OutputExist & Overwrite) delete(outputname); end
eval(com);

disp(['Converted to PDF file: ',outputname]);
cd(currentpath);

% % delete " marks
% outputname(strfind(outputname,'"'))='';
% pdfnames(strfind(pdfnames,'"'))='';

D = dir(outputname); OutputExist = size(D,1)>0;
if OutputExist & deletefiles
    disp('Deleting single pdf pages');
    D = dir(pdfnames); D = {D.name}';                   % get pdf file list
    filestodelete = D(find(~strcmp(D,outputname)));     % exclude the combined pdf
    oldstate = recycle; recycle on;                     % turn on recycle bin, just in case
    for i=1:length(filestodelete); 
        delete(char(filestodelete(i)));                 % delete pdfs
    end; 
    recycle(oldstate);                                  % reset recycle bin state
end