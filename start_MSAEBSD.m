function MSAEBSD=start_MSAEBSD(varargin)
% Load MSAEBSD

%set up the path
local_path = fileparts(mfilename('fullpath'));

fprintf('Loading MSAEBSD');

if MATLABverLessThan('9.0')
    error('MSAEBSD needs at least version 9 of Matlab')
end

%read MSA version from version file - adapted from MTEX
try
  fid = fopen('VERSION','r');
  MSAVersion = fgetl(fid);
  fclose(fid);
  fprintf([' ' MSAVersion '  ']);
catch
  MSAVersion = 'MSAEBSD';
end

%initalise folders as needed


%read the current folder
msaloc=strfind(local_path,'MSAofEBSD');

%check that w are in a subfolder
if isempty(msaloc)
    error('Path must be run first from a subfolder of ''MSAEBSD''');
end

%get the file path
MSA_FP=local_path(1:msaloc+8);

%build folders for things if this is a first run
if exist([MSA_FP '\testing'],'dir') ~= 7
    mkdir([MSA_FP '\testing']);
end
if exist([MSA_FP '\outputs'],'dir') ~= 7
    mkdir([MSA_FP '\outputs']);
end
if exist([MSA_FP '\decks'],'dir') ~= 7
    mkdir([MSA_FP '\decks']);
end

%build the paths

addpath([MSA_FP '\bin']);
addpath([MSA_FP '\testing']);
addpath([MSA_FP '\decks']);
addpath([MSA_FP '\outputs']);
addpath([MSA_FP]);

disp('MSAEBSD file paths loaded');

end

% check MATLAB version - borrowed from MTEX loader
% --------------------
function result = MATLABverLessThan(verstr)

MATLABver = ver('MATLAB');

toolboxParts = getParts(MATLABver(1).Version);
verParts = getParts(verstr);

result = (sign(toolboxParts - verParts) * [1; .1; .01]) < 0;

end

function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d')';
if length(parts) < 3
  parts(3) = 0; % zero-fills to 3 elements
end
end

function p()
if isempty(lasterr)
  fprintf('.');
end
end

function pathadd(path)
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(Folder, pathCell));
else
  onPath = any(strcmp(Folder, pathCell));
end
end
%eof