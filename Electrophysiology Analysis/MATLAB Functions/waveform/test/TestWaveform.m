function TestWaveform(whichTest)

if nargin < 1
  whichTest = 'all';
end

thisFile = mfilename('fullpath');
dataDir = [thisFile(1:(end-12)), 'Datasets/'];

if strcmp(whichTest, 'Jon') || strcmp(whichTest, 'all')
  testWaveform_Jon(dataDir)
end

if strcmp(whichTest, 'Tilman') || strcmp(whichTest, 'all')
  testWaveform_Tilman(dataDir)
end

if strcmp(whichTest, 'Maria') || strcmp(whichTest, 'all')
  testWaveform_Maria(dataDir)
end

if strcmp(whichTest, 'Lamont') || strcmp(whichTest, 'all')
  testWaveform_Lamont(dataDir)
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveform_Jon(dataDir)
testDir = [dataDir, 'Jon'];

dirList = dir(testDir);
for n = 1:length(dirList)
  name_n = dirList(n).name;
  if strncmp(name_n, '.', 1)
    continue;
  end
  disp(name_n)
  clear t v;
  load([testDir, '/', name_n]);
  
  An = AnalyzeWaveform(t, v, 'plotSubject', name_n, 'bracketWidth', 4.0); %#ok<NASGU>
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveform_Tilman(dataDir)
testDir = [dataDir, 'Tilman'];

% Tilman 792_10: 0 - 122
%   Tons of missed spikes ?
%   Many bursts detected, none present
% Tilman 792_10: 122 - 300
%   Bursts detected (none real)
% Tilman Bursting:
%   Missing spikes on end of burst (2.6s, 3.6 s, etc)  --> fixed
% Tilman STA data 01: OK
% Tilman STA data 02, OK
%   Missing spike at t = 66.8, 101.6, ...

dirList = dir(testDir);

v = [];
for n = 1:length(dirList)
  name_n = dirList(n).name;
  if strncmp(name_n, '.', 1)
    continue;
  end
  disp(name_n)
  
  title = '';
  load([testDir, '/', name_n]);
  if isempty(title)
    title = name_n;
  end
  
  if length(v) * dt / 1000 > 45
    fprintf('  -> skipping long file: %s\n', name_n)
    continue
  end
  
  if strcmp(name_n, 'Tilman_792_10.mat')
    % divide up trace into first part (extreme noise) and second
    ind = round(122000 / dt);
    % skip this one
    %beginTitle = [title, ': 0 - 122'];
    %An = AnalyzeWaveform(dt, v(1:ind), 'plotSubject', beginTitle); %#ok<NASGU>
    endTitle = [title, ': 122 - 300'];
    An = AnalyzeWaveform(dt, v((ind+1):end), 'plotSubject', endTitle); %#ok<NASGU>
  else
    An = AnalyzeWaveform(dt, v, 'plotSubject', title, 'debugPlots', true); %#ok<NASGU>
    %An = GetSpikes(dt, v, 'plotSubject', title, 'debugPlots', false); %#ok<NASGU>
  end
  
  clear v dt title
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveform_Maria(dataDir)

testDir = [dataDir, 'Maria'];

dirList = dir(testDir);

for n = 1:length(dirList)
  name_n = dirList(n).name;
  if strncmp(name_n, '.', 1)
    continue;
  end
  disp(name_n)
  
  abfStruct = LoadAbf([testDir, '/', name_n]);
  
  if any(strcmp(fieldnames(abfStruct.units), 'MCVm'))
    abfStruct.units.MCVm = 'pA';
  end
  
  GetFICurve(abfStruct, 'plotSubject', name_n);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveform_Lamont(dataDir)

testDir = [dataDir, 'Lamont'];
extracellularFile = [testDir, '/', '753_079_0036_pdn_extra.abf'];

[t, pdn] = GetExtracellular(extracellularFile, 'pdn');

extra = AnalyzeExtracellular(t, pdn(1,:), 'plotSubject', 'pdn', ...
                             'debugPlots', true);
extra.burst
extra.burst.durations

return