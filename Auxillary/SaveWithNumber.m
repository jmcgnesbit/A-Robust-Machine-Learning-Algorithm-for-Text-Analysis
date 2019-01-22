function SaveWithNumber(FileName, Data)
[fPath, fName, fExt] = fileparts(FileName);
if isempty(fExt)  % No '.mat' in FileName
  fExt     = '.mat';
  FileName = fullfile(fPath, [fName, fExt]);
end
if exist(FileName, 'file')
  % Get number of files:
  fDir     = dir(fullfile(fPath, [fName, '*', fExt]));
  fStr     = lower(sprintf('%s*', fDir.name));
  fNum     = sscanf(fStr, [fName, '%d', fExt, '*']);
  newNum   = max(fNum) + 1;
  FileName = fullfile(fPath, [fName, sprintf('%d', newNum), fExt]);
end
save(FileName, 'Data');