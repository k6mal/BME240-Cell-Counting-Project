function stems = readSplit(dataDir, splitName)
%READSPLIT  File stems for one official BBBC039 split ('training','validation','test').
txt   = fileread(fullfile(dataDir, 'metadata', [splitName '.txt']));
lines = strtrim(splitlines(txt));
lines = lines(~cellfun(@isempty, lines));
stems = erase(lines, '.png');
end
