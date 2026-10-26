function P = projectPaths()
%PROJECTPATHS  Absolute paths of the repository folders, so scripts work from any folder.
P.root    = fileparts(fileparts(mfilename('fullpath')));
P.data    = fullfile(P.root, 'data', 'BBBC039');
P.results = fullfile(P.root, 'results', 'matlab');
P.figures = fullfile(P.root, 'figures');
P.pyResults = fullfile(P.root, 'results');      % output of the Python reference
if ~exist(P.results, 'dir'), mkdir(P.results); end
if ~exist(P.figures, 'dir'), mkdir(P.figures); end
end
