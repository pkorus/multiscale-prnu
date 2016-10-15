verbose = false;
extra_paths = strsplit(genpath('3rd-party/'),':');

paths = {
	'commons'
	'demos'
	'detectors'
	'training'
	'fusion'
	};

paths = [paths ; extra_paths';];

for i=1:length(paths)
    if numel(paths{i}) > 0 && exist(paths{i}, 'dir')
        % Skip directories without useful files
        if numel(dir([paths{i} '/*.m'])) + numel(dir([paths{i} '/*.' mexext()])) > 0
            if verbose; fprintf('Adding to path: %s\n', paths{i}); end;
            addpath(paths{i});
        else
            if verbose; fprintf('Skipping: %s (no M-files or MEX files found)\n', paths{i}); end
        end
    end
end

clearvars;

global scenario_path
scenario_path = './data-test-scenarios/';

if ~exist('cache', 'dir')
    mkdir('cache');
end