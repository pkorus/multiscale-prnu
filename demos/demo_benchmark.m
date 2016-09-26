%% Multi-scale Fusion Benchmark
% 
% A simple benchmark that compares performance of multi-scale CRF fusion
% with a standard single-scale decision heuristic.
%
% Requires the multi-scale fusion dataset. Either use:
%
%   # ./configure data:maps
%
% or download manually from http://kt.agh.edu.pl/~korus 

%% Setup and list test images

data_path = 'data-tifs-2016-maps';

if ~exist(data_path, 'dir')    
    fprintf('Data-set not found (%s)! Use ./configure data:maps or see readme.md for more info.\n', data_path);
end

% Setup thresholds
thresholds = getThresholds(24);

% List all test images
cameras = {};
files = {};
sub_dirs = dir(data_path);
for sname = sub_dirs'
    if numel(sname.name) > 2 && sname.isdir
        candidates = dir(sprintf('%s/%s/*.TIF', data_path, sname.name));
        files = {files{:} candidates(:).name};
        cams = repmat({sname.name}, [numel(candidates) 1]);
        cameras = {cameras{:} cams{:}};
    end
end

fprintf('Found %d test files in %s directory\n', numel(files), data_path);

%% Run simple benchmark
results = struct();
results.fusion = cell(numel(files), numel(thresholds));
results.single = cell(numel(files), numel(thresholds));

last_progress = -1;
total = numel(files);

for file_id = 1:total
    
    camera_name = cameras{file_id};
    image_name = strrep(files{file_id}, '.TIF', '');

    progress = round(100*(file_id/total)/2);
    
    % Display progress bar
    if progress ~= last_progress
        displayProgress(sprintf('Image   %6d / %6d', file_id, total), progress, 50);
    end
    
    last_progress = progress;
    
    % Load the tampered image (down-sampled) and ground truth maps
    pixmap = imread(sprintf('%s/%s/%s.TIF', data_path, camera_name, image_name));
    ground_truth = imread(sprintf('%s/%s/%s_mask.PNG', data_path, camera_name, image_name));
    ground_truth = imresize(ground_truth, [135 240]);

    % Load multi-scale response maps
    response_maps = loadMaps(camera_name, image_name, 'central');
    response_maps = postprocessMaps(response_maps, @(x) imresize(x, 0.5, 'bilinear'));
    
    for t_id = 1:numel(thresholds)
        
        % CRF fusion
        map_fusion = fuseCRF(response_maps, pixmap, thresholds(t_id), [-1 0.55 5.6 25 0.18]);
        
        % Single-scale heuristic decision for the 128 px window
        map_single = response_maps{5}.candidate > thresholds(t_id);
        map_single = mapCleanup(map_single, 64);
        map_single = imdilate(map_single, strel('disk', 1, 8));
        
        % Scoring
        results.fusion{file_id, t_id} =  scoreLocalization(map_fusion, ground_truth);
        results.single{file_id, t_id} =  scoreLocalization(map_single, ground_truth);
    end
    
end

%% Show results
tpr_fusion = mean(extractField(results.fusion, 'tpr'));
tnr_fusion = mean(extractField(results.fusion, 'tnr'));

tpr_single = mean(extractField(results.single, 'tpr'));
tnr_single = mean(extractField(results.single, 'tnr'));

f1_fusion = mean(extractField(results.fusion, 'f1'));
f1_single = mean(extractField(results.single, 'f1'));

subplot(2,1,1);
plot(1 - tnr_fusion, tpr_fusion, 'k-o'); hold on;
plot(1 - tnr_single, tpr_single, 'b-s'); hold off;
legend('CRF fusion', 'single scale', 'Location', 'SouthEast');
xlim([0 0.2]);
ylim([0 1]);
xlabel('false positive rate');
ylabel('true positive rate');

subplot(2,1,2);
plot(thresholds, f1_fusion, 'k-o'); hold on;
plot(thresholds, f1_single, 'b-s'); hold off;
legend('CRF fusion', 'single scale');
xlim([0 1]);
ylim([0 0.7]);
xlabel('threshold');
ylabel('F_1 score');
