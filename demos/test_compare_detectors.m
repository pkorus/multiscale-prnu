%% Compare Detectors
%
% Runs different variants of tampering localization and compares the
% resulting correlation fields and tampering probability maps.

%% Load an example image

image_name = 'DSC06083'; 

gt = imread(sprintf('data/sample_images/%s_mask.PNG', image_name));
pixmap = imread(sprintf('data/sample_images/%s.TIF', image_name));

cam = load('data/camera_models/Sony_A57_tiff.mat');

%% Detect forgeries using all available algorithms

responses = struct();
timing = struct();

start_time = clock();
responses.ss = detectForgeryPRNUCentral(pixmap, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8, 'image_padding', true));
timing.ss = etime(clock(), start_time);

start_time = clock();
responses.sg = detectForgeryPRNUCentral(pixmap, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8, 'segmentwise_correlation', true, ...
    'image_padding', true));
timing.sg = etime(clock(), start_time);

start_time = clock();
responses.aw = detectForgeryPRNUAdaptiveWnd(pixmap, cam.camera_model, ...
    struct('verbose', true, 'stride', 8, 'image_padding', true));
timing.aw = etime(clock(), start_time);

start_time = clock();
responses.std = detectForgeryPRNUMultiscale(pixmap, cam.camera_model, ...
    129, 8, 'standard', true);
timing.std = etime(clock(), start_time);

start_time = clock();
responses.fwa = detectForgeryPRNUMultiscale(pixmap, cam.camera_model, ...
    128, 8, 'blockwise', true);
timing.fwa = etime(clock(), start_time);

start_time = clock();
responses.fc = detectForgeryPRNUMultiscale(pixmap, cam.camera_model, ...
    129, 8, 'central', true);
timing.fc = etime(clock(), start_time);

start_time = clock();
responses.fb = detectForgeryPRNUMultiscale(pixmap, cam.camera_model, ...
    129, 8, 'boxcar', true);
timing.fb = etime(clock(), start_time);

start_time = clock();
responses.fg = detectForgeryPRNUMultiscale(pixmap, cam.camera_model, ...
    129, 8, 'guided', true);
timing.fg = etime(clock(), start_time);

%% Display correlation fields and probability maps 

detectors = {'fwa', 'ss', 'fc', 'fb', 'fg', 'aw', 'sg'};

labels = {'Standard FWA detector', ...
    'Standard CPA detector', 'Filter-based CPA detector', ...
    'Filter-based CPA detector (boxcar post.)', 'Guided filter-based CPA detector', ...
    'Adaptive-window CPA detector', 'Segmentation-guided CPA detector'};

% Print timing statistics:
for i = 1:numel(detectors)
    runtime = timing.(detectors{i});
    runtime_vis = repmat('#', [1, round(runtime / 5)]);
    fprintf('%50s : %5.1f sec : %s\n', labels{i}, runtime, runtime_vis);
end

figure('Name', 'Correlation Fields');
set(0,'DefaultFigureWindowStyle','docked');
subplot(4,2,1); imsc(pixmap, 'Tampered image');
for i = 1:numel(detectors)
    subplot(4,2,i+1); imsc(responses.(detectors{i}).map_cor, labels{i});
end

figure('Name', 'Tampering Probability Maps');
set(0,'DefaultFigureWindowStyle','docked');
subplot(4,2,1); imsc(pixmap, 'Tampered image');
for i = 1:numel(detectors)
    subplot(4,2,i+1); imsc(responses.(detectors{i}).map_prp, labels{i});
end
