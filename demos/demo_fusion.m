%% Multi-scale Fusion Demo
% Demonstrates fusion strategies that combine multi-scale candidate
% tampering probability maps into a single map. 

%% Choose image for display

image_name = 'DSC06083';
camera_name = 'Sony_A57';

threshold = 0.4;

%% Load cached data

if exist('data-tifs-2016-maps', 'dir')
    data_path = 'data-tifs-2016-maps';
else
    data_path = 'data/sample_maps';
end

% Load the tampered image (down-sampled) and ground truth maps
pixmap = imread(sprintf('%s/%s/%s.TIF', data_path, camera_name, image_name));
ground_truth = imread(sprintf('%s/%s/%s_mask.PNG', data_path, camera_name, image_name));
ground_truth = imresize(ground_truth, [135 240]);

% Load multi-scale response maps
response_maps = loadMaps(camera_name, image_name, 'central');
response_maps = postprocessMaps(response_maps, @(x) imresize(x, 0.5, 'bilinear'));

% Multi-scale CRF-based fusion
fusion_labeling = fuseCRF(response_maps, pixmap, threshold, [-1 0.5 5.6 25 0.18]);

% Adaptive-window strategy (implemented as a fusion method)
fusion_aw = fuseAdaptWnd(response_maps, 0.1, 0.25);
% Here fuseCRF is used as an improved decision strategy capable of modeling
% adaptive neighborhood interactions.
localization_aw = fuseCRF(fusion_aw, pixmap, threshold, [-1 1 4 25]);

% Standard, 128 px window
standard = response_maps{5}.candidate > threshold;
standard = mapCleanup(standard, 64);
standard = imdilate(standard, strel('disk', 2));

%% Show results

subplot(3,2,1);
imsc(pixmap, 'tampered image')

subplot(3,2,2);
imsc(response_maps, 'candidate maps : %s')

subplot(3,2,3);
imsc(colorCodeMap(standard, ground_truth), 'Single-scale 128 px window')

subplot(3,2,4);
imsc(colorCodeMap(fusion_labeling, ground_truth), 'Multi-scale CRF fusion')

subplot(3,2,5);
imsc(fusion_aw, 'Adaptive-window (tamp. prob.)')

subplot(3,2,6);
imsc(colorCodeMap(localization_aw, ground_truth), 'Adaptive-window (decision)')
