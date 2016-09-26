
%% Choos and load examples

image_name = 'DSC06083'; % 'DSC05748';

gt = imread(sprintf('data/sample_images/%s_mask.PNG', image_name));
pixmap = imread(sprintf('data/sample_images/%s.TIF', image_name));

cam = load('data/camera_models/Sony_A57_tiff.mat');

%% Detect forgeries

response_ss = detectForgeryPRNUCentral(pixmap, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8, 'image_padding', true));

response_sg = detectForgeryPRNUCentral(pixmap, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8, 'segmentwise_correlation', true, ...
    'image_padding', true));

response_aw = detectForgeryPRNUAdaptiveWnd(pixmap, cam.camera_model, ...
    struct('verbose', true, 'stride', 8, 'image_padding', true));

%% Post-process maps

decision_ss = response_ss.map_prp > 0.5;
decision_ss = mapCleanup(decision_ss, 64);
decision_ss = imdilate(decision_ss, strel('disk', 2));

decision_sg = fuseCRF(response_sg.map_prp, pixmap, 0.5, [-0.5 1 3 25 0]);
decision_aw = fuseCRF(response_aw.map_prp, pixmap, 0.5, [-1 1 4 25 0]);

%% Display results

clf;
subplot(2,4,1);
imsc(pixmap, 'tampered image');

subplot(2,4,5);
imsc(gt, 'ground truth');

subplot(2,4,2);
imsc(response_ss.map_prp, 'tamp. prob. (single-scale)');

subplot(2,4,6);
imsc(colorCodeMap(decision_ss, gt), 'decision (single-scale)');

subplot(2,4,3);
imsc(response_sg.map_prp, 'tamp. prob. (segmentation-guided)');

subplot(2,4,7);
imsc(colorCodeMap(decision_sg, gt), 'decision (segmentation-guided)');

subplot(2,4,4);
imsc(response_aw.map_prp, 'tamp. prob. (adaptive-window)');

subplot(2,4,8);
imsc(colorCodeMap(decision_aw, gt), 'decision (adaptive-window)');