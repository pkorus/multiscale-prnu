function validateDataset(camera_name, thresholds)
% validateDataset(camera_name, thresholds)
%
% Simple dataset + camera model validation script. Validation involves
% iteration of the tampering localization procedure over the pristine and 
% PRNU estimation sub-sets. Tampering localization maps are obtained by
% simple thresholding followed by standard heuristic cleaning (small CC
% removal). 
%
% The script shows statistics and examples of problematic images with
% excessive false positive rates. 
%
% Params:
%   camera_name   - self explanatory,
%   thresholds    - 2D vector with:
%                   a) threshold for tampering localization, def: 0.75
%                   b) threshold for acceptable FP rate per image, def: 0.1
%
% Validation results are cached in ./data-input/(camera_name)/validation.mat 
% for future reuse. Delete the file to recompute. 
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu-localization-toolbox
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: September 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------

    if nargin < 2 || isempty(thresholds)
        thresholds = [0.75, 0.1];
    end

    % Basic setup
    dirname = './data-input/';
    camera_model_filename = sprintf('./data/camera_models/%s.mat', camera_name);
    m = load(camera_model_filename);
    camera_model = m.camera_model;

    % Estimation of camera's PRNU -----------------------------------------
    
    fprintf('\n############################# VALIDATING DATA SET ################################\n\n')
    
    validation_filename = sprintf('data-input/%s/validation.mat', camera_name);

    if ~exist(validation_filename, 'file')

        validation_stats = struct();
        validation_stats.camera_name         = camera_name;

        % Pristine images -------------------------------------------------

        data_set_query = [dirname camera_name '/pristine/' '*.TIF'];
        filenames = dir(data_set_query);

        validation_stats.response_maps       = cell(numel(filenames),1);
        validation_stats.false_positive_rate = zeros(numel(filenames), 1);
        validation_stats.mean_response       = zeros(numel(filenames), 1);
        validation_stats.total_images        = numel(filenames);
        validation_stats.dataset_path        = data_set_query;
        validation_stats.camera_name         = camera_name;
        validation_stats.pristine_filenames  = filenames;

        last_progress = -1;
        total = size(filenames,1);

        for i = 1:size(filenames,1)

            counter = i;
            progress = round(100*(counter/total)/2);
            % Display progress bar
            if progress ~= last_progress
                displayProgress(sprintf('Pristine   % 5d', i), progress, 50);
            end
            last_progress = progress;

            image = imread([dirname camera_name '/pristine/' filenames(i).name]);
            response_map = detectForgeryPRNUFilter(image, camera_model, 129, struct('stride', 32, 'verbose', false));
            validation_stats.mean_response(i) = mean2(response_map.map_prp);
            validation_stats.response_maps{i} = response_map.map_prp;
        end

        % PRNU images -------------------------------------------------
        data_set_query = [dirname camera_name '/prnu/' '*.TIF'];
        filenames = dir(data_set_query);

        validation_stats.response_maps_prnu       = cell(numel(filenames),1);
        validation_stats.false_positive_rate_prnu = zeros(numel(filenames), 1);
        validation_stats.mean_response_prnu       = zeros(numel(filenames), 1);
        validation_stats.total_images_prnu        = numel(filenames);
        validation_stats.dataset_path_prnu        = data_set_query;
        validation_stats.prnu_filenames           = filenames;

        last_progress = -1;
        total = size(filenames,1);

        for i = 1:size(filenames,1)

            counter = i;
            progress = round(100*(counter/total)/2);
            % Display progress bar
            if progress ~= last_progress
                displayProgress(sprintf('PRNU       % 5d', i), progress, 50);
            end
            last_progress = progress;

            image = imread([dirname camera_name '/prnu/' filenames(i).name]);
            response_map = detectForgeryPRNUFilter(image, camera_model, 129, struct('stride', 32, 'verbose', false));
            validation_stats.mean_response_prnu(i) = mean2(response_map.map_prp);
            validation_stats.response_maps_prnu{i} = response_map.map_prp;
        end

        save(validation_filename, 'validation_stats', '-v7.3');
    else
        load(validation_filename);
    end

%%  Calculate stats

    for i = 1:validation_stats.total_images_prnu
        validation_stats.false_positive_rate_prnu(i) = mean2(mapCleanup(validation_stats.response_maps_prnu{i} > thresholds(1), 16*16));
    end

    for i = 1:validation_stats.total_images
        validation_stats.false_positive_rate(i) = mean2(mapCleanup(validation_stats.response_maps{i} > thresholds(1), 16*16));
    end
    
%%  Print and display stats

    fprintf('Camera name        : %s\n', validation_stats.camera_name);
    fprintf('Decision threshold : %.2f\n', thresholds(1));
    fprintf('Tested images      : %d (prnu)\n', validation_stats.total_images_prnu);
    fprintf('Tested images      : %d (pristine)\n', validation_stats.total_images);
    fprintf('\n');

    figure(1);
    subplot(2,1,1);
    hist(validation_stats.false_positive_rate_prnu, 26);
    stat = sum(validation_stats.false_positive_rate_prnu > thresholds(2));
    label = sprintf('%3d (%3.0f%%) invalid PRNU images     (false alarm rate > %.3f)\n', stat, 100*stat / validation_stats.total_images, thresholds(2));
    title(label);
    xlabel('false alarm rate');
    fprintf('%s', label);
    xlim([0 1.0]);
    
    subplot(2,1,2);
    hist(validation_stats.false_positive_rate, 26);
    stat = sum(validation_stats.false_positive_rate > thresholds(2));
    label = sprintf('%3d (%3.0f%%) invalid pristine images (false alarm rate > %.3f)\n', stat, 100*stat / validation_stats.total_images, thresholds(2));
    title(label);
    xlabel('false alarm rate');
    fprintf('%s', label);
    xlim([0 1.0]);

%%  List problematic images

    [ranking, indices] = sort(validation_stats.false_positive_rate, 'descend');
    fprintf('\n15 pristine images with the worst FP rates:\n');
    for i = 1:min(15, numel(ranking))
        fprintf(' %2d. %s --> false alarm rate = %.2f\n', i, validation_stats.pristine_filenames(indices(i)).name, ranking(i));
    end

    % show them...
    nRows = 135;
    nCols = 240;
    worst_images = zeros(nRows, nCols, 15);
    for i = 1:min(15, numel(ranking))
        img = imread([dirname camera_name '/pristine/' validation_stats.pristine_filenames(indices(i)).name]);
        worst_images(:,:,i) = rgb2gray(imresize(img, [nRows, nCols]));
    end
    
    figure(2); clf;
    imagesc(generateThumbnails(worst_images, [3 5])); cll;
    title('Pristine images with the worst FP rates');
    colormap gray;

%%  Show the invalid maps
    
    nRows = size(validation_stats.response_maps{1},1);
    nCols = size(validation_stats.response_maps{1},2);
    invalid_pristine_maps = zeros(nRows, nCols, stat);
%     selection = find(validation_stats.false_positive_rate > thresholds(2));
    
    if stat > 0
        for i = 1:stat
            invalid_pristine_maps(:,:,i) = validation_stats.response_maps{indices(i)};
        end

        figure(3); clf;
        imagesc(generateThumbnails(invalid_pristine_maps, stat)); cll;
        title('Response maps for pristine images with excessive FPs');
        colormap gray;
    end

%%  All images worse than 0.5

    indices = find(validation_stats.false_positive_rate > 0.5);
    
    fprintf('\n%d images with false alarm rate > 0.5:\n', numel(indices));
    
    if numel(indices) > 0    
        for i = 1:numel(indices)
            fprintf('%s\n', validation_stats.pristine_filenames(indices(i)).name);
        end

        % show them...
        nRows = 135;
        nCols = 240;
        worst_images = zeros(nRows, nCols, numel(indices));
        for i = 1:numel(indices)
            img = imread([dirname camera_name '/pristine/' validation_stats.pristine_filenames(indices(i)).name]);
            worst_images(:,:,i) = rgb2gray(imresize(img, [nRows, nCols]));
        end

        figure(4); clf;
        imagesc(generateThumbnails(worst_images, numel(indices))); cll;
        title(sprintf('%d images with FP rate > 0.5', numel(indices)));
        colormap gray;
    end
end