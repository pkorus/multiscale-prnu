function importDataset(camera_name, directory_path, image_counts, resolution)
% importDataset(camera_name, directory_path, image_counts, resolution)
%
% Imports images from filesystem to the testing environment. The process involves:
%  - creation of a corresponding data set in ./data-input directory
%  - division of input images into three sub-sets for:
%    a) PRNU estimation
%    b) predictor training
%    c) further testing
%  
% The function will randomly select two sub-sets for predictor training and further 
% testing. Then, it will select preferable remaining images (bright, low-texture) for
% PRNU estimation.
%  
% Images will be cropped to the middle fragment of a given size (1920 x 1080 by default).
%
% The structure of the created directories is as follows:
%
%   ./data-input/(camera_name)/{predictor,pristine,prnu,tampered-realistic}
%
% The function should recognize popular image formats like JPG, PNG, TIFF
% and PGM. All output images will be saved in the TIFF format.
%
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: September 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------

    if nargin < 3 || isempty(image_counts)
        prnu_images = 200;
        pred_images = 50;
        test_images = 400;
    else
        prnu_images = image_counts(1);
        pred_images = image_counts(2);
        test_images = image_counts(3);
    end
    
    if nargin < 4 || isempty(resolution)
        resolution = [1080 1920];
    end
    
    all_images = prnu_images + pred_images + test_images;
    
    if directory_path(end) ~= '/'
        directory_path = [directory_path '/'];
    end

    filenames = dir([directory_path '*']);    
    
    selection = false(1, numel(filenames));
    for fi = 1:numel(filenames)
        if regexp(filenames(fi).name, '\.(jpg|jpeg|tif|tiff|png|pgm|JPG|JPEG|TIF|TIFF|PNG|PGM)$')
            selection(fi) = true;
        end
    end
    filenames = filenames(selection);
    
    if numel(filenames) < all_images
        fprintf('Not enough (%d of %d needed) images found in %s\n', numel(filenames), all_images, directory_path);
        return;
    end
    
    % Remove everything that is not an image
    
    permutation = randperm(numel(filenames));    
    filenames = filenames(permutation);
    
    out_dirname = sprintf('./data-input/%s/', camera_name);
    
    fprintf('\n################################ DATASET IMPORT ####################################\n\n')
    fprintf('Importing %d random images out of %d found in %s\n', all_images, numel(filenames), directory_path);
    fprintf('Dataset composition:\n - %d images for predictor training\n - %d images for subsequent testing\n - %d best remaining images for PRNU estimation\n', pred_images, test_images, prnu_images);
    fprintf('The images will be saved as %d x %d TIFF bitmaps in %s\n', resolution(2), resolution(1), out_dirname);
    fprintf('\nPress any key to continue\n');
    pause();
    
    mkdir(out_dirname);
    mkdir([out_dirname 'prnu']);
    mkdir([out_dirname 'pristine']);
    mkdir([out_dirname 'predictor']);
    mkdir([out_dirname 'tampered-realistic']);
    
    % 1. Select images for 
    fi = 1;
    i = 1;
    skip = 0;
    while fi <= test_images + pred_images + skip
        fprintf('Processing image %d of %d (predictor and test images)\n', i, test_images + pred_images);
        image = imread([directory_path filenames(fi).name]);

        if size(image,1) > size(image,2)
            fprintf('Warning: %s not horizontal! - skipping\n', filenames(fi).name);
            skip = skip + 1;
            fi = fi + 1;
            continue;
        end
        
        if size(image,3) ~= 3
            fprintf('%s is not RGB - skipping\n', filenames(fi).name);
            skip = skip + 1;
            fi = fi + 1;
            continue;
        end

        image = imcropmiddle(image, resolution);
        
        [~, new_name, ~] = fileparts(filenames(fi).name);
        new_name = sprintf('%s.TIF', new_name);

        if i <= pred_images
            imwrite(image, [out_dirname 'predictor/' new_name]);
        else
            imwrite(image, [out_dirname 'pristine/' new_name]);
        end
        i = i + 1;
        fi = fi + 1;
    end
    
    % 2. Out of the remaining images, try to find the best ones for PRNU
    % estimation
    
    predictor_correlation = [-0.0046 0.0731 0.0657 -0.0273 0.4125 ...
    -0.0767 0.4435 0.0037 0.0177 -0.0485 0.1258 -1.3250 0.0082 ...
    -0.2177 0.8077]';

    predictions = nan(1, numel(filenames));
    
    last_fi = fi - 1;
    while fi <= numel(filenames)
        fprintf('Analyzing image %d of remaining %d (PRNU)\n', fi - last_fi, numel(filenames) - last_fi);
        image = imread([directory_path filenames(fi).name]);

        if size(image,1) > size(image,2)
            fprintf('Warning: %s not horizontal! - skipping\n', filenames(fi).name);
            fi = fi + 1;
            continue;
        end
        
        if size(image,3) ~= 3
            fprintf('%s is not RGB - skipping\n', filenames(fi).name);
            fi = fi + 1;
            continue;
        end
    
        image = imcropmiddle(image, resolution);
        f = extractPatchFeatures(image);
        predictions(fi) = f * predictor_correlation;
        fi = fi + 1;
    end

    if sum(~isnan(predictions)) < prnu_images
        fprintf('Warning: not enough images for PRNU estimation! using %d images\n', sum(~isnan(predictions)));
    end
    
    [~, sorted_indices] = sort(predictions, 'descend');
    % remove nans
    sorted_indices(isnan(predictions(sorted_indices))) = [];
    sorted_indices = sorted_indices(1:prnu_images);        
    i = 1;
    for fi = sorted_indices
        fprintf('Saving image %d of %d (PRNU)\n', i, prnu_images);
        image = imread([directory_path filenames(fi).name]);
        image = imcropmiddle(image, resolution);
        [~, new_name, ~] = fileparts(filenames(fi).name);
        new_name = sprintf('%s.TIF', new_name);
        imwrite(image, [out_dirname 'prnu/' new_name]);
        i = i + 1;
    end    

end