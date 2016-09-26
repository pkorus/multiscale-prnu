function null_model = trainNullModel(camera_model_filename, patch_size, advanced)
% trainNullModel(camera_model_filename, patch_size, patches_per_image)
%
% Trains a null model of PRNU correlation (when the PRNU is absent) for a
% given camera and analysis window size. The absence model will be based
% on different patches from the same image.
%
% The structure of the null model for an example camera:
%           patch_size: 64
%              normfit: [-5.4382e-05 0.0107]
%    number_of_patches: 25000
%             raw_data: [25000x1 double]
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

    % Helper functions / variables
    sum2 = @(x) sum(x(:));
    dirname = './data-input/';

    % Load the camera model
    m = load(camera_model_filename);
    [~, camera_name, ~] = fileparts(camera_model_filename);
    camera_model = m.camera_model;
    
    s = defaultPredictorSettings();
    if exist('advanced', 'var')
        s = applySettings(s, advanced);
    end   
    
    if s.seed <= 0
        s.seed = round(prod(clock()));
    end
    rand('twister', s.seed);

    % Loop setup
    patch_id = 1;    

    % List images available for predictor training
    filenames = dir([dirname camera_name '/predictor/' '*.TIF']);

    % Data collection variables
    target = zeros(s.patches_per_image * numel(filenames),1);

    for k = 1:numel(filenames)

        % Load image
        image = imread([dirname camera_name '/predictor/' filenames(k).name]);

        % Get its PRNU - either calculate or load from cache
        if strcmpi(s.noise_method, 'bm3d')        
            image_noise = getPRNU(image, s.noise_sigma);
        elseif strcmpi(s.noise_method, 'wavelet')
            image_noise = NoiseExtractFromImage(image, s.noise_sigma);
        else
            throw(MException('prnu:InvalidParameter', 'Unsupported denoising filter!'));
        end
        image_noise = WienerInDFT(image_noise, std2(image_noise));

        % Collect correlation data for patches at different locations
        for i = 1:s.patches_per_image
            % Pick random location
            bx = 1 + round(rand() * (size(image,2) - patch_size - 1));
            by = 1 + round(rand() * (size(image,1) - patch_size - 1));
            
            % Pick a different location for PRNU mismatch
            px = bx;
            py = by;
            
            while px == bx
                px = 1 + round(rand() * (size(image,2) - patch_size - 1));
            end
            while py == by
                py = 1 + round(rand() * (size(image,1) - patch_size - 1));
            end

            % Extract the image patch and its PRNU at location 1
            patch   = image(by:by+patch_size-1, bx:bx+patch_size-1, :);
            i_noise = image_noise(by:by+patch_size-1, bx:bx+patch_size-1, :);

            % Extract the patch of the camera's PRNU at location 2
            c_prnu = camera_model.prnu(py:py+patch_size-1, px:px+patch_size-1, :);
            % Multiply the PRNU by luminance from location 1 - the location
            % that will be verified
            c_prnu = rgb2gray1(c_prnu) .* double(rgb2gray(patch));

            % Make sure the noise is 0-mean
            i_noise = i_noise - mean2(i_noise);
            c_prnu = c_prnu - mean2(c_prnu);
            
            % Compute normalization factors
            n1 = sqrt(sum2(i_noise .* i_noise));
            n2 = sqrt(sum2(c_prnu .* c_prnu));

            target(patch_id) = sum2(i_noise .* c_prnu) / n1 / n2;

            patch_id = patch_id + 1;
        end

    end

    % Fit a Gaussian
    [mu, sig] = normfit(target);

    % Populate the null model
    null_model = struct();
    null_model.patch_size = patch_size;
    null_model.normfit = [mu, sig];
    null_model.number_of_patches = s.patches_per_image * numel(filenames);
    null_model.raw_data = target;
end
