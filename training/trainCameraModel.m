function trainCameraModel(camera_name, analysis_windows, rebuild_predictors, options)
% trainCameraModel(camera_name, analysis_windows, rebuild_predictors, options)
%
% Trains a PRNU model of a specified camera that contains:
%  - estimate of the PRNU,
%  - a correlation predictor and signature presence model,
%  - signature absence model.
%
% The resulting camera model will be stored in a 'camera_model' structure and saved to 
% file ./data/camera_models/{camera_name}.mat
%
% Input parameters:
%   - camera_name        - self explanatory, must exist in ./data-input directory with the
%                          following sub-dir structure: {predictor,prnu}. See function
%                          `importDataset` for more information on importing data sets.
%   - analysis_windows   - vector with desired analysis window sizes, e.g., [64, 128]
%   - rebuild_predictors - if resuming camera model training, set to true to retrain all 
%                          predictors (not only the missing ones)
%   - options            - structure with advanced training options:
%
%      'patches_per_image' - number of patches extracted from a single image for predictor
%                            training
%      'seed'              - initializes the seed for random image patch selection
%      'jpeg_quality'      - if specified, the given quality level(s) will be used for 
%                            compressing the images before predictor training (PRNU 
%                            estimation remains unaffected). Possible values:
%                             a) scalar - use specific JPEG quality
%                             b) 2D vector - use random quality from the given range.
%                            Notes:
%                             I) if random quality is used, the training procedure will go
%                                through the data set 3 times to include the same image with
%                                different quality levels. This behavior can be customized by 
%                                hacking the code (look for field `nIters`).
%                             II) JPEGs are saved using Matlab's imwrite, so expect chrominance
%                                sub-sampling in the resulting images.
%      'noise_sigma'       - standard deviation for noise estimation
%      'noise_method'      - supported values:
%                            'mihcak' - standard, wavelet-based denoising filter
%                            'bm3d'   - BM3D (you need to have BM3D routines in your Matlab path)
%      'predictor_type'    - type of the correlation predictor:
%                            'lr' - linear regression
%                            'nn' - neural network
%      'feature_set'       - image patch features used for correlation prediction:
%                            0 - standard, 15-D feature sets
%                            1 - standard + normalized JPEG feature
%                            2 - only first-order features + normalized JPEG feature
%
% See also: viewCameraModel, importDataset
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

    % Set default parameter values
    if nargin < 2 || isempty(analysis_windows)
        analysis_windows = [32 48 64 96 128 192 256];
    end
    
    if nargin < 3 || isempty(rebuild_predictors)
        rebuild_predictors = false;
    end
    
    s = defaultPredictorSettings();
    if exist('options', 'var')
        s = applySettings(s, options);
    end
    
    % Setup correlation predictor
    s.patches_per_image = 500;
    s.seed = 2456;
    % If mixed JPEG compression is requested, increase the number of
    % iterations over the training set
    if isempty(s.jpeg_quality) || isscalar(s.jpeg_quality)
        s.nIters = 1;
    else
        s.nIters = 3;
    end
        
    % Basic setup
    dirname = './data-input/';
    camera_model_filename = sprintf('./data/camera_models/%s.mat', camera_name);

    % Estimation of camera's PRNU -----------------------------------------
    
    fprintf('\n############################# BUILDING CAMERA MODEL ################################\n\n')

    % If the specified model does not exist
    if ~exist(camera_model_filename, 'file')
        fprintf('Building camera model to %s\n', camera_model_filename);
        fprintf('Looking for images for PRNU estimation: %s\n', [dirname camera_name '/prnu/' '*.TIF']);
        filenames = dir([dirname camera_name '/prnu/' '*.TIF']);
        
        if numel(filenames) == 0
            fprintf('No images found, quitting...\n');
            return;
        end

        % List images available for PRNU estimation
        clear Images

        for i = 1:size(filenames,1)
            Images(i).name = [dirname camera_name '/prnu/' filenames(i).name];
        end

        % Populate basic structure of the camera model
        camera_model = struct();
        camera_model.name = camera_name;
        camera_model.stats = struct();        
        camera_model.stats.noise_sigma = s.noise_sigma;
        camera_model.stats.denoising = s.noise_method;
        
        if strcmpi(s.predictor_type, 'lr')
            camera_model.predictor_type = 'linear regression';
        elseif strcmpi(s.predictor_type, 'nn')
            camera_model.predictor_type = 'neural network';
        else
            throw(MException('prnu:InvalidParameter', 'Unsupported predictor type!'));
        end
        
        if strcmpi(s.noise_method, 'bm3d')            
            camera_model.prnu = estimateSPN(Images, 'mle', 'bm3d', s.noise_sigma);
        elseif strcmpi(s.noise_method, 'wavelet')            
            camera_model.stats.denoising = 'wavelet';
            camera_model.prnu = estimateSPN(Images, 'mle', 'mihcak', s.noise_sigma);
        else
            throw(MException('prnu:InvalidParameter', 'Unsupported denoising method!'));
        end
        
        camera_model.feature_set_code = s.feature_set;
        if s.feature_set == 0
            camera_model.feature_set = 'standard (15-D)';            
        elseif s.feature_set == 1
            camera_model.feature_set = 'standard + JPEG quality (16-D)';
        elseif s.feature_set == 2
            camera_model.feature_set = 'compact + JPEG quality (5-D)';
        else
            camera_model.feature_set = '(unknown)';
        end
        
        camera_model.stats.prnu_input_count = size(filenames,1);
        camera_model.stats.prnu_crop_mode = 'middle';
        if isempty(s.jpeg_quality)
            camera_model.stats.jpeg_preprocessing = 'none';
        elseif isscalar(s.jpeg_quality)
            camera_model.stats.jpeg_preprocessing = sprintf('jpeg(%d)', s.jpeg_quality);
        else
            camera_model.stats.jpeg_preprocessing = sprintf('jpeg(%d-%d)', s.jpeg_quality);
        end

        fprintf('Saving camera model to %s\n', camera_model_filename);
        save(camera_model_filename, 'camera_model');
    else
        fprintf('Using existing camera model from %s\n', camera_model_filename);
    end

    % Train Predictors ----------------------------------------------------

    for patch_size = analysis_windows
        % Check if the camera model already has this predictor
        [h1] = searchCameraModel(camera_name, patch_size);
        % If not, or if rebuilding is forced, rebuild
        if ~h1 || rebuild_predictors
            fprintf('Adding predictor for scale %d (%s) : ', patch_size, camera_model_filename);            
            predictor = trainPredictor(camera_model_filename, patch_size, s);
            fprintf('R2 = %.2f, max corr = %.2f\n', predictor.rsquare, max(predictor.raw_data));
            addPredictorToModel(camera_model_filename, predictor);
        end
    end

    % Train null model ----------------------------------------------------

    for patch_size = analysis_windows
        % Check if the camera model already has this predictor
        [~, h0] = searchCameraModel(camera_name, patch_size);
        % If not, or if rebuilding is forced, rebuild
        if ~h0 || rebuild_predictors
            fprintf('Adding null model for scale %d (%s) : ', patch_size, camera_model_filename);
            model = trainNullModel(camera_model_filename, patch_size, s);
            fprintf('max corr = %.2f\n', max(model.raw_data));
            addDistToModel(camera_model_filename, model);
        end
    end
end

function [foundPredictor, foundNullModel] = searchCameraModel(camera_model, patch_size)
    foundNullModel = false;
    foundPredictor = false;
    
    if ischar(camera_model)
        if ~exist(camera_model, 'file')
            if regexp(camera_model, '\.mat$') <= 0
                camera_model = sprintf('%s.mat', camera_model);
            end
            camera_model = sprintf('./data/camera_models/%s', camera_model);
        end 
        load(camera_model);
    end
    
    if isfield(camera_model, 'predictors')    
        for i = 1:numel(camera_model.predictors)
            if camera_model.predictors{i}.patch_size == patch_size
                foundPredictor = true;
                break
            end
        end
    end
    
    if isfield(camera_model, 'absence_models')
        for i = 1:numel(camera_model.absence_models)
            if camera_model.absence_models{i}.patch_size == patch_size
                foundNullModel = true;
                break
            end
        end
    end
end

function addPredictorToModel(camera_model_filename, predictor)
    load(camera_model_filename);
    if ~isfield(camera_model, 'predictors')
        camera_model.predictors = cell(1,1);
        camera_model.predictors{1} = predictor;
    else
        found = false;
        for i = 1:numel(camera_model.predictors)
            if camera_model.predictors{i}.patch_size == predictor.patch_size
                camera_model.predictors{i} = predictor;
                found = true;
                break
            end
        end
        if found == false
            camera_model.predictors{numel(camera_model.predictors)+1} = predictor;
        end
    end
    save(camera_model_filename, 'camera_model');
end

function addDistToModel(camera_model_filename, predictor)
    load(camera_model_filename);
    if ~isfield(camera_model, 'absence_models')
        camera_model.absence_models = cell(1,1);
        camera_model.absence_models{1} = predictor;
    else
        found = false;
        for i = 1:numel(camera_model.absence_models)
            if camera_model.absence_models{i}.patch_size == predictor.patch_size
                camera_model.absence_models{i} = predictor;
                found = true;
                break
            end
        end
        if found == false
            camera_model.absence_models{numel(camera_model.absence_models)+1} = predictor;
        end
    end
    save(camera_model_filename, 'camera_model');
end