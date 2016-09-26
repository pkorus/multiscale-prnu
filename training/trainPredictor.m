function predictor = trainPredictor(camera_model_filename, patch_size, advanced)
% trainPredictor(camera_model_filename, patch_size, patches_per_image, advanced)
%
% Trains a PRNU correlation & PCE predictor for a given camera model and analysis
% window size. The predictor uses linear regression on features extracted
% from the analysis window using extractPatchFeatures function.
%
% For details about the supported advanced settings, see function
% defaultPredictorSettings
%
% The structure of the predictor for an example camera:
%           patch_size: 64                        basic info
%          correlation: [15x1 double]             lin. regr. coeffs.
%              rsquare: 0.781                     prediction quality R^2
%              normfit: [-2.4157e-12 0.0359]      params of a Gaussian fit
%    number_of_patches: 25000                     stats
%             raw_data: [25000x1 double]          raw data for debugging
%      raw_predictions: [25000x1 double]          raw data for debugging
%
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

    % Load the camera model & extract camera name
    m = load(camera_model_filename);
    [~, camera_name, ~] = fileparts(camera_model_filename);
    camera_model = m.camera_model;
    
    s = defaultPredictorSettings();
    if exist('advanced', 'var')
        s = applySettings(s, advanced);
    end    
            
    % Setup the pseudo-random number generator
    
    if s.seed <= 0
        s.seed = round(prod(clock()));
    end
    rand('twister', s.seed);
    
    % Loop setup
    patch_id = 1;    

    % List images available for predictor training
    filenames = dir([dirname camera_name '/predictor/' '*.TIF']);

    % Get camera's PRNU
    prnu = rgb2gray1(camera_model.prnu);
    prnu = WienerInDFT(prnu,std2(prnu));

    % Data collection variables
    if s.feature_set == 0
        nFeatures = 15;
    elseif s.feature_set == 1
        nFeatures = 16;
    elseif s.feature_set == 2
        nFeatures = 5;
    else
        throw(MException('prnu:InvalidParameter', 'Unsupported feature set!'));
    end
    
    batch_size = s.patches_per_image * numel(filenames);
    features = zeros(s.nIters * batch_size, nFeatures);
    target = zeros(s.nIters * batch_size,1);

    for iter = 1:s.nIters
    for k = 1:numel(filenames)

        % Load the image
        image = imread([dirname camera_name '/predictor/' filenames(k).name]);

        % Compress the image, if desired
        if ~isempty(s.jpeg_quality)
            if isscalar(s.jpeg_quality)
                tmpname = ['/tmp/' filenames(k).name '.jpg'];
                imwrite(image, tmpname, 'Quality', s.jpeg_quality)
                image = imread(tmpname);
            else
                jpq_local = min(s.jpeg_quality) + randi(max(s.jpeg_quality) - min(s.jpeg_quality));
                tmpname = ['/tmp/' filenames(k).name '.jpg'];
                imwrite(image, tmpname, 'Quality', jpq_local)
                image = imread(tmpname);
            end
        end
        
        if exist('jpq_local', 'var')
            jpeg_feature = (100 - jpq_local) / 50;
        else
            jpeg_feature = 0;
        end
        
        % Extract PRNU estimate from the image
        if strcmpi(s.noise_method, 'bm3d')        
            image_noise = getPRNU(image, s.noise_sigma);
        elseif strcmpi(s.noise_method, 'wavelet')
            image_noise = NoiseExtractFromImage(image, s.noise_sigma);
        else
            throw(MException('prnu:InvalidParameter', 'Unsupported denoising filter!'));
        end
        image_noise = ZeroMeanTotal(image_noise);
        image_noise = WienerInDFT(image_noise, std2(image_noise));

        % Collect data for predictor training
        for i = 1:s.patches_per_image
            % Choose random location
            bx = 1 + round(rand() * (size(image,2) - patch_size - 1));
            by = 1 + round(rand() * (size(image,1) - patch_size - 1));
            
            % Extract the patch of the image and its PRNU
            patch   = image(by:by+patch_size-1, bx:bx+patch_size-1, :);
            i_noise = image_noise(by:by+patch_size-1, bx:bx+patch_size-1, :);
            
            % Extract the corresponding patch of the camera PRNU
            c_prnu = prnu(by:by+patch_size-1, bx:bx+patch_size-1, :);
            c_prnu = c_prnu .* double(rgb2gray(patch));

            % Make sure the noise is 0-mean
            i_noise = i_noise - mean2(i_noise);
            c_prnu = c_prnu - mean2(c_prnu);
            
            % Compute normalization factors
            n1 = sqrt(sum2(i_noise .* i_noise));
            n2 = sqrt(sum2(c_prnu .* c_prnu));

            % Extract features of the patch
            if s.feature_set == 0
                f = extractPatchFeatures(patch);                        
            elseif s.feature_set == 1
                f = [extractPatchFeatures(patch) jpeg_feature];                                    
            elseif s.feature_set == 2
                f = [extractPatchFeatures(patch) jpeg_feature];                                    
                f = f([2,3,4,5,16]);
            end
            features(patch_id,:) = f;
          
            % Compute normalized correlation
            target(patch_id) = sum2(i_noise .* c_prnu) / n1 / n2;
            
            % Go to the next patch
            patch_id = patch_id + 1;
        end
    end
    end
    
    % Populate the predictor model
    predictor = struct();
    predictor.patch_size = patch_size;
    predictor.type = s.predictor_type;
    
    % Prediction of Correlation Scores - Least Squares Fit
    if strcmpi(s.predictor_type, 'lr')
        params = inv(features'*features)*features'*target;
        predictions = features * params;
        predictor.correlation = params;
    elseif strcmpi(s.predictor_type, 'nn')
        net = feedforwardnet(5);
        net = configure(net, features', target');
        net.trainParam.showWindow = 0;
        net = train(net, features', target');
        predictions = net(features')';
        predictor.correlation = net;
    end

    % Fit a Gaussian
    [mu, sig] = normfit(target - predictions);

    % Populate the predictor model    
    predictor.rsquare = rsquare(target, predictions);
    predictor.normfit = [mu, sig];
    predictor.number_of_patches = s.patches_per_image * numel(filenames) * s.nIters;
    predictor.raw_data = target;
    predictor.raw_predictions = predictions;
end