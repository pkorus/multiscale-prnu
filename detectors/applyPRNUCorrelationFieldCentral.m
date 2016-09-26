function [map_cfar, map_prp, est_cor, map_cor] = applyPRNUCorrelationFieldCentral(image, camera_model, corr_field, bs, bsks, rates, verbose)
% applyPRNUCorrelationFieldCentral(image, camera_model, corr_field, bs, bsks, rates)
%
% Estimates correlation with a PRNU signature using the camera's predictor and 
% generates a tampering probability map and a conventional CFAR decision map based
% on the provided correlation field.
%
% This function is not intended for individual use. It is a helper function for 
% localization algorithms implemented based on image filtering. See docs for
% function `detectForgeryPRNUFilter` for more information.
% 
% Params:
%  - image           - RGB image, uint8
%  - camera_model    - camera model structure
%  - corr_field      - previously obtained actual correlation field
%                      the correlation is expected not to be normalized;
%                      due to problems with global normalization,
%                      local normalization will be applied here
%  - bs              - block size, e.g., 128 (only some values are allowed)
%  - bsks            - strides, separately for correlation estimation and for
%                      response computation, e.g., [8 4]
%  - rates           - 2D vector with CFAR thresholds: 
%                      1. tampering acceptance threshold (number of 
%                         standard deviations), default 2
%                      2. false alarm rate for CFAR map, default 0.01
%
% Output:
%   map_cfar         - standard CFAR map
%   map_prp          - tampering probability map
%   est_cor          - estimated correlation map
%   map_cor          - normalized correlation field
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
    
    if mod(bs,2) ~= 1
        error('prnu:InvalidParameter', 'Invalid block size - needs to be odd!');
    end

    if nargin < 6 || isempty(rates)
        rates = [2 1e-2];
    end
    
    if isscalar(rates)
      rates = [2 rates];
    end
    
    if nargin < 7 || isempty(verbose)
        verbose = false;
    end
    
    if isscalar(bsks)
        bsks = [bsks bsks];
    end

    if size(image,3) > 1
        image = rgb2gray(image);
    end

    sum2 = @(x) sum(x(:));

    % Prepare filtered images beforehand to speed up feature extraction
    std_image = stdfilt3(image);
    hpf_std_image = stdfilt3(highPass(image));
        
    % Extract noise from the image
    if strcmp(camera_model.stats.denoising, 'bm3d')
    	noise = getPRNU(image, camera_model.stats.noise_sigma);
    else
    	noise = NoiseExtractFromImage(image, camera_model.stats.noise_sigma);
    end
    noise = WienerInDFT(noise, std2(noise));
    
    % Read settings from the camera model    
    if ~isfield(camera_model, 'feature_set_code') || camera_model.feature_set_code == 0
        feature_set = 0;
    elseif camera_model.feature_set_code == 1
        feature_set = 1;
        jpeg_feature = (100 - estimateQFactor(image)) / 50;
    elseif camera_model.feature_set_code == 2
        feature_set = 2;
        jpeg_feature = (100 - estimateQFactor(image)) / 50;
    else
        throw(MException('PRNU:BadParameter', 'Unsupported feature set in the camera model!'))
    end  

    if ~isfield(camera_model, 'predictor_type') || strcmpi(camera_model.predictor_type, 'linear regression')
        nn = false;
    elseif strcmpi(camera_model.predictor_type, 'neural network')
        nn = true;
    else
        throw(MException('PRNU:BadParameter', 'Unsupported predictor type!'))
    end
    
    % Fetch camera's PRNU
    image = double(image);
    prnu = rgb2gray1(camera_model.prnu);
    prnu = WienerInDFT(prnu,std2(prnu));
    prnu = prnu .* image;

    % Set stride for correlation estimation
    bsk = bsks(1);

    bs2 = floor(bs/2);
    maxbx = size(image,2) - bs;
    maxby = size(image,1) - bs;
            
    maxfx = ceil((size(image,2) - bs + 1)/bsk);
    maxfy = ceil((size(image,1) - bs + 1)/bsk);
            
    % Check if the camera has an appropriate predictor
    predictor = [];
    using_predictor = false;
    if isfield(camera_model, 'predictors')
        for p = 1:numel(camera_model.predictors)
            if camera_model.predictors{p}.patch_size == bs
                predictor = camera_model.predictors{p};
                using_predictor = true;
            end
        end
    end
    if ~using_predictor
        available = zeros(1, numel(camera_model.predictors));
        for p = 1:numel(camera_model.predictors)
            available(p) = camera_model.predictors{p}.patch_size;
        end
        dists = abs(available - bs);
        [md, p] = min(dists);
        if md < 8
            predictor = camera_model.predictors{p};
            if md > 4
                fprintf('WARNING No suitable predictor found (%s, %d px) substituting with %d\n', camera_model.name, bs, available(p));
            end
        else
            fprintf('WARNING Missing predictor (%s, %d px)', camera_model.name, bs);
        end
    end

    % Check if the camera has an appropriate PRNU absence model
    absence_model = [];
    if isfield(camera_model, 'absence_models')
        available = zeros(1, numel(camera_model.predictors));
        for p = 1:numel(camera_model.absence_models)
            available(p) = camera_model.absence_models{p}.patch_size;
        end
        dists = abs(available - bs);
        [md, p] = min(dists);
        if md < 8
            absence_model = camera_model.absence_models{p};
            if md > 4
                fprintf('WARNING No suitable predictor found (%s, %d px) substituting with %d\n', camera_model.name, bs, available(p));
            end
        else
            fprintf('WARNING Missing PRNU absence model (%s, %d px)\n', camera_model.name, bs);
            absence_model = struct();
            absence_model.normfit = [0 0.01];
        end

    end        

    % Allocate output structure
    est_cor = nan(size(image,1), size(image,2));
    div_factor = nan(size(image,1), size(image,2));
    map_cfar = nan(size(image,1), size(image,2));
    map_prp = nan(size(image,1), size(image,2));

    if verbose
        tic;
    end

    bx = 1;
    last_progress = -1;
    total = maxfx * maxfy;
    fi = 1;
    bsl = round(bs / 8) * 8;
    while bx <= maxbx
        by = 1;
        while by <= maxby
            if verbose
                % Simple progress display
                counter = fi;
                progress = round(100*(counter/total)/2);
                % Display progress bar
                if progress ~= last_progress
                    displayProgress(sprintf('Estimating correlation (window %6d / %6d)', fi, total), progress, 50);
                end
                last_progress = progress;       
            end

            % Estimate correlation
            l_patch       =         image((by-1)+1:(by+bsl-1), (bx-1)+1:(bx+bsl-1));
            std_patch     =     std_image((by-1)+1:(by+bsl-1), (bx-1)+1:(bx+bsl-1));
            hpf_std_patch = hpf_std_image((by-1)+1:(by+bsl-1), (bx-1)+1:(bx+bsl-1));            
            if feature_set == 0
                features = extractPatchFeaturesFast(l_patch, std_patch, hpf_std_patch);
            elseif feature_set == 1
                features = [extractPatchFeaturesFast(l_patch, std_patch, hpf_std_patch) jpeg_feature];
            elseif feature_set == 2
                features = [extractPatchFeaturesFast(l_patch, std_patch, hpf_std_patch) jpeg_feature];                
                features = features([2 3 4 5 16]);
            end                        
            % Predict correlation
            if nn == true
                pred_cor = max(0, predictor.correlation(features'));
            else
                pred_cor = max(0, features * predictor.correlation);
            end
                
            est_cor(by + bs2, bx + bs2) = pred_cor;
            
            % Estimate local normalization factors
            l_noise = noise((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));
            l_prnu = prnu((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));

            % Make sure the noise is 0-mean
            l_noise = l_noise - mean2(l_noise);
            l_prnu = l_prnu - mean2(l_prnu);

            % Compute normalization factors
            n1 = sqrt(sum2(l_noise .* l_noise));
            n2 = sqrt(sum2(l_prnu .* l_prnu));

            % Normalized correlation
            div_factor(by + bs2, bx + bs2) =  n1 * n2;

            % Move window
            by = by + bsk;
            fi = fi + 1;
        end
        bx = bx + bsk;
    end
    
    if verbose
        t_last   = toc; fprintf(' done in %.2f min', t_last/60);
    end

    % Compact missing pixels due to stride > 1
    if bsk > 1
        est_cor  = est_cor(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end);
        div_factor  = div_factor(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end);
    end

    % Fill missing border values
    est_cor = fillBorders(est_cor);
    est_cor = imresize(est_cor, size(corr_field));
    div_factor = fillBorders(div_factor);
    div_factor = imresize(div_factor, size(corr_field));

    if verbose
        tic;
    end

    % Set stride for response map computation
    bs_org = bs;
    bs = 2*bsks(2) - 1;
    bsk = bsks(2);

    bs2 = floor(bs/2);
    maxbx = size(image,2) - bs;
    maxby = size(image,1) - bs;

    maxfx = ceil((size(image,2) - bs + 1)/bsk);
    maxfy = ceil((size(image,1) - bs + 1)/bsk);

    bx = 1;
    last_progress = -1;
    total = maxfx * maxfy;
    fi = 1;
    while bx <= maxbx
        by = 1;
        while by <= maxby
            if verbose
                % Simple progress display
                counter = fi;
                progress = round(100*(counter/total)/2);
                % Display progress bar
                if progress ~= last_progress
                    displayProgress(sprintf('Computing decisions    (window %6d / %6d)', fi, total), progress, 50);
                end
                last_progress = progress;       
            end

            pred_cor = est_cor(by + bs2, bx + bs2);
            ncorr = corr_field(by + bs2, bx + bs2) / div_factor(by + bs2, bx + bs2);

            % Standard CFAR-based detection (as in original paper)
            threshold = rates(1) * absence_model.normfit(2);
            label = 0;
            if ncorr < threshold
                % Check the H1 distribution to see if fp limit is not
                % exceeded                
                if normcdf(ncorr, pred_cor, predictor.normfit(2)) < rates(2)
                    label = 2;
                else
                    label = 1;
                end                                
            end
            map_cfar(by + bs2, bx + bs2) = label;

            % Tampering probability
            mu_0 = 0;
            mu_1 = pred_cor;
            si_0 = absence_model.normfit(2);
            si_1 = predictor.normfit(2);
                        
            map_prp(by + bs2, bx + bs2) = corr2Response(ncorr, mu_0, si_0, mu_1, si_1);            

            % Move window
            by = by + bsk;
            fi = fi + 1;
        end
        bx = bx + bsk;
    end
    
    if verbose
        t_last   = toc; fprintf(' done in %.2f min\n', t_last/60);
    end

    % Compact results and fill missing
    if bsk > 1
        map_cfar   = map_cfar(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end);
        map_prp    = map_prp(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end);
    end
    map_cfar = fillBorders(map_cfar);
    map_prp = fillBorders(map_prp);    

    % Return compensated correlation field
    map_cor = corr_field ./ div_factor;

    % Post-process the CFAR decision map
    map_cfar = imresize(map_cfar, [size(image,1), size(image,2)], 'nearest');
    map_cfar = mapCleanup(map_cfar, (bs_org/2)^2);
    map_cfar = imdilate(map_cfar, strel('square', 20));
    map_cfar = imresize(map_cfar, size(est_cor), 'nearest');
    map_cfar = im2bw(map_cfar);

    % Compact maps for smaller storage requirements
    est_cor = single(est_cor(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end));
    map_prp = single(map_prp);
end