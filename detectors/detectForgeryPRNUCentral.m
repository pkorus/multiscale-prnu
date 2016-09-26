function response_map = detectForgeryPRNUCentral(image, camera_model, bs, advanced)
% response_map = detectForgeryPRNUCentral(image, camera_model, bs, advanced)
%
% Tampering localization based on PRNU analysis. The function implements a
% sliding-window detector with central pixel attribution. Optionally, the
% correlation statistics can be limited to the central image segment. For
% more information, see [1].
% 
% Params: 
%  - image           - RGB image, uint8
%  - camera_model    - camera model structure
%  - bs              - odd block size, e.g., 129
%  - advanced        - advanced parameters:
%     stride                  - stride, e.g., 8
%     cfar_std_mul            - tampering acceptance threshold (number of 
%                               standard deviations), default 2
%     cfar_fp_rate            - false alarm rate for CFAR map, default 0.01
%     verbose                 - boolean
%     calc_pce                - include also PCE and p-value maps
%     calc_cfar               - include also a CFAR map (3-label, see below)
%     segmentwise_correlation - corr only over similar pixels
%     similarity_threshold    - threshold for pixel similarity
%     minimum_window_size     - minimum segment size (in pixels)
%     image_padding           - mirror-pad image before processing
%
% The response_map structure will contain the following fields:
%  - map_cor         - correlation field
%  - est_cor         - predicted correlation field
%  - map_prp         - tampering probability with real-valued scores [0,1]
%  - map_cfar        - binary decision map according to the two-threshold CFAR method
% and optionally:
%  - map_pva         - p-value based response map
%  - map_pce         - PCE response map
%
% Notes:
%
% The optional CFAR map is generated with 3 labels: 0 - authentic pixels; 1 - potentially
% tampered pixels that do not meet the constraint on the desired false positive rate;
% 2 - potentially tampered pixels that satisfy the FAR constraint.
%
% References:
%
% [1] P. Korus, J. Huang, Multi-scale Analysis Strategies in PRNU-based 
%     Tampering Localization, Submitted to IEEE Tran. Information Forensics
%     and Security
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

    % Sanitize parameters
    if mod(bs,2) ~= 1
        error('prnu:InvalidParameter', 'Invalid block size - needs to be odd!');
    end
    
    % Unpack advanced parameters
    if isfield(advanced, 'cfar_fp_rate')
        cfar_fp_rate = advanced.cfar_fp_rate;
    else
        cfar_fp_rate = 1e-2;
    end

    if isfield(advanced, 'cfar_std_mul')
        cfar_std_mul = advanced.cfar_std_mul;
    else
        cfar_std_mul = 2;
    end
    
    if isfield(advanced, 'verbose')
        verbose = advanced.verbose;
    else
        verbose = true;
    end

    if isfield(advanced, 'calc_pce')
        calc_pce = advanced.calc_pce;
    else
        calc_pce = false;
    end

    if isfield(advanced, 'calc_cfar')
        calc_cfar = advanced.calc_cfar;
    else
        calc_cfar = false;
    end
    
    if isfield(advanced, 'segmentwise_correlation')
        segmentwise_correlation = advanced.segmentwise_correlation;
    else
        segmentwise_correlation = false;
    end

    if isfield(advanced, 'image_padding')
        image_padding = advanced.image_padding;
    else
        image_padding = false;
    end
    
    if isfield(advanced, 'stride')
        bsk = advanced.stride;
    else
        bsk = 16;
    end
    
    if isfield(advanced, 'similarity_threshold')
        similarity_threshold = advanced.similarity_threshold;
    else
        similarity_threshold = 20;
    end
    
    if isfield(advanced, 'minimum_window_size')
        minimum_window_size = advanced.minimum_window_size;
    else
        minimum_window_size = 64 * 64;
    end          
    
    % Start processing
    if verbose
      fprintf('Extracting image PRNU... '); tic;
    end
    
    if ischar(image)
        image = imread(image);
    end
    
    if image_padding
        image = padarray(image, floor(bs/2)*[1 1], 'symmetric', 'both');
    end
    
    % Estimate the PRNU for the investigated image
    if strcmp(camera_model.stats.denoising, 'bm3d')
        noise = getPRNU(image, camera_model.stats.noise_sigma);
    else
        noise = NoiseExtractFromImage(image, camera_model.stats.noise_sigma);
    end
    noise = ZeroMeanTotal(noise);
    noise = WienerInDFT(noise, std2(noise));
    
    if ~isfield(camera_model, 'feature_set_code') || camera_model.feature_set_code == 0
        feature_set = 0;
    elseif camera_model.feature_set_code == 1
        feature_set = 1;
        jpeg_feature = (100 - estimateQFactor(image)) / 50;
    elseif camera_model.feature_set_code == 2
        feature_set = 2;
        jpeg_feature = (100 - estimateQFactor(image)) / 50;
    else
        throw(MException('prnu:InvalidParameter', 'Unsupported feature set in the camera model!'))
    end
    
    if ~isfield(camera_model, 'predictor_type') || strcmpi(camera_model.predictor_type, 'linear regression')
        nn = false;
    elseif strcmpi(camera_model.predictor_type, 'neural network')
        nn = true;
    else
        throw(MException('prnu:InvalidParameter', 'Unsupported predictor type!'))
    end
    
    % Get camera PRNU from the model
    image = double(rgb2gray(image));
    std_image = stdfilt3(image);
    hpf_std_image = stdfilt3(highPass(image));
    prnu = rgb2gray1(camera_model.prnu);
    prnu = WienerInDFT(prnu,std2(prnu));
    if image_padding
        prnu = padarray(prnu, floor(bs/2)*[1 1], 'symmetric', 'both');
    end
    prnu = prnu .* image;

    if verbose
      t_last = toc; fprintf('done in %.2f min\n', t_last/60);
    end

    % Helper functions / variables
    sum2 = @(x) sum(x(:));
            
    bs2 = floor(bs/2);
    maxbx = size(image,2) - bs;
    maxby = size(image,1) - bs;
            
    maxfx = ceil((size(image,2) - bs + 1)/bsk);
    maxfy = ceil((size(image,1) - bs + 1)/bsk);    
    
    % Structural elements for segmentation-guided strategy
    seg_strel = strel('disk', 8);
    inc_strel = strel('disk', 4);
        
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
            using_predictor = true;
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

    % Extract some basic properties of the models for successive analysis 
    % windows - will be used when guiding the correlation with image
    % content
    w = zeros(1, numel(camera_model.predictors));
    v = zeros(1, numel(camera_model.predictors));
    b = zeros(1, numel(camera_model.predictors));
    for i = 1:numel(w); 
        b(i) = camera_model.predictors{i}.patch_size;
        w(i) = camera_model.predictors{i}.normfit(2); 
        v(i) = camera_model.absence_models{i}.normfit(2); 
    end;
    

    % Setup which maps should be included
    map_list = {'map_cor', 'map_prp', 'est_cor'};
    
    if calc_pce
      map_list = { map_list{:} 'map_pce', 'map_pva'};
    end
    
    if calc_cfar
        map_list = { map_list{:} 'map_cfar'};
    end
        
    % Allocate output structure
    response_map = struct();
    for map_name = map_list
        response_map.(map_name{1}) = nan(size(image,1), size(image,2));
    end
    
    if verbose
      fprintf('Detecting forgeries: '); tic;
    end

    % Start processing
    bx = 1;
    fi = 1;

    last_progress = -1;
    total = maxfx * maxfy;

    while bx <= maxbx
        by = 1;
        while by <= maxby
            if verbose
                counter = fi;
                progress = round(100*(counter/total)/2);
                % Display progress bar
                if progress ~= last_progress
                    displayProgress(sprintf('Window   %6d / %6d', fi, total), progress, 50);
                end
                last_progress = progress;
            end
        
            % Extract patch features
            l_patch = image((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));
            std_patch = std_image((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));
            hpf_std_patch = hpf_std_image((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));
                        
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
            
            % Extract local patches for verification
            l_noise = noise((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));
            l_prnu = prnu((by-1)+1:(by+bs-1), (bx-1)+1:(bx+bs-1));

            % Make sure the noise is 0-mean
            l_noise = l_noise - mean2(l_noise);
            l_prnu = l_prnu - mean2(l_prnu);
                        
            % Peak Correlation Energy
            if calc_pce
                Corr = crosscorr(l_noise, l_prnu);
                C = PCE(Corr);            
                response_map.map_pce(by + bs2, bx + bs2) = C.PCE;
                response_map.map_pva(by + bs2, bx + bs2) = C.pvalue;
            end
            
            % Normalized correlation
            if segmentwise_correlation
                % choose pixels similar to the central pixel
                similar_pixels = get_similarpixels(l_patch, similarity_threshold, seg_strel);

                % If the similar region is too small, expand it
                min_eff_blk_size = minimum_window_size;
                % Optionally, consider adapting based on the estimated corr               
                % Computed from 2*sigma rule for H0 hypothesis
%                 min_eff_blk_size = max(minimum_window_size, min(bs-1, 2 / pred_cor))^2;
                
                counter = 1;
                while sum(similar_pixels(:)) < min_eff_blk_size
                    similar_pixels = imdilate(similar_pixels, inc_strel);
                    counter = counter + 1;
                end
                
                % Contract the patches to pixels of interest only
                l_noise = l_noise(similar_pixels);
                l_prnu  =  l_prnu(similar_pixels);
                                
                % Compute the effective size of the window (in pixels)
                eff_bs = sqrt(sum(similar_pixels(:)));
            end
            
            % Compute normalization factors
            n1 = sqrt(sum2(l_noise .* l_noise));
            n2 = sqrt(sum2(l_prnu .* l_prnu));

            % Compute the normalized correlation
            ncorr = sum2(l_noise .* l_prnu) / n1 / n2;
            
            % Populate output structures
            response_map.map_cor(by + bs2, bx + bs2) = ncorr;
            response_map.est_cor(by + bs2, bx + bs2) = pred_cor;

            % Standard CFAR-based detection (as in the original approach)
            % the first threshold (here, fixed at 2 * standard deviation - approx 1% of false miss rate
            % the second threshold is dynamic - guarantees certain fale acceptance rate limit
            if calc_cfar
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
                response_map.map_cfar(by + bs2, bx + bs2) = label;
            end
            
            % Compute tampering probability
            mu_0 = 0;
            mu_1 = pred_cor;
            si_0 = absence_model.normfit(2);
            si_1 = predictor.normfit(2);
            
            % Adapt variances if the correlation scope has been contracted
            if segmentwise_correlation
                si_1 = spline(b, w, eff_bs);
                si_0 = spline(b, v, eff_bs);
                % Alternatively, use si_0 = 1/eff_bs;
            end
                        
            response_map.map_prp(by + bs2, bx + bs2) = corr2Response(ncorr, mu_0, si_0, mu_1, si_1);
            
            % Move window
            by = by + bsk;
            fi = fi + 1;
        end
        bx = bx + bsk;
    end

    if verbose
      t_last = toc; fprintf(' done in %.2f min\n', t_last/60);
    end
    
    % Compact missing pixels due to stride > 1
    if bsk > 1
        for map_name = map_list
            response_map.(map_name{1}) = response_map.(map_name{1})(mod(bs2-1,bsk)+2:bsk:end, mod(bs2-1,bsk)+2:bsk:end);
            response_map.(map_name{1}) = fillBorders(response_map.(map_name{1}));
            response_map.(map_name{1}) = single(response_map.(map_name{1}));
            if image_padding
                margin = floor(bs/2) / bsk;
                map = response_map.(map_name{1});
                map = map(margin:end-margin, margin:end-margin);                
                response_map.(map_name{1}) = map;
            end
        end
    end

    % Post-process the CFAR-based map
    if calc_cfar
        map_cfar = response_map.map_cfar == 2;
        map_cfar = imresize(map_cfar, [size(image,1), size(image,2)], 'nearest');
        map_cfar = mapCleanup(map_cfar, (bs/2)^2);
        map_cfar = imdilate(map_cfar, strel('square', 20));
        response_map.map_cfarb = imresize(map_cfar, size(response_map.map_cfar), 'nearest');
        response_map.map_cfar = uint8(response_map.map_cfar);
    end
    response_map.block = bs;
end

function similar_pixels = get_similarpixels(l_patch, similarity_threshold, seg_strel)

    bs = size(l_patch,1);
    bs2 = floor(bs / 2);
    cpx = l_patch(bs2 + 1 ,bs2 + 1, :);
    similar_pixels = mean(abs(repmat(cpx, [bs, bs]) - l_patch),3) < similarity_threshold;
    similar_pixels = imclose(similar_pixels, seg_strel);

    % include only components that contain the central pixel
    [labels, nLabels] = bwlabel(similar_pixels);
    if nLabels > 1
        i = 1;
        while i <= nLabels
            mask = labels == i;
            if mask(bs2 + 1 ,bs2 + 1) == 1
                similar_pixels = mask;
            end
            i = i + 1;
        end
    end
end