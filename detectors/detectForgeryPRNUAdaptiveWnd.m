function response_map = detectForgeryPRNUAdaptiveWnd(image, camera_model, advanced)
% response_map = detectForgeryPRNUAdaptiveWnd(image, camera_model, advanced)
%
% Tampering localization based on PRNU analysis with adaptive window size and
% central pixel attribution. Initially, small windows are used, and window size
% is increased if smaller scale analysis turns out to be unreliable. For more
% information, see [1].
%
% 
% Params: 
%  - image           - RGB image, uint8
%  - camera_model    - camera model structure
%  - advanced        - advanced parameters:
%     stride                     - stride, e.g., 8
%     verbose                    - boolean
%     calc_pce                   - include also PCE and p-value maps
%     image_padding              - mirror-pad image before processing
%     primary_confidence_level   - confidence level for immediate result 
%                                  acceptance, default 0.1 
%     secondary_confidence_level - confidence level for reinforced result
%                                  acceptance, default 0.25
%
% The response_map structure will contain the following fields:
%  - map_cor         - correlation field
%  - est_cor         - predicted correlation field
%  - map_prp         - normalized response map with real-valued scores [0,1]
%  - wnd_size        - utilized window size map
% and optionally:
%  - map_pva         - p-value based response map
%  - map_pce         - PCE response map
%
% References:
%
% [1] P. Korus, J. Huang, Multi-scale Analysis Strategies in PRNU-based 
%     Tampering Localization, Submitted to IEEE Tran. Information Forensics
%     and Security
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

    if ~exist('advanced', 'var')
        advanced = struct();
    end

    % Unpack advanced parameters
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

    if isfield(advanced, 'primary_confidence_level')
        primary_confidence_level = advanced.primary_confidence_level;
    else
        primary_confidence_level = 0.1;
    end       

    if isfield(advanced, 'secondary_confidence_level')
        secondary_confidence_level = advanced.secondary_confidence_level;
    else
        secondary_confidence_level = 0.25;
    end    
    
    if ischar(image)
        image = imread(image);
    end    
    
    if image_padding
        image = padarray(image, floor(129/2)*[1 1], 'symmetric', 'both');
    end
    
    % Start processing
    if verbose
      fprintf('Extracting image PRNU... '); tic;
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
        prnu = padarray(prnu, floor(129/2)*[1 1], 'symmetric', 'both');
    end    
    prnu = prnu .* image;

    if verbose
      t_last = toc; fprintf('done in %.2f min\n', t_last/60);
    end

    sum2 = @(x) sum(x(:));
                    
    if calc_pce
      map_list = {'map_cor', 'map_prp', 'est_cor', 'wnd_size', 'map_pce', 'map_pva'};
    else
      map_list = {'map_cor', 'map_prp', 'est_cor', 'wnd_size'};
    end
        
    % Allocate output structure
    response_map = struct();
    for map_name = map_list
        response_map.(map_name{1}) = nan(size(image,1), size(image,2));
    end
    
    if verbose
      fprintf('Detecting forgeries: '); tic;
    end

    bs = camera_model.predictors{1}.patch_size;
    
    bs2_min = floor(bs/2);
    bx = bs2_min + 1;
    fi = 1;
    
    maxbx = size(image,2) - bs2_min;
    maxby = size(image,1) - bs2_min;
            
    maxfx = ceil((size(image,2) - bs + 1)/bsk);
    maxfy = ceil((size(image,1) - bs + 1)/bsk);

    last_progress = -1;
    total = maxfx * maxfy;

    while bx <= maxbx
        by = bs2_min + 1;
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
            
            window_size_id = 1;
            lastScore = 0.5;
            
            while window_size_id <= numel(camera_model.predictors) && abs(lastScore - 0.5) < 0.5 - primary_confidence_level
                
                % Check if the current window can fit in this location
                bs = camera_model.predictors{window_size_id}.patch_size;
                bs2 = floor(bs/2);
                if bx - bs2 <= 0 || by - bs2 <= 0 || by + bs2 > size(image,1) || bx + bs2 > size(image,2)
                    window_size_id = numel(camera_model.predictors) + 1;
                    continue
                end
                
                predictor = camera_model.predictors{window_size_id};
                absence_model = camera_model.absence_models{window_size_id};
                
                % Extract patch features
                l_patch       = image((by-bs2):(by+bs2-1), (bx-bs2):(bx+bs2-1));
                std_patch     = std_image((by-bs2):(by+bs2-1), (bx-bs2):(bx+bs2-1));
                hpf_std_patch = hpf_std_image((by-bs2):(by+bs2-1), (bx-bs2):(bx+bs2-1));
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
                l_noise = noise((by-bs2):(by+bs2-1), (bx-bs2):(bx+bs2-1));
                l_prnu  =  prnu((by-bs2):(by+bs2-1), (bx-bs2):(bx+bs2-1));

                % Make sure the noise is 0-mean
                l_noise = l_noise - mean2(l_noise);
                l_prnu = l_prnu - mean2(l_prnu);

                % Compute normalization factors
                n1 = sqrt(sum2(l_noise .* l_noise));
                n2 = sqrt(sum2(l_prnu .* l_prnu));

                % Peak Correlation Energy
                if calc_pce
                    Corr = crosscorr(l_noise, l_prnu);
                    C = PCE(Corr);            
                    response_map.map_pce(by, bx) = C.PCE;
                    response_map.map_pva(by, bx) = C.pvalue;
                end

                % Normalized correlation
                ncorr = sum2(l_noise .* l_prnu) / n1 / n2;

                % Compute tampering probability
                mu_0 = 0;
                mu_1 = pred_cor;
                si_0 = absence_model.normfit(2);
                si_1 = predictor.normfit(2);
                
                pScore = corr2Response(ncorr, mu_0, si_0, mu_1, si_1);

                % If the new score is better than the previous one
                if abs(pScore - 0.5) > abs(lastScore - 0.5)
                    
                    % Update current value in the maps
                    response_map.map_prp(by, bx) = pScore;
                    response_map.wnd_size(by, bx) = window_size_id;
                    response_map.map_cor(by, bx) = ncorr;
                    response_map.est_cor(by, bx) = pred_cor;
                    
                    % If the new score reinforces a reasonably confident 
                    % old score, skip further processing                    
                    if abs(lastScore - 0.5) > 0.5 - secondary_confidence_level && (pScore - 0.5)*(lastScore-0.5) > 0 
                        window_size_id = numel(camera_model.predictors);
                    end
                    
                    lastScore = response_map.map_prp(by, bx);
                end
               
                % Increase window size
                window_size_id = window_size_id + 1;
            end
                        
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
            response_map.(map_name{1}) = response_map.(map_name{1})(mod(bs2_min-1,bsk)+2:bsk:end, mod(bs2_min-1,bsk)+2:bsk:end);
            response_map.(map_name{1}) = fillBorders(response_map.(map_name{1}));
            if image_padding
                margin = floor(129/2) / bsk;
                map = response_map.(map_name{1});
                map = map(margin:end-margin, margin:end-margin);                
                response_map.(map_name{1}) = map;
            end            
        end
    end
    
end