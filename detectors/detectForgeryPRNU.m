function response_map = detectForgeryPRNU(image, camera_model, bs, bsk, verbose, calc_pce)
% response_map = detectForgeryPRNU(image, camera_model, bs, bsk, verbose, calc_pce)
%
% Tampering localization based on PRNU analysis. Analysis is performed with full window
% attribution. For central pixel attribution, see detectForgeryPRNUCentral.
% 
% Params: 
%  - image           - RGB image, uint8
%  - camera_model    - camera model structure
%  - bs              - block size, e.g., 128 (only some values are allowed)
%  - bsk             - stride, e.g., 8
%  - verbose         - boolean
%  - calc_pce        - include also PCE and p-value
%
% The response_map structure will contain the following fields:
%  - map_cor         - correlation field
%  - est_cor         - predicted correlation field
%  - map_prp         - normalized response map with real-valued scores [0,1]
% and optionally:
%  - map_pva         - p-value based response map
%  - map_pce         - PCE response map
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
    if ~ismember(bs, [32 48 64 96 128 192 256])
        error('prnu:InvalidParameter', 'Invalid analysis window size, accepted values include [32 48 64 96 128 192 256].');
    end

    bs = bs / 8;
    bsk = bsk / 8;

    if mod(bs, bsk) ~= 0
        error('prnu:InvalidParameter', 'Unsupported block skip - needs to divide block size');
    end

    if bsk > bs
        error('prnu:InvalidParameter', 'Block skip cannot exceed block size');
    end
    
    if nargin < 5 || isempty(verbose)
        verbose = false;
    end

    if nargin < 6 || isempty(calc_pce)
        calc_pce = false;
    end

    % ---------------------------------------------------------------------

    % Helper variables
    maxbx = floor(size(image,2)/8);
    maxby = floor(size(image,1)/8);
    
    maxfx = ceil((maxbx - bs + 1)/bsk);
    maxfy = ceil((maxby - bs + 1)/bsk);

    % Optional other maps 'map_plm', 'map_cnm', 'map_cso'    
    if calc_pce
        field_mapping = {'map_cor', 'est_cor', 'map_prp', 'map_pce', 'map_pva'};    
    else
        field_mapping = {'map_cor', 'est_cor', 'map_prp'};
    end

    scores = zeros(maxfx*maxfy, numel(field_mapping));
    
    if verbose
        fprintf('Running (%s, %d px) - ', camera_model.name, bs*8);
        fprintf('In total %d windows to process.\n', size(scores,1));        
        fprintf('  Extracting image PRNU... '); tic;
    end

    % Extract PRNU from the image
    if strcmp(camera_model.stats.denoising, 'bm3d')
        noise = getPRNU(image, camera_model.stats.noise_sigma);
    else
        noise = NoiseExtractFromImage(image, camera_model.stats.noise_sigma);
    end
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
        throw(MException('PRNU:BadParameter', 'Unsupported feature set in the camera model!'))
    end
    
    if verbose
       t_last = toc; fprintf('done in %.2f min\n', t_last/60);
    end
    
    % Fetch camera's PRNU
    image = double(rgb2gray(image));
    std_image = stdfilt3(image);
    hpf_std_image = stdfilt3(highPass(image));
    prnu = rgb2gray1(camera_model.prnu);
    prnu = WienerInDFT(prnu,std2(prnu));
    prnu = prnu .* image;
    
    % Check if the camera has an appropriate predictor
    predictor = [];
    using_predictor = false;
    if isfield(camera_model, 'predictors')
        for p = 1:numel(camera_model.predictors)
            if camera_model.predictors{p}.patch_size == bs*8
                predictor = camera_model.predictors{p};
                using_predictor = true;
            end
        end
    end
    if ~using_predictor
        fprintf('WARNING Missing predictor (%s, %d px)\n', camera_model.name, bs*8);
    end
    
    % Check if the camera has an appropriate null model
    absence_model = [];
    if isfield(camera_model, 'absence_models')
        for p = 1:numel(camera_model.absence_models)
            if camera_model.absence_models{p}.patch_size == bs*8
                absence_model = camera_model.absence_models{p};
            end
        end
    end
    if isempty(absence_model)
        fprintf('WARNING Missing PRNU absence model (%s, %d px)\n', camera_model.name, bs*8);
        absence_model = struct();
        absence_model.normfit = [0 0.01];
    end
    
    sum2 = @(x) sum(x(:));
    
    if verbose
        fprintf('  Detecting forgeries'); tic;
    end
    bx = 1;
    fi = 1;
    last_progress = -1;
    total = maxfx * maxfy;

    while bx <= maxbx - bs + 1
        by = 1;
        while by <= maxby - bs + 1
            
            % Simple progress display
            if verbose
                counter = fi;
                progress = round(100*(counter/total)/2);
                % Display progress bar
                if progress ~= last_progress
                    displayProgress(sprintf('Window   %6d / %6d', fi, total), progress, 50);
                end
                last_progress = progress;
            end
            
            % Run predictor
            l_patch       = image((by-1)*8+1:(by+bs-1)*8, (bx-1)*8+1:(bx+bs-1)*8);               
            std_patch     = std_image((by-1)*8+1:(by+bs-1)*8, (bx-1)*8+1:(bx+bs-1)*8);
            hpf_std_patch = hpf_std_image((by-1)*8+1:(by+bs-1)*8, (bx-1)*8+1:(bx+bs-1)*8);
            if feature_set == 0
                features = extractPatchFeaturesFast(l_patch, std_patch, hpf_std_patch);
            else
                features = [extractPatchFeaturesFast(l_patch, std_patch, hpf_std_patch) jpeg_feature];
            end            
            pred_cor = features * predictor.correlation;
            
            % Extract local patches for verification
            l_noise = noise((by-1)*8+1:(by+bs-1)*8, (bx-1)*8+1:(bx+bs-1)*8);
            l_prnu  = prnu((by-1)*8+1:(by+bs-1)*8, (bx-1)*8+1:(bx+bs-1)*8);

            % Make sure they are 0-mean
            l_noise = l_noise - mean2(l_noise);
            l_prnu  = l_prnu - mean2(l_prnu);
            
            % Compute normalization
            n1 = sqrt(sum2(l_noise .* l_noise));
            n2 = sqrt(sum2(l_prnu .* l_prnu));
            
            % Prevent numerical instabiliy
            if n2 == 0; n2 = 1e-5; end            
            if n1 == 0; n1 = 1e-5; end
            
            % Traditional correlation field -------------------------------
            scores(fi,1) = sum2(l_noise .* l_prnu) / n1 / n2;

            % Just the predictors -----------------------------------------
            scores(fi,2) = pred_cor;            

            % Remapped response maps (to [0,1]) ---------------------------
            mu_0 = 0;
            mu_1 = pred_cor;
            si_0 = absence_model.normfit(2);
            si_1 = predictor.normfit(2);
            
            scores(fi,3) = corr2Response(scores(fi,1), mu_0, si_0, mu_1, si_1);

            if calc_pce
                Corr = crosscorr(l_noise, l_prnu);
                C = PCE(Corr);
                scores(fi,4) = C.PCE / pred_pce;
                scores(fi,5) = C.pvalue;
            end

            % -------------------------------------------------------------
            
            % Move window
            by = by + bsk;
            fi = fi + 1;
        end
        bx = bx + bsk;
    end

    if verbose
        t_last = toc; fprintf(' done in %.2f min\n', t_last/60);    
        fprintf('  Fusing candidate scores...'); tic;
    end
    
    % Allocate output structure
    fused_maps = zeros(maxby, maxbx, numel(field_mapping));
    
    ptypes = ceil(bs/bsk)*ceil(bs/bsk);
    prototype = zeros(ptypes,2);
    si = 1;
    for s1 = 0:bsk:bs-1
        for s2 = 0:bsk:bs-1
            prototype(si,:) = [s1 s2];
            si = si + 1;
        end
    end
    for bx = 1:maxbx
        for by = 1:maxby
            % Determine indices of relevant examples - for overlapping windows
            cands = repmat([bx by], ptypes, 1) - prototype;
            cands = ceil(cands/bsk);            
            v = zeros(1,numel(field_mapping));
            n = 0;
            for si = 1:ptypes
                if cands(si,1) > 0 && cands(si,2) > 0 && cands(si,1) <= maxfx && cands(si,2) <= maxfy
                    fi = (cands(si,1)-1)*maxfy + cands(si,2);
                    for vali = 1:numel(field_mapping)
                        v(vali) = v(vali) + scores(fi,vali);
                    end
                    n = n + 1;
                end
            end
            for vali = 1:numel(field_mapping)
                fused_maps(by,bx,vali) = v(vali)/n;
            end
        end
    end
    
    if verbose
        t_last = toc; fprintf('done in %.2f s\n', t_last);
    end

    response_map = struct();
    for vali = 1:numel(field_mapping)
        response_map.(field_mapping{vali}) = fused_maps(:,:,vali);
    end
    
    if bsk > 1 && any(any(isnan(response_map.map_cor)))
        for vali = 1:numel(field_mapping)
            response_map.(field_mapping{vali}) = fillMissing(response_map.(field_mapping{vali}));
        end
    end    

end