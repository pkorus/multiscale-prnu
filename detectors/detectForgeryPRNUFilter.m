function guided_response = detectForgeryPRNUFilter(image, camera_model, bs, advanced)
% guided_response = detectForgeryPRNUFilter(image, camera_model, bs)
%
% Detect forgeries by means of PRNU analysis. The sliding window approach
% is reformulated as an image filtering problem, either conventional, or
% guided by the image content.
%
% Params:
%  - image           - RGB image, uint8
%  - camera_model    - camera model structure
%  - bs              - block size, e.g., 128 (only some values are allowed)
%  - advanced        - advanced parameters:
%     mode                 - filtering mode: 'central', 'boxcar' or 'guided'
%     stride               - stride for probability computation, e.g., 8
%     predictor_stride     - stride for correlation estimation, e.g., 8
%     cfar_std_mul         - tampering acceptance threshold (number of standard deviations), default 2
%     cfar_fp_rate         - false alarm rate for CFAR map, default 0.01
%     verbose              - boolean
%
% Output structure:
% The response_map structure will contain the following fields:
%  - map_cor         - correlation field
%  - est_cor         - predicted correlation field
%  - map_cfar        - binary decisions according to two-threshold CFAR rule
%  - bsk             - stride for corr field estimation, e.g., 8
%  - rates           - false alarm rate for CFAR detection
%  - mode            - 'central', 'boxcar', or 'guided'
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: October 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------


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

    if isfield(advanced, 'mode')
        mode = advanced.mode;
    else
        mode = 'central';
    end
    
    if isfield(advanced, 'stride')
        bsk = advanced.stride;
    else
        bsk = 4;
    end

    if isfield(advanced, 'predictor_stride')
        pbsk = advanced.predictor_stride;
    else
        pbsk = 8;
    end

    % Extract noise from the image
    if strcmp(camera_model.stats.denoising, 'bm3d')
        noise = getPRNU(image, camera_model.stats.noise_sigma);
    else
        noise = NoiseExtractFromImage(image, camera_model.stats.noise_sigma);
    end
    noise = WienerInDFT(noise, std2(noise));
    
    % Fetch camera's PRNU
    image = double(rgb2gray(image));
    prnu = rgb2gray1(camera_model.prnu);
    prnu = WienerInDFT(prnu,std2(prnu));
    prnu = prnu .* image;

    % Make sure the noises are 0-mean
    noise = noise - mean2(noise);
    prnu  = prnu  - mean2(prnu);

    sum2 = @(x) sum(x(:));
    
    % Normalization factors, compensated for the analysis window size
    n1 = sqrt(sum2(noise .* noise)); n1 = n1 * sqrt(bs*bs/numel(noise));
    n2 = sqrt(sum2(prnu .* prnu));   n2 = n2 * sqrt(bs*bs/numel(noise));

    filter_me = noise .* prnu;    
    
    matched_mode = 'central';

    % Peform actual filtering
    if strfind('guided', mode) == 1
        matched_mode = 'guided';
        response = imguidedfilter(bs * bs * filter_me, (image)/255, 'NeighborhoodSize', bs, 'DegreeOfSmoothing', 0.05);
    else
        response = filter2(ones(bs, bs), filter_me);
        response = response ./ filter2(ones(bs, bs), ones(size(filter_me))) * bs * bs; % Compensate for zero-padding at the borders    

        % If requested, the response map can be post-processed by an averaging filter
        if strfind('boxcar', mode) == 1
            matched_mode = 'boxcar';
            h = 1/(bs*bs)*ones(bs);
            response = imfilter(response, h, 'symmetric');        
        end
    end
    
    % Estimate the correlation field and apply predictor
    [decision_map, map_prp, est_field, map_cor] = applyPRNUCorrelationFieldCentral(image, camera_model, response, bs, [pbsk bsk], [cfar_std_mul cfar_fp_rate], verbose);

    % Populate output struct
    guided_response = struct();
    guided_response.window_size = bs;
    guided_response.filter_mode = matched_mode;
    guided_response.stride = bsk;
    guided_response.cfar_threshold_false_alarm = cfar_fp_rate;
    guided_response.cfar_threshold_multiplier = cfar_std_mul;
    guided_response.map_cfar = decision_map;
    guided_response.map_cor = map_cor;
    guided_response.map_prp = map_prp;
    guided_response.est_cor = est_field;
    
end