function detection_results = detectForgeryPRNUMultiscale(image, camera_model, scales, stride, detector, verbose)
% detection_results = detectForgeryPRNUMultiscale(image, camera_model, scales, overlap, detector, verbose)
% 
% Wrapper function that unifies various localization variants with a single API.
%
% Parameters:
%   - image         - input image, either uint8 RGB or string with image path
%   - camera_model  - structure with camera's PRNU model (see trainCameraModel)
%   - scales        - desired scales of analysis, e.g., [65, 97, 129]
%   - stride        - sliding window stride, fractional values can be used to
%                     use stride relative to window size, e.g.:
%                       >=1  - use directly as stride
%                       0.5 - stride = half window size
%                       0   - no overlap, windows next to each other
%                     Note:
%                      - values will be floored to the powers of 2
%                      - for FWA analysis, minimal stride is 8
%                        (localization resolution)
%   - detector      - analysis strategy:
%                     'standard'  - standard implementation with 
%                                   central-pixel attribution (CPA)
%                     'blockwise' - alternative implementation with
%                                   full-window attribution (FWA)
%                     'central'   - faster implementation based on image filtering, should be roughly equivalent to 'standard'
%                     'boxcar'    - same as 'central' but the correlation field will undergo an additional boxcar filtering
%                     'guided'    - same as 'central' but used guided filtering instead of a box filter
%   - verbose       - set to true to display a progress bar
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

    if nargin < 6 || isempty(verbose)
        verbose = true;
    end

    maps = cell(1, numel(scales));    
    if ischar(image)    
        image = imread(image);
    end
    parfor si = 1:numel(scales)
        s = scales(si);
        if stride >= 1
            local_stride = round(stride);
        else
            local_stride = max(1, power(2,floor(log2((1-stride)*s))));
        end
        if strfind('blockwise', detector) == 1        
            maps{si} = detectForgeryPRNU(image, camera_model, s, max(8, local_stride), verbose);
        elseif strfind('boxcar', detector) == 1
            maps{si} = detectForgeryPRNUFilter(image, camera_model, s, struct('stride', local_stride, 'rates', 0.01, 'mode', 'b', 'verbose', verbose));
        elseif strfind('guided', detector) == 1
            maps{si} = detectForgeryPRNUFilter(image, camera_model, s, struct('stride', local_stride, 'rates', 0.01, 'mode', 'g', 'verbose', verbose));
        elseif strfind('central', detector) == 1
            maps{si} = detectForgeryPRNUFilter(image, camera_model, s, struct('stride', local_stride, 'rates', 0.01, 'mode', 'c', 'verbose', verbose));
        elseif strfind('standard', detector) == 1
            maps{si} = detectForgeryPRNUCentral(image, camera_model, s, struct('stride', local_stride, 'rates', 0.01, 'verbose', verbose, 'calc_pce', true));
        else
            exception = MException('prnu:InvalidParameter', sprintf('ERROR Unknown detector (%s)!', detector));
            throw(exception);
        end
            
    end
    
    detection_results = struct(); 
    for si = 1:numel(scales)
        if isempty(maps{si})
            continue;
        end        
        for name = fieldnames(maps{1})'
            detection_results.(name{1}){si} = maps{si}.(name{1});
        end        
    end
    detection_results.block = scales;
    
end

