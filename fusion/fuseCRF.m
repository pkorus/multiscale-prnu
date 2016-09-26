function [MAP] = fuseCRF(response_map, image, thresh, w, o)
% [MAP] = fuseCRF(response_map, image, thresh, w, o)
%
% Fusion of candidate response maps based on conditional random fields.
% Initially developed for multi-scale fusion, but can also be used for 
% multi-modal fusion. For more information, see [1,2].
%
% Input parameters:
% 
%  - response_map      - cell array with input maps: structures with fields
%                        'candidate' and 'reliability' for tampering
%                        probability maps and their corresponding
%                        reliability information.
%
%  - image             - the tampered image to use for content guidance;
%                        (will be downsampled if necessary), leave empty if
%                        not needed
%
%   - thresh           - initial quasi-threshold, default 0.5
%
%   - w                - parameters of the conditional random field (see
%                        [1,2] for more details). 5-D vector:
%                          w(1) - alpha (decision bias)
%                          w(2) - base interaction strength
%                          w(3) - content-dependent interaction strength
%                          w(4) - pixel similarity attenuation
%                          w(5) - threshold drift
%                        The parameter w(4) may need to be adjusted
%                        depending on whether your image is in (0,255) or
%                        (0,1). In the former case, values around 25
%                        typically lead to good performance. 
%
%   - o                - advanced parameters (see defaultFusionSettings):
%
%     neighborhood_mode - 4/8 - use either 4 or 8 nearest neighbours.
%
%     min_potential    - potential clipping threshold, prevents infinite energies.
%
%     discard_unreliable_maps - whether to reject candidate maps based on their
%                        expected utility:
%                        0 - use all maps,
%                        1 - discard empty maps,
%                        2 - discard maps similar to Gaussian noise.
%
%     cand_map_filter_threshold - a threshold for the above candidate map filtering
%
%     unreliable_score - candidate score values that will be used inside of
%                        unreliable image regions: values close to 0 will prefer
%                        labelling such regions as authentic; values around the
%                        decision threshold will facilitate easier propagation
%                        of decisions from the areas' neighbourhood.
%
%     unreliable_score_strength - can be used to partially retain original scores
%                        in unreliable areas (by simple weighted average).
%
%     threshold_saturation_gap - distance from 0/1 below/above which it is no
%                        longer possible to shift the quasi-threshold.
%
% References:
%
% [1] P. Korus, J. Huang, Multi-scale Analysis Strategies in PRNU-based 
%     Tampering Localization, Submitted to IEEE Tran. Information Forensics
%     and Security
%
% [2] P. Korus, J. Huang, Multi-Scale Fusion for Improved Localization of 
%     Malicious Tampering in Digital Images, IEEE Tran. Image Processing,
%     March 2016, 10.1109/TIP.2016.2518870
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

    % Read candidate maps from the provided structure
    if ~isfloat(response_map) && ~iscell(response_map)
        response_map = double(response_map);
    end
    
    if ~iscell(response_map)        
        if size(response_map,3) == 1        
            cell_map = cell(1);
            cell_map{1}.candidate = response_map;
            cell_map{1}.reliability = ones(size(cell_map{1}.candidate));
        else
            cell_map = cell(1, size(response_map,3));
            for n = 1:numel(cell_map)
                cell_map{n}.candidate = response_map(:,:,n);
                cell_map{n}.reliability = ones(size(cell_map{n}.candidate));
            end
        end
        
        response_map = cell_map;
    end
    
    [nRows, nCols] = size(response_map{1}.candidate);

    if ~exist('w', 'var') || isempty(w)
        w = 0;
    end    
    
    % Pad with zeros, if unspecified
    if numel(w) < 5
        w = [w zeros(1, 5 - numel(w))];
    end       
    
    if ~exist('o', 'var') || isempty(o)
        o = defaultFusionSettings();
    end
    
    if ~exist('thresh', 'var') || isempty(thresh)
        thresh = 0.5;
    end
    
    if isscalar(o.min_potential)
        o.min_potential = o.min_potential * [1 1];
    end
    
    % Pre-process the tampered image (if given)
    if ~isempty(image)
        if ischar(image)
            image = imread(image);
        end

        if size(image,1) ~= nRows || size(image,2) ~= nCols
            image = imresize(image, [nRows, nCols]);
        end

        % Reshape the input image for easier access
        image = double(image);
        image_v = reshape(image, [nCols*nRows, size(image,3)]);        
    end
      
    % Remove useless candidate maps
    if o.discard_unreliable_maps > 0
        
        if o.discard_unreliable_maps == 1
            % Empty map removal
            [valid_indices, map_distances] = filterEmptyMaps(response_map, o.cand_map_filter_threshold);
        elseif o.discard_unreliable_maps == 2            
            % Remove maps based on similarity to Gaussian noise
            [valid_indices, map_distances] = filterValidMaps(response_map, o.cand_map_filter_threshold);
        end               
                        
        % If there are no valid maps, use the best available one
        if all(~valid_indices)
            [~, best_map_ind] = max(map_distances);
            valid_indices(best_map_ind) = true;
        end
        
        response_map = response_map(valid_indices);
    end
        
    % Construct the undirected graphical model
    nStates = 2;    
    nNodes = nRows * nCols;        
    edgeStruct = createNeighborhood(nRows, nCols, o.neighborhood_mode);
    
    % Remember per-node thresholds (for threshold drift)
    thresholds = thresh * ones(numel(response_map{1}.candidate), 1);
    
    % Initialize node potentials
    nodePot = ones(nNodes,nStates);
    
    % Setup resetting of unreliable scores
    % If the score is negative, set it to the current threshold (mul. by
    % abs value of the unreliable score).
    if o.unreliable_score < 0
        o.unreliable_score = thresh * abs(o.unreliable_score);
    end

    for i = 1:numel(response_map)
        
        Y = response_map{i}.reliability(:);
        X = response_map{i}.candidate(:);
                
        % Set a predefined value in unreliable regions.
        % If this is not desirable, setting o.unreliable_score_strength = 0 will
        % disable this step.
        X = X .* Y + (o.unreliable_score*o.unreliable_score_strength + (1-o.unreliable_score_strength) .* X) .* (1-Y);

        % Threshold drift
        if any(w(5)) > 0
            if i == 1
                thresholds = thresh * ones(size(X));
            else
                sel = lastX > last_thresholds;                
                
                % Shift the thresholds according to current "decisions"
                new_thresholds = thresholds;
                
                % Small scale has not been detected
                new_thresholds(~sel) = min(1 - o.threshold_saturation_gap, thresholds(~sel) + w(5));
                
                % Small scale has been detected
                new_thresholds(sel) = max(o.threshold_saturation_gap, thresholds(sel) - w(5));
                
                % Apply the new thresholds only in reliable areas (soft)
                thresholds = new_thresholds .* Y + thresholds .* (1 - Y);                
            end

            lastX = X;
            last_thresholds = thresholds;
        end
                
        potFun1 = @(X) max(o.min_potential(1), 1 - (1./2./thresholds) .* X);
        potFun2 = @(X) max(o.min_potential(2), (1./2./(1-thresholds)) .* X + 1 - (1./2./(1-thresholds)));
        
        nodePot(:,1) = nodePot(:,1) .* potFun1(X);
        nodePot(:,2) = nodePot(:,2) .* potFun2(X);
    end

    % Normalize by the number of candidate maps
    nodePot = nodePot .^ (1/numel(response_map));
    
    % Penalty for tampered regions
    nodePot(:,2) = nodePot(:,2) * exp(-w(1));
    
    % Setup pairwise potentials
    if ~isempty(image) && (w(3) > 0 || nargout > 1)
        % Compute adaptive interactions
        edgePot = zeros(nStates, nStates, edgeStruct.nEdges);
        l2_dist_cache = zeros(1, edgeStruct.nEdges);
        fade_strength = 0.5 / (w(4)^2);
        
        for ei = 1:edgeStruct.nEdges
            nodes = edgeStruct.edgeEnds(ei,:);
            n1 = nodes(1);
            n2 = nodes(2);
            col_distance = image_v(n1,:) - image_v(n2,:);
            l2_dist_cache(ei) = sum(col_distance.^2);
            beta_local = w(2) + w(3) * exp(-fade_strength*(l2_dist_cache(ei)));
            A = [0 beta_local; beta_local 0];
            edgePot(:,:,ei) = exp(-A);
        end
        
    else
        % If not possible, use only default interactions
        A = [0 w(2); w(2) 0];
        A = exp(-A);
        edgePot = repmat(A,[1 1 edgeStruct.nEdges]);
    end
    
    % Decoding with Graph Cuts
    optimalDecoding = UGM_Decode_GraphCut(nodePot,edgePot,edgeStruct);
    
    % Convert the solution to a tampering map
    MAP = im2bw(double(reshape(optimalDecoding,nRows,nCols))-1);    
end
 
function [valid_indices, distance] = filterEmptyMaps(cand_maps, th)
% Remove nearly empty candidate maps
    
    if nargin < 2
      th = 0.1;
    end

    distance = zeros(1, numel(cand_maps));
    
    for i = 1:numel(cand_maps)
        map = cand_maps{i}.candidate;
        mx = quantile(map(:), 0.99);
        distance(i) = mx;
    end
    
    valid_indices = distance > th;
end

function [valid_indices, distance] = filterValidMaps(cand_maps, threshold, histogram_bins)
% Removes unreliable candidate maps by thresholding the KL distance
% between the empirical distribution and a MLE-based Gaussian fit.

    if nargin < 2
        threshold = 0.1;
    end
    if nargin < 3
        histogram_bins = 25;
    end

    distance = zeros(1, numel(cand_maps));
    
    if isfield(cand_maps{1}, 'reliability')
        reliability_field = 'reliability';
    elseif isfield(cand_maps{1}, 'saturation')
        reliability_field = 'saturation';
    else
        reliability_field = [];
    end
    
    for i = 1:numel(cand_maps)
        
        data = cand_maps{i}.candidate;
        if ~isempty(reliability_field)
            satp = cand_maps{i}.(reliability_field);
        else
            satp = ones(size(data));
        end
        data(satp < 0.025) = NaN;
        data = data(~isnan(data));
        
        if numel(data) > 64
            [pi, xi] = hist(data,linspace(0,1,histogram_bins));
            pi = pi/(sum(pi)*(xi(2)-xi(1)));
            [u, s] = normfit(data);
            ni = normpdf(xi, u, s);
            
            if (u < 0.05 || u > 0.95) && s < 0.05
                distance(i) = -1;
            elseif abs(u - 0.5) < 0.05 && s < 0.04
                distance(i) = -2;
            else
                distance(i) = kldiv(pi, ni);
            end
        end
    end

    valid_indices = distance > threshold;
end