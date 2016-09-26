function [decision, indices] = fuseAdaptWnd(response_maps, primary_conf_level, secondary_conf_level)
% [decision, indices] = fuseAdaptWnd(response_maps, primary_conf_level, secondary_conf_level)
%
% Adaptive-window localization strategy implemented as a fusion of
% multi-scale candidate tampering maps. 
%
%  - primary_conf_level   - confidence level for immediate acceptance of
%                           small scale score, defaults to 0.1
%
%  - secondary_conf_level - confidence level for accepting a candidate
%                           score, if the next-scale candidate agrees with 
%                           the decision
%
% References:
%
% [1] P. Korus, J. Huang, Multi-scale Analysis Strategies in PRNU-based 
%     Tampering Localization, Submitted to IEEE Tran. Information Forensics
%     and Security
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

    if ~exist('primary_conf_level', 'var')
        primary_conf_level = 0.1;
    end

    if ~exist('secondary_conf_level', 'var')
        secondary_conf_level = 0.25;
    end    
    
    response = extractCellMaps(response_maps);
    [nRows, nCols, nMaps] = size(response);

    indices = zeros(nRows, nCols);
    decision = zeros(nRows, nCols);

    for i = 1:nRows
        for j = 1:nCols

            index = 2;        
            accept_index = 1;

            while index <= nMaps && abs(response(i, j, accept_index) - 0.5) < 0.5 - primary_conf_level

                if abs(response(i, j, index) - 0.5) > abs(response(i, j, accept_index) - 0.5)

                    accept_index = index;
                    if secondary_conf_level > 0 && abs(response(i, j, index - 1) - 0.5) > 0.5 - secondary_conf_level && (response(i, j, index - 1) - 0.5)*(response(i, j, index) - 0.5) > 0
                        break;
                    end                        
                end

                index = index + 1;                                        
            end

            indices(i, j) = accept_index;        
        end
    end

    for i = 1:nRows
        for j = 1:nCols
            index = indices(i, j);
            decision(i, j) = response_maps{index}.candidate(i, j);
        end
    end

end