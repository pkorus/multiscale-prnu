function maps = extractCellMaps(candidate_maps, fieldname, selector)
% maps = extractCellMaps(candidate_maps, fieldname, selector)
%
% Extracts tampering localization maps from a cell array and converts
% them into a 3-D array. Successive candidate maps are stacked along
% the 3rd dimension. 

    if isstruct(candidate_maps{1})
        % If the input is a cell array of structures

        if ~exist('fieldname', 'var')
            if isfield(candidate_maps{1}, 'candidate')
                fieldname = 'candidate';
            else
                all_fields = fieldnames(candidate_maps{1});
                is_array = false(1, numel(all_fields));
                for i = 1:numel(all_fields);
                    temp = candidate_maps{1}.(all_fields{i});
                    if size(temp,1) > 1 && size(temp,2) > 1
                        is_array(i) = true;
                    end
                end
                fieldname = all_fields(1);
            end            
        end

        first_map = candidate_maps{1}.(fieldname);
        maps = zeros(size(first_map,1), size(first_map, 2), numel(candidate_maps));

        if ~exist('selector', 'var')
            selector = 1:numel(candidate_maps);
        else
            if islogical(selector)
                selector = find(selector);
            end
            if any(selector > numel(candidate_maps)) || any(selector < 1)
                throw(MException('prnu:InvalidParameter', 'Unsupported or invalid map selector!'));
            end
        end

        for i = selector
            temporary = candidate_maps{i}.(fieldname);
            if size(first_map,1) ~= size(temporary, 1) || size(first_map,2) ~= size(temporary, 2)
                temporary = imresize(temporary, [size(first_map,1), size(first_map, 2)], 'bilinear');
            end
            maps(:,:,i) = temporary;
        end
    
    else
        % If the input is a cell array of arrays
        [nRows, nCols, nDim] = size(candidate_maps{1});
        maps = zeros(nRows, nCols, numel(candidate_maps));
        for n = 1:numel(candidate_maps)
            temp = candidate_maps{n};
            if nDim > 1
                temp = rgb2gray(temp);
            end
            if size(temp,1) ~= nRows || size(temp, 2) ~= nCols
                temp = imresize(temp, [nRows, nCols]);
            end
            maps(:,:,n) = temp;
        end        
    end
end