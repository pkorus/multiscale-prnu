function x = fillBorders(x, zeros)
% x = fillBorders(x, zeros)
% 
% Simple procedure that fills NaNs along image edges with the last known
% non-NaN values.
    
    if nargin < 2
        zeros = false;
    end

    to_fill = false(1,size(x,2));
    for i = 1:size(x,2);
        if zeros
            to_fill(i) = all(x(:,i) == 0);
        else
            to_fill(i) = all(isnan(x(:,i)));
        end
    end

    indices = find(to_fill);
    cands = [find(~to_fill, 1, 'first') find(~to_fill, 1, 'last')];
    for col = indices
        dist = abs(col - cands);
        [~, candidate] = min(dist);
        x(:, col) = x(:, cands(candidate));
    end


    to_fill = false(1,size(x,1));
    for i = 1:size(x,1);
        if zeros
            to_fill(i) = all(x(i,:) == 0);
        else
            to_fill(i) = all(isnan(x(i,:)));
        end
    end

    indices = find(to_fill);
    cands = [find(~to_fill, 1, 'first') find(~to_fill, 1, 'last')];
    for col = indices
        dist = abs(col - cands);
        [~, candidate] = min(dist);
        x(col, :) = x(cands(candidate), :);
    end

end