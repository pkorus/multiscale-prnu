function [map, removed_mass, removed_ccs] = mapCleanup(map, cc_thresh)
% [map, removed_mass, removed_ccs] = mapCleanup(map, cc_thresh)
% 
% Simple tampering map cleanup heuristic - removes connected components
% smaller than given size. Optionally reports stats on the number of 
% removed pixels / objects.

    ccs = bwconncomp(map);

    % Remove components smaller than N blocks
    removed_ccs = 0;
    removed_mass = 0;
    for l = 1:ccs.NumObjects
        stm = numel(ccs.PixelIdxList{l});
        if stm < cc_thresh
            removed_mass = removed_mass + stm;            
            removed_ccs = removed_ccs + 1;
            map(ccs.PixelIdxList{l}) = 0;
        end
    end
    removed_ccs = removed_ccs / ccs.NumObjects;
    removed_mass = removed_mass / sum(map(:));
end