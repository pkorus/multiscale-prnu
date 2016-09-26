function [edgeStruct, adj] = createNeighborhood(nRows, nCols, nMode)
% [edgeStruct, adj] = createNeighborhood(nRows, nCols, nMode)
% 
% Creates an adjacency structure for a grid-neighborhood system.
% 
% Parameters:
%   - nRows, nCols    - image size
%   - nMode           - neighborhood mode, either 4 or 8
%                       4 - includes top,down,left,right neighbors
%                       8 - same as 4 + includes diagonal neighbors
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

    if ~exist('nMode', 'var') 
        nMode = 4;
    end

    if nMode ~= 4 && nMode ~= 8
        throw(MException('crf:InvalidParameter', 'Invalid neighborhood mode (4 and 8 are supported)!'));
    end

    % Construct the Undirected Graphical Model
    nNodes = nRows*nCols;
    nStates = 2;

    cache_name = sprintf('cache/mrf_edge_struct_cache_%d_%d_%d.mat', nRows, nCols, nMode);

    if ~exist(cache_name, 'file')
        adj = sparse(nNodes,nNodes);
        
        dist = 1;

        exclude_lr = zeros(dist, nCols);
        exclude_fr = zeros(dist, nCols);
        exclude_lc = zeros(dist, nRows);

        for d = 1:dist
            % No Down  edge for last  row
            exclude_lr(d,:) = sub2ind([nRows nCols],repmat(nRows - d + 1,[1 nCols]),1:nCols);
            % No Up    edge for first row
            exclude_fr(d,:) = sub2ind([nRows nCols],repmat(d,[1 nCols]),1:nCols);
            % No Right edge for last  column
            exclude_lc(d,:) = sub2ind([nRows nCols],1:nRows,repmat(nCols - d + 1,[1 nRows]));
        end

        % Add Down Edges
        ind = 1:nNodes;
        ind = setdiff(ind,exclude_lr);
        adj(sub2ind([nNodes nNodes],ind,ind+dist)) = 1;

        % Add Right Edges
        ind = 1:nNodes;
        ind = setdiff(ind,exclude_lc);
        adj(sub2ind([nNodes nNodes],ind,ind+nRows*dist)) = 1;

        if nMode == 8
            % Add Diagonal Down Edges
            ind = 1:nNodes;
            ind = setdiff(ind,exclude_lc);
            ind = setdiff(ind,exclude_lr);
            adj(sub2ind([nNodes nNodes],ind,ind+nRows*dist+dist)) = 1;    

            % Add Diagonal Up Edges
            ind = 1:nNodes;
            ind = setdiff(ind,exclude_lc);
            ind = setdiff(ind,exclude_fr);
            adj(sub2ind([nNodes nNodes],ind,ind+nRows*dist-dist)) = 1;    
        end

        % Add Up/Left Edges
        adj = adj+adj';
        edgeStruct = UGM_makeEdgeStruct(adj,nStates);
        save(cache_name, 'adj', 'edgeStruct');
    else
        load(cache_name, 'adj', 'edgeStruct');
    end

end