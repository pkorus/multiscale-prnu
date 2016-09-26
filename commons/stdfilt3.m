function J = stdfilt3(I)
% Matlab's stdfilt2 stripped down to its bare essentials (for faster processing).

    h = [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1];
    n = 25;

    if (~isa(I,'double'))
        I = double(I);
    end

    % J = algstdfilt2(I, h);
    n1 = n - 1;
    if exist('imfilter_mex') == 3
        conv1 = imfilter3(I.^2, h/n1);
        conv2 = imfilter3(I, h).^2 / (n*n1);
    else
        conv1 = imfilter(I.^2, h/n1);
        conv2 = imfilter(I, h).^2 / (n*n1);        
    end
    J = sqrt(max((conv1 - conv2),0));

end

function b = imfilter3(a, h)
    boundary = 'symmetric';
    sameSize = 1;
    convMode = 0;

    [finalSize, pad] = computeSizes(a, h, sameSize);

    % Zero-pad input based on dimensions of filter kernel.
    a = padarray(a, pad, boundary, 'both');
    b = filterPartOrWhole(a, finalSize, h, pad, sameSize, convMode);
end

%--------------------------------------------------------------
function ippFlag = useIPPL(a,outSize,h,nonzero_h)

    %We are disabling the use of IPP for double precision inputs on win32.
    if strcmp(computer('arch'),'win32')
        supportedTypes = {'uint8','uint16','int16','single'};
    else
        supportedTypes = {'uint8','uint16','int16','single','double'};
    end

    if ~any(strcmp(class(a),supportedTypes))
        ippFlag = false;
        return;
    end

    if numel(nonzero_h)/numel(h) > 0.05
        densityFlag = true;
    else
        densityFlag = false;
    end

    hDimsFlag = ismatrix(h);

    % Determine if the image is big depending on datatype
    if (isfloat(a))      
        tooBig = any(outSize>32750);
    else
        tooBig = any(outSize>65500);
    end

    ippFlag = densityFlag && hDimsFlag && (~tooBig);
end
%--------------------------------------------------------------
function [finalSize, pad] = computeSizes(a, h, sameSize)

    rank_a = ndims(a);
    rank_h = ndims(h);

    % Pad dimensions with ones if filter and image rank are different
    size_h = [size(h) ones(1,rank_a-rank_h)];
    size_a = [size(a) ones(1,rank_h-rank_a)];

    if (sameSize)
      %Same output
      finalSize = size_a;

      %Calculate the number of pad pixels
      filter_center = floor((size_h + 1)/2);
      pad = size_h - filter_center;
    else
      %Full output
      finalSize = size_a+size_h-1;
      pad = size_h - 1;
    end                        
end
%--------------------------------------------------------------
function result = filterPartOrWhole(a, outSize, h, start, sameSize, convMode)

% Create connectivity matrix.  Only use nonzero values of the filter.
    conn = h~=0;
    nonzero_h = h(conn);

    ippFlag  = useIPPL(a,outSize,h,nonzero_h);

    result = imfilter_mex(a, outSize, h, nonzero_h,...
                        conn, start, sameSize, convMode, ippFlag);
end
