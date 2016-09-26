function y = catChannels(x)
    [nRows, nCols, nChannels] = size(x);
    y = zeros(nRows, nCols * nChannels);
    for n = 1:nChannels
        y(:, 1+(n-1)*nCols:n*nCols) = x(:,:,n);
    end
end