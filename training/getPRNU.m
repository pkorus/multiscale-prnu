function noise = getPRNU(image, sigma, sepChannels)

if ischar(image)
    image = imread(image);
end

image = double(image);

if ~exist('sepChannels', 'var')
    sepChannels = false;
end

% Extract basic noise
if size(image,3) == 1
    noise = extractNoise(image, sigma^2);
elseif size(image,3) == 3
    noise = zeros(size(image,1), size(image,2));
    for i = 1:3
        noise(:,:,i) = extractNoise(image(:,:,i), sigma^2);
    end
    if ~sepChannels
        noise = rgb2gray1(noise);
    end
else
    throw(MException('PRNU:InputDataError', 'Unsupported color mode'));
end

% Zero mean row and column-wise
noise = ZeroMeanTotal(noise);

% Wiener filter to remove other visually identifiable patterns
noise = WienerInDFT(noise,std2(noise));

end

function noise = extractNoise(image, sigma)
    [~, noise] = BM3D([], image, sigma);
    noise = image - 255*noise;
end
