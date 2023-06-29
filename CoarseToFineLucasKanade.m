function [u, v] = CoarseToFineLucasKanade(im1, im2, windowSize, numLevels)

    % Create image pyramid
    pyramid1 = createPyramid(im1, numLevels);
    pyramid2 = createPyramid(im2, numLevels);
    
    % Initialize optical flow
    u = zeros(size(im1));
    v = zeros(size(im1));
    
    % For each level of the pyramid
    for level = numLevels:-1:1
        
        % Get images at current pyramid level
        currIm1 = pyramid1{level};
        currIm2 = pyramid2{level};
        
        % Resize the flow field to the current level
        currU = imresize(u, size(currIm1));
        currV = imresize(v, size(currIm1));
        
        % Warp im2 to the current estimate of the flow
        [x, y] = meshgrid(1:size(currIm2, 2), 1:size(currIm2, 1));
        warpedIm2 = interp2(currIm2, x + currU, y + currV);
        
        % Compute Lucas-Kanade optical flow at the current pyramid level
        [du, dv] = LucasKanade(currIm1, warpedIm2, windowSize);
        
        % Update the flow field
        u = currU + du;
        v = currV + dv;
    end
end

function pyramid = createPyramid(im, numLevels)
    pyramid = cell(1, numLevels);
    pyramid{1} = im;
    for i = 2:numLevels
        pyramid{i} = impyramid(pyramid{i-1}, 'reduce');
    end
end


