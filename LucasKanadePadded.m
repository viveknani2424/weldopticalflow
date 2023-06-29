function [du, dv] = LucasKanadePadded(im1, im2, windowSize)

    % Pad images to reducing the edge smoothening in lucas kanade algorithm
    im1 = padarray(im1, [windowSize windowSize], 'replicate');
    im2 = padarray(im2, [windowSize windowSize], 'replicate');

    % Compute image gradients
    [Ix, Iy] = gradient(double(im1));
    It = double(im2) - double(im1);

    % Initialize du and dv
    du = zeros(size(im1));
    dv = zeros(size(im1));

    % For each pixel in the image
    for y = windowSize+1:size(im1, 1)-windowSize
        for x = windowSize+1:size(im1, 2)-windowSize

            % Get window around pixel
            xRange = (x-windowSize):(x+windowSize);
            yRange = (y-windowSize):(y+windowSize);
            IxWindow = Ix(yRange, xRange);
            IyWindow = Iy(yRange, xRange);
            ItWindow = It(yRange, xRange);

            % Compute elements of the system of linear equations
            IxIx = sum(sum(IxWindow.^2));
            IyIy = sum(sum(IyWindow.^2));
            IxIy = sum(sum(IxWindow.*IyWindow));
            IxIt = sum(sum(IxWindow.*ItWindow));
            IyIt = sum(sum(IyWindow.*ItWindow));

            % Solve system of linear equations
            A = [IxIx IxIy; IxIy IyIy];
            b = [-IxIt; -IyIt];
            flow = pinv(A) * b;

            % Update flow
            du(y, x) = flow(1);
            dv(y, x) = flow(2);
        end
    end

    % Remove padding from du and dv
    du = du(windowSize+1:end-windowSize, windowSize+1:end-windowSize);
    dv = dv(windowSize+1:end-windowSize, windowSize+1:end-windowSize);
end
