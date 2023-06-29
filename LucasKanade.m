function [du, dv] = LucasKanade(im1, im2, windowSize)

    % Compute image gradients
    [Ix, Iy] = gradient(double(im1));
    It = double(im2) - double(im1);

    % Initialize du and dv
    du = zeros(size(im1));
    dv = zeros(size(im1));

    % For each pixel in the image
    for y = 1:size(im1, 1)
        for x = 1:size(im1, 2)

            % Get window around pixel
            xRange = max(1, x-windowSize):min(size(im1, 2), x+windowSize);
            yRange = max(1, y-windowSize):min(size(im1, 1), y+windowSize);
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
end
