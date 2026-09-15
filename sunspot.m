% Read in fits image of spotless sun
%spotlessimg=fitsread("SSAII\hmi.Ic_noLimbDark_720s.20191203_010000_TAI.3.continuum.fits","image");

% Preallocating bin size
%bins = 0:0.01:1.3;

% Histogram of log10 intensity
%[spotlesscount,bins] = hist(spotlessimg(spotlessimg > 0), bins);
%semilogy(bins,(spotlesscount), Color="g"); hold on;

% Read in fits image of no limb dark sun for testing
img=fitsread("SSAII\hmi.Ic_noLimbDark_720s.20240504_040000_TAI.3.continuum.fits","image");
%img=fitsread("corrected.fits","image");
%[spotcount,bins] = hist(img(img > 0), bins);
%semilogy(bins,(spotcount), Color="b"); hold on;

% Generating mode and interval for seperating sun from spots
%[intensityCentre,binsCentre] = max((spotcount));
%halfMax = intensityCentre/2;

% Tried using fwhm for the automatic detection of sun values but pivotted
% to sigma clipping instead
%lIndex = (bins<=1);
%rIndex = (bins>=1);
%lxIntercept = interp1(spotcount,halfMax);
%yline(halfMax,'--r');

% Finding sd of data
%sigma = std(img(img>0));

% Setting n sigma as value thereshold
% Commented out while testing different method.
%n = 3;
%lBound = bins(binsCentre) - (n*sigma);
%rBound = bins(binsCentre) + (n*sigma);
%yBounds = (bins(binsCentre) >= lBound) & (bins(binsCentre) <= rBound);
%yclipped = spotcount(yBounds);
%figure(2);plot(yclipped); 

% Using rmoutliers function to set bounds of data for sigma clipping around
% mode
%yClipped = rmoutliers(spotcount,"median","ThresholdFactor=3");
% A lot of errors so possibly not possible with that inbuilt function

% Attempting to create own function for clipping with help from AI
function [maskedImg, mask] = sigmaClip(data, nSigma, maxI)
% Make all negative values to be NaN for ease of finding 0 in data
maskedImg = logical(data);
maskedImg(:) = 0;
data(data<0) = NaN;
% Dedicated map of where NaN's are
nanMask = isnan(data);
% Start mask (true = valid, false = clipped)
mask = true(size(data));
for i=1:maxI
    % Make current unclipped values as 1D vector
    validData = data(mask & ~nanMask);
    if isempty(validData), break; end
    % Running statistics on valid data and ommiting clipped values
    centreVal = mode(validData);
    stdVal = std(validData);
    if stdVal == 0, break; end
    % Defining boundaries
    lBound = centreVal - nSigma * stdVal;
    % Use bounds to generate new mask 
    newMask = mask & (data >= lBound);
    if isequal(newMask,mask), break; end
    mask = newMask;
end
% Replacing clipped values to keep same matrix dimensions
maskedImg(~mask & ~nanMask) = 1;
end

% Testing function
%clippedImg = sigmaClip(img, 5, 10);
% Using imagesc to compare to DS9 fits image
%imagesc(clippedImg), colorbar, clim([-2 0]);
% Image was flipped vertically so just flipped rows of the matrix to make
% the same as the DS9 fits image
%flippedImg = flipud(clippedImg);
% Altered the imagesc function instead
%imagesc(clippedImg);axis xy, colorbar, clim([-2 0]);
% Imagesc produces the same as DS9 showing the same sunspots

% For sunspot group recognition need to show a difference between penumbra
% and umbra, will test merging two images one with higher sigma and one
% with lower. Tried to add values directly but eventually came to the
% conclusion that adding the masks would work

[penumbraImg,penumbraMask] = sigmaClip(img,6,1);
[umbraImg,umbraMask] = sigmaClip(img,12,1);
mergedMask = (penumbraMask+umbraMask);
imagesc(mergedMask);axis xy, colorbar;

% Sunspot number equation is R=k(10g+f) where f is total number of
% individual sunspots, g is the number of clusters and k is the scaling
% factor

% researched flood fill algorithms and trying to code a BFS queue type for
% identifying how many umbra and penumbra are present in created mask

function [penumbraCount, umbraCount] = floodFillCount(matrix)
    % setting counts to 0 for penumbra and umbra to -1 to account for the
    % background value
    penumbraCount = 0;
    umbraCount = -1;
    % getting size of input for boundary for function
    [numRows, numCols] = size(matrix);
    
    % creating copy of matrix to change values to not search values more than
    % once
    copyMatrix = matrix;
    
    % setting directions for flooding
    directions = [1 0; -1 0; 0 1; 0 -1];
    
    % for loop to search each pixel in the matrix
    for r = 1:numRows
        for c = 1:numCols
            if copyMatrix(r,c) == 0
                umbraCount = umbraCount + 1;
    
                % starting queue for penumbra
                umbraQueue = [r, c];
    
                % setup loop to run while there are values in queue
                while ~isempty(umbraQueue)
    
                    % setting up current coordinates
                    currentCell = umbraQueue(1,:);
                    umbraQueue(1,:) = [];
                    currentRow = currentCell(1);
                    currentCol = currentCell(2);
    
                    % changes current coordinate in loop to -1 to avoid value in loop again
                    if copyMatrix(currentRow,currentCol) == 0
                        copyMatrix(currentRow,currentCol) = 1;
    
                        % for loop for searching neighbour values
                        for d = 1:size(directions,1)
                            nextRow = currentRow + directions(d,1);
                            nextCol = currentCol + directions(d,2);
    
                            % checking if neighbour values are within matrix bounds
                            if nextRow >= 1 && nextRow <= numRows && nextCol >= 1 && nextCol <= numCols
                                % if statement for when equal to 1 add to the queue to
                                % search further from there
                                if copyMatrix(nextRow,nextCol) == 0
                                    umbraQueue = [umbraQueue; nextRow, nextCol];
                                end
                            end
                        end
                    end
                end
            end
            if copyMatrix(r,c) == 1
                penumbraCount = penumbraCount + 1;
    
                % starting queue for penumbra
                penumbraQueue = [r, c];
    
                % setup loop to run while there are values in queue
                while ~isempty(penumbraQueue)
    
                    % setting up current coordinates
                    currentCell = penumbraQueue(1,:);
                    penumbraQueue(1,:) = [];
                    currentRow = currentCell(1);
                    currentCol = currentCell(2);
    
                    % changes current coordinate in loop to -1 to avoid value in loop again
                    if copyMatrix(currentRow,currentCol) == 1
                        copyMatrix(currentRow,currentCol) = -1;
    
                        % for loop for searching neighbour values
                        for d = 1:size(directions,1)
                            nextRow = currentRow + directions(d,1);
                            nextCol = currentCol + directions(d,2);
    
                            % checking if neighbour values are within matrix bounds
                            if nextRow >= 1 && nextRow <= numRows && nextCol >= 1 && nextCol <= numCols
                                % if statement for when equal to 1 add to the queue to
                                % search further from there
                                if copyMatrix(nextRow,nextCol) == 1
                                    penumbraQueue = [penumbraQueue; nextRow, nextCol];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
    
%[penumbraCount, umbraCount] = floodFillCount(mergedMask);