% Read in fits image of spotless sun
%spotlessimg=fitsread("SSAII\hmi.Ic_noLimbDark_720s.20191203_010000_TAI.3.continuum.fits","image");

% Preallocating bin size
bins = 0:0.01:1.3;

% Histogram of log10 intensity
%[spotlesscount,bins] = hist(spotlessimg(spotlessimg > 0), bins);
%semilogy(bins,(spotlesscount), Color="g"); hold on;

% img=fitsread("SSAII\hmi.Ic_noLimbDark_720s.20240504_040000_TAI.3.continuum.fits","image");
%img=fitsread("corrected.fits","image");
[spotcount,bins] = hist(img(img > 0), bins);
%semilogy(bins,(spotcount), Color="b"); hold on;

% Generating mode and interval for seperating sun from spots
[intensityCentre,binsCentre] = max((spotcount));
halfMax = intensityCentre/2;

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
maskedImg = data;
maskedImg(~mask & ~nanMask) = -1;
end

% Testing function
clippedImg = sigmaClip(corrected_img, 3, 1);
% Using imagesc to compare to DS9 fits image
%imagesc(clippedImg), colorbar, clim([-2 0]);
% Image was flipped vertically so just flipped rows of the matrix to make
% the same as the DS9 fits image
flippedImg = flipud(clippedImg);
imagesc(flippedImg), colorbar, clim([-2 0]);
% Imagesc produces the same as DS9 showing the same sunspots