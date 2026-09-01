uncompressedimg = gunzip("")
img = fitsread("k4v72090930t165755.fits");
imginfo = fitsinfo("k4v72090930t165755.fits");

img = imgaussfilt(img, 2); % removing background  noise
fitsdisp(img, "image");

fluxthreshold = 50;

% 1 image per 3 months? so over a 44 year period that is 136 images to
% verify hales law
% ^ increment may need to be reduced by every 1 month but we'll see
% 1) find bipolar magnetic field regions
% note: white = lines out of sun, black = lines into the sun, grey =
% no/little lines 
