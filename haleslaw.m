% 1 image per 3 months? so over a 44 year period that is 136 images to
% verify hales law
% ^ increment may need to be reduced by every 1 month but we'll see

% note: white = lines out of sun, black = lines into the sun, grey =
% no/little lines 

% process
% 1) smooth the image (gauss filtering to remove noise while preserving the
% bipolar active regions on the magnetogram).
% 

% NOTE: raw archive files from the 512-channel spectrograph from kitt peak space telescope give
% fits files in 2048x2048x5 (3D array). the image is the same in all 5
% arrays
img = fitsread("512c_eo000_C1_19760116_1620.fits");
% raw image
figure(1); imagesc(img(:,:,1)); colormap("bone"); colorbar;

imgframe = img(:,:,1);
img = medfilt2(imgframe,[3 3]);

% note: need image processing toolbox to run imgaussfilt
medianimg = medfilt2(img, [3 3]);
gaussimg = imgaussfilt(img, 2);

fluxthreshold = 50; % to be adjusted depending on avg magnetic field strength (TBD)
% maybe set a flux threshold depending on the images time of
% observation within the Hale's Law cycle. 
minimumarea = 50; % num of pixels to keep in the patch of bipolar 

% apply the threshold to identify field line regions
%positivemask = smoothedimg > fluxthreshold; % white regions of strong positive polarity
%negativemask = smoothedimg < -fluxthreshold; % black regions of strong negative polarity

% gauss filtering
figure(2); imagesc(gaussimg); colormap("bone"); colorbar;

% median filtering
figure(3); imagesc(medianimg); colormap("bone"); colorbar;
