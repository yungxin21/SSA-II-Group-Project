img = fitsread("512c_eo000_C1_19760116_1620.fits");
fitsdisp("512c_eo000_C1_19760116_1620.fits");
% note: need image processing toolbox to run imgaussfilt
smoothedimg = imgaussfilt(img, 2);

fluxthreshold = 50; % to be adjusted depending on magnetic field strength (TBD)
% apply the threshold to identify field line regions
positivemask = smoothedimg > fluxthreshold;
negativemask = smoothedimg < -fluxthreshold;



% 1 image per 3 months? so over a 44 year period that is 136 images to
% verify hales law
% ^ increment may need to be reduced by every 1 month but we'll see
% 1) find bipolar magnetic field regions
% note: white = lines out of sun, black = lines into the sun, grey =
% no/little lines 
