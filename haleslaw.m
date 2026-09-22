% note: image processing toolbox is required
% 1 image per 3 months? so over a 44 year period that is 136 images to
% verify hales law
% ^ increment may need to be reduced by every 1 month but we'll see

% note: white = lines out of sun, black = lines into the sun, grey =
% no/little lines 

% process
% 1) smooth the image (gauss filtering to remove noise while preserving the
% bipolar active regions on the magnetogram).
% 2) fit an elipsoid to the image to find bipolar flux regions

% NOTE: raw archive files from the 512-channel spectrograph from kitt peak space telescope give
% fits files in 2048x2048x5 (3D array). use 1st channel (:,:,1) by default
img = fitsread("svsm_8542l_S2_20130419_1605.fits");
continuumimg = fitsread("hmi.Ic_720s.20130419_160000_TAI.1.continuum.fits", "image");
continuumsizedimg = imresize(continuumimg,0.5,"bilinear");
img = img(:,:,1);
% raw image
figure(1); imagesc(img); colormap("bone"); colorbar;

% 1. identify the solar disc (code courtesy of the limbdarkness branch)
% assuming 
center_x = 1024;
center_y = 1024;
% find the first and last index with positive values larger than 1000
rowvals = img(center_y, :); % row at centre_y
colvals = img(: , center_x); % col at center_x
xmin = find(rowvals/1000 >= 0, 1, 'first'); 
xmax = find(rowvals/1000 >= 0, 1, 'last');
ymin = find(colvals/1000 >= 0, 1, 'first');
ymax = find(colvals/1000 >= 0, 1, 'last');

diameter_x = xmax - xmin + 1;
radius_x = diameter_x / 2;
diameter_y = ymax - ymin + 1;
radius_y = diameter_y / 2;
radius = floor(mean([radius_x,radius_y]));
radius = 935;

% 2. making background = NaN so there is differentiation
% between background and solar disc
[cols, rows] = meshgrid(1:size(img,2), 1:size(img,1));
maskimg = sqrt((cols - center_x).^2 + (rows - center_y).^2) <= radius;
% median filtering
medianimg = medfilt2(img, [3 3]);
% apply the mask to the median filtered image to isolate the solar disc
medianimg(~maskimg) = NaN; 
figure(2); imagesc(medianimg); colormap("bone"); colorbar; axis xy;

% 3. proceeding with median filtering (clearest filtering):
fluxthreshold = 120; % to be adjusted depending on avg magnetic field strength (TBD)
% maybe set a flux threshold depending on the images time of
% observation within the Hale's Law cycle. 
minimumarea = 75; % num of pixels to keep in the patch of flux
% apply the threshold to identify field line regions
positivemask = medianimg > fluxthreshold; % white regions of strong lines out of disc
negativemask = medianimg < -fluxthreshold; % black regions of strong lines into disc 
% remove areas of fewer than minimumarea pixels of "flux"
positivemask = bwareaopen(positivemask, minimumarea);
negativemask = bwareaopen(negativemask, minimumarea);
% positive mask & negative mask put together
maskedimg = positivemask | negativemask;

% 4. fitting an elipse
% find connected regions (bipolar regions)
cc = bwconncomp(maskedimg);
% gives array of potential centroids of bipolar regions to fit elipses at
S = regionprops(cc, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', ...
    'Area');
% comparing with continuum img taken same time
% note: as hmi images only goes back so far, there is a limited range which we can
% compare with magnetograms.
figure(3); imagesc(continuumsizedimg); colormap("bone"); colorbar; hold on; axis xy;
theta = linspace(0, 2*pi, 200);
for k = 1:numel(S)
    cx = 2049 - S(k).Centroid(1); % centroid x coord
    cy = 2049 - S(k).Centroid(2); % centroid y coord 
    a = S(k).MajorAxisLength/2;
    b = S(k).MinorAxisLength/2;
    angle = deg2rad(-S(k).Orientation);
    % making a parametric elipse:
    x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
    y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
    plot(x, y, 'r', 'LineWidth', 1);
end

% 5. sunspot group/individual spot classification using elipses
hold off;
