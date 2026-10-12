close all
clear all

% find how close sunspots have to be to count as a group using Zurich
% classification.
% A, B, and C types are assumed to be between 1-10 degrees of solar
% longitude.
% D: <=10 Ls
% E: 10-15 Ls
% F: >15 Ls
% G: 1 to >15 Ls
% H: 1 sunspot that is >2.5 Ls (1 to >15Ls)
% J: 1 sunspot that is <2.5 Ls

% Method: use an ellipse to group sunspots. 
% ellipse function: x^2/a^2 + y^2/b^2 = 1
% a is semi-major axis (horizontal)
% b is semi-minor axis (vertical)
% find the right eccentricity.

% ellipse in polar 
% r^2 = ((r^2)*(cos(theta))^2)/a^2 + ((r^2)*(sin(theta))^2)/b^2;

[img,radius_sun] = limbDarkness("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits",0,0);
% need to clip edges of the Sun to exclude background space (producing
% intensity values of ~0)

info = fitsinfo("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");
keywords = info.Image.Keywords; % retrieve info from HDU file.
valIdx1 = strcmp(keywords(:,1), 'CDELT1');
valIdx2 = strcmp(keywords(:,1), 'CRPIX1');
valIdx3 = strcmp(keywords(:,1), 'CRPIX2');
pixelScale = keywords{valIdx1, 2}; % arcsec/pixel redundant?
% x and y centers of the Sun. 
center_x = round(keywords{valIdx2, 2});
center_y = round(keywords{valIdx3, 2});

% %% 
% % checking for sigma values
% bins = 0:0.00001:1.3;
% [count,bins] = hist(img(img > 0), bins);
% figure(6); semilogy(bins,(count), Color="g");
%%

% create a mask that "clips" edges of the Sun.
[X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
% Limitation: Limb edges are discarded, so if a sunspot is present there,
% it is omitted.
sunMask = sqrt((X - center_x).^2 + (Y - center_y).^2) < (0.99*radius_sun);
% Restrict segmentation to the solar disk and remove background pixels.
img(~sunMask) = NaN;
clippedImg = sigmaClip(img,0,1); % full image.
[penumbraImg, penumbraMask] = sigmaClip(img,1,1);
[umbraImg, umbraMask] = sigmaClip(img,2,1);

% bwareaopen is currently not working as expected
penumbraMask = bwareaopen(penumbraMask, 75);
umbraMask = bwareaopen(umbraMask, 225);
mergedMask = (penumbraMask | umbraMask);

%%

% figure(1); imagesc(penumbraImg+umbraImg); colormap; colorbar; set(gca, "YDir", "normal");
figure(2); imagesc(clippedImg); colormap; colorbar; set(gca, "YDir", "normal"); hold on;

% hold off; % temp line
% locate places to place ellipsoids, retrieved from magnetograms branch
cc = bwconncomp((penumbraMask == 0) & sunMask);
S = regionprops(cc, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Area');

% draw individual sunspot ellipsoids
% current problem: centroids within other ellipsoids still produce an
% ellipsoid. Since bwconncomp is just BFS, isolated penumbra will be
% affected.

% solution: use rangesearch to find neighboring centroids within
% major/minor axis length, and essentially shifting the centroid to produce a new ellipsoid.


theta = linspace(0,2*pi,1000);
for k = 1:numel(S)
    cx = S(k).Centroid(1); % centroid x coord
    cy = S(k).Centroid(2); % centroid y coord
    a = S(k).MajorAxisLength/2;
    b = S(k).MinorAxisLength/2;
    angle = deg2rad(-S(k).Orientation);
    % making a parametric ellipse:
    x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
    y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
    plot(x, y, 'r-', 'LineWidth', 1);
end

% create group ellipsoids
% at the moment minGroupArea is constant and very static. Not good!
minGroupArea = 1000;
% idea is to find neighboring centroids within a certain distance, then
% merging them into one group and redrawing new ellipsoids on top of
% existing ones
g = ;

hold off;

% Wolf Number
% still need to figure out group number g
f=cc.NumObjects;
sunspotNum = 10*g + f;
disp(sunspotNum);