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

info = fitsinfo("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");
img = limbDarkness("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits",0,0);
keywords = info.Image.Keywords; % retrieve info from HDU file.
valIdx1 = strcmp(keywords(:,1), 'CDELT1');
valIdx2 = strcmp(keywords(:,1), 'CRPIX1');
valIdx3 = strcmp(keywords(:,1), 'CRPIX2');
pixelScale = keywords{valIdx1, 2}; % Extract arcsec/pixel conversion for later
% x and y centers of the Sun.
center_x = round(keywords{valIdx2, 2});
center_y = round(keywords{valIdx3, 2});

% checking for sigma values
bins = 0:0.01:1.3;
[count,bins] = hist(img(img > 0), bins);
figure(2); semilogy(bins,(count), Color="g");

%%
clippedImg = sigmaClip(img,5,10);
disp("done clippedImg");
[penumbraImg, penumbraMask] = sigmaClip(img,6,1);
disp("done penumbra");
[umbraImg, umbraMask] = sigmaClip(img,17,1);
disp("done umbra");
mergedMask = (penumbraMask+umbraMask);

figure(1); imagesc(clippedImg); colormap; colorbar; set(gca, "YDir", "normal"); 
figure(3); imagesc(mergedMask); axis xy, colorbar;

% locate places to place ellipsoids, retrieved from magnetograms branch
cc = bwconncomp(penumbraMask);
S = regionprops(cc, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Area');

theta = linspace(0,2*pi,1000);
for k = 1:numel(S)
    cx = 4096-S(k).Centroid(1); % centroid x coord
    cy = 4096-S(k).Centroid(2); % centroid y coord
    a = S(k).MajorAxisLength/2;
    b = S(k).MinorAxisLength/2;
    angle = deg2rad(-S(k).Orientation);
    % making a parametric ellipse:
    x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
    y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
    plot(x, y, 'r-', 'LineWidth', 1);
end

% Wolf Number
% still need to figure out f
f=0;
sunspotNum = 10*numel(S) + f;