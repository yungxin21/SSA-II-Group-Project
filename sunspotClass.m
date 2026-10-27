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

function [individuals, GE] = sunspotClass(img, radius_sun, center_x, center_y, pixelScale)
% [img,radius_sun] = limbDarkness("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits",0,0);
% % need to clip edges of the Sun to exclude background space (producing
% % intensity values of ~0)

% moved to main.m
% info = fitsinfo("hmi.Ic_720s.20240504_040000_TAI.3.continuum.fits");
% keywords = info.Image.Keywords; % retrieve info from HDU file.
% valIdx1 = strcmp(keywords(:,1), 'CDELT1');
% valIdx2 = strcmp(keywords(:,1), 'CDELT2');
% % valIdx3 = strcmp(keywords(:,1), 'CRPIX1');
% % valIdx4 = strcmp(keywords(:,1), 'CRPIX2');
% pixelScale = keywords{valIdx1, 2}; 
% pixelScale2 = keywords{valIdx2, 2};
% % % x and y centers of the Sun. 
% % center_x = round(keywords{valIdx3, 2});
% % center_y = round(keywords{valIdx4, 2});

% checking for sigma values (for debugging)
% bins = 0:0.00001:1.3;
% [count,bins] = hist(img(img > 0), bins);
% figure(6); semilogy(bins,(count), Color="g");

% create a mask that "clips" edges of the Sun.
% [X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
% Limitation: Limb edges are discarded, so if a sunspot is present there,
% it is omitted.
% sunMask = sqrt((X - center_x).^2 + (Y - center_y).^2) < (0.99*radius_sun);

% transform version
imgsz = size(img,1);
% sunMask = zeros(imgsz);
% Sporer's law states 95% of sunspots lie within 45 deg N and S
% heliographic latitudes.
% rows = floor(1/2*radius_sun+(imgsz/2)-radius_sun):1:ceil(3/2*radius_sun+(imgsz/2)-radius_sun);
% sunMask(rows,175:3921) = 1; % HARD CODED NEED TO CHANGE

% rounded rectangle (to ignore intensity falloff at corners)
% Center and dimensions of the rounded rectangle
w = 3700;  % total width (HARDCODED)
h = radius_sun;  % total height
r = 500;   % corner radius (HARDCODED)
% debugging
% disp(w+" "+h);

[X, Y] = meshgrid(1:imgsz, 1:imgsz);

% Distance from center
dx = abs(X - center_x) - (w/2 - r);
dy = abs(Y - center_y) - (h/2 - r);

% A pixel is inside the rounded rectangle if its:
% 1. inside the central cross (flat edges)
% 2. inside the 4 corner circles
insideFlatX = (abs(X - center_x) <= w/2) & (abs(Y - center_y) <= (h/2 - r));
insideFlatY = (abs(X - center_x) <= (w/2 - r)) & (abs(Y - center_y) <= h/2);
insideCorners = (max(dx, 0).^2 + max(dy, 0).^2) <= r^2;

% Combine to form the binary mask
sunMask = insideFlatX | insideFlatY | insideCorners;
% debugging
% figure(4); imagesc(sunMask);

% Restrict segmentation to the solar disk and remove background pixels.
clippedImg = img;
clippedImg(~sunMask) = NaN;
% clippedImg = sigmaClip(img,0,1); % full image.
[~, penumbraMask] = sigmaClip(clippedImg,1,1);
[~, umbraMask] = sigmaClip(clippedImg,2,1);

% need to change this hardcoding fir minimum area.
% 2019-04-13: 50, 100.
penumbraMask = ~bwareaopen(penumbraMask==0, 10); % HARDCODED
% umbraMask = ~bwareaopen(umbraMask==0, 5); % HARDCODED

% testing
% mergedMask = xor(penumbraMask, umbraMask);
% debugging
% figure(5); imagesc(mergedMask); colormap; colorbar; set(gca, 'YDir', 'normal');


% %%
% locate places to place ellipsoids, code retrieved from magnetograms branch
connectedRegions = bwconncomp((penumbraMask == 0) & sunMask);
individuals = regionprops(connectedRegions, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Area');
% big issue with sunspots with connected penumbra but 2 or more separated
% umbra being identified as one.

% find group number g
% create group ellipsoids
% idea is to find neighboring centroids within a certain distance, then
% merging them into one group and redrawing new ellipsoids on top of
% existing ones
% this becomes redundant if using transform
% ang_dist = (pixelScale * radius_sun)/90;
% distPx = (ang_dist*15)/pixelScale;
theta = linspace(0,2*pi,1000);
% Ls = (1.496*10^11) * sind(1/3600);
% disp("1: " + distPx + newline + "2: " + Ls);

% % Merge nearby regions into groups using centroid distance and area threshold.
% centroids = vertcat(individuals.Centroid);
% distances = pdist2(centroids, centroids);
% upperTriangle = triu(distances);
% % agglomerative hierarchical cluster tree (data structure)
% Z = linkage(upperTriangle, 'centroid'); % test methods here (2nd param)
% groups = cluster(Z, Cutoff=distPx, Criterion='distance');
% groupLabels = bwlabel(~penumbraMask);

% % groups = regionprops(groupLabels, 'Centroid', 'MajorAxisLength', ...
% %     'MinorAxisLength', 'Orientation', 'Area');
% uniqueGroups = unique(groups);
% newCentroids = zeros(length(uniqueGroups), size(centroids, 2));
% 
% for i = 1:length(uniqueGroups)
%     idx = (groups == uniqueGroups(i));
%     newCentroids(i, :) = mean(centroids(idx, :), 1);
% end

% for i = 1:numel(groups)
%     cx2 = groups(i).Centroid(1); % centroid x coord
%     cy2 = groups(i).Centroid(2); % centroid y coord
%     a2 = groups(i).MajorAxisLength/2;
%     b2 = groups(i).MinorAxisLength/2;
%     angle2 = deg2rad(-groups(i).Orientation);
%     % making a parametric ellipse:
%     x2 = cx2 + a2*cos(theta)*cos(angle2) - b2*sin(theta)*sin(angle2);
%     y2 = cy2 + a2*cos(theta)*sin(angle2) + b2*sin(theta)*cos(angle2);
%     plot(x2, y2, 'b-', 'LineWidth', 1);
% end




% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% how to prevent hardcoding here?
% answer: you dont
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
minsize_arcsec = 2;
groupseperation_arcsec = 80; 
distPx = groupseperation_arcsec / pixelScale; 
mingrouparea = 1;
pixarea = pixelScale * pixelScale;
% this anonymous function handles the parameters for every group ellipse
ellipsepoints = @(c,a,b,an,theta) [c(1) + a*cos(theta)*cos(an) - b*sin(theta)*sin(an); ...
    c(2) + a*cos(theta)*sin(an) + b*sin(theta)*cos(an)];
% % for reference
% function p = ellipsepoints(c, a, b, an, theta)
%     p = [c(1) + a*cos(theta)*cos(an) - b*sin(theta)*sin(an); ...
%         c(2) + a*cos(theta)*sin(an) + b*sin(theta)*cos(an)];
% end

% span of an ellipse (semi-axes a,b, angle an) along x and y, in arcsec
span = @(a,b,an) [2*sqrt((a*cos(an))^2 + (b*sin(an))^2)*abs(pixelScale), ...
    2*sqrt((a*sin(an))^2 + (b*cos(an))^2)*abs(pixelScale)];
Iqs = median(clippedImg(isfinite(clippedImg) & clippedImg > 0)); % median brightness of whole image
spotmask = isfinite(clippedImg) & clippedImg > 0 & clippedImg < Iqs;
keep = find(cellfun(@(p) any(spotmask(p)), connectedRegions.PixelIdxList)); % keep flux regions that contain a sunspot (removes anomalies)

% drop ellipses with size less than minimum from the max of major/minor axis or orientation
bigenough = arrayfun(@(s) max(span(s.MajorAxisLength/2, s.MinorAxisLength/2, ...
    deg2rad(-s.Orientation))) >= minsize_arcsec, individuals); % logical array which checks if every ellipse is large enough
keep = keep(bigenough(keep)); % new array of ellipses without the ones that were too small unlike cc array
keptmask = false(size(penumbraMask));

% pixels that belong to an ellipse
keptmask(vertcat(connectedRegions.PixelIdxList{keep})) = true;

% link the individual ellipses whose edges are within groupseparation_px of each other
GL = labelmatrix(bwconncomp(bwdist(keptmask) <= distPx/2));
gid = zeros(1, connectedRegions.NumObjects);
for k = keep
    gid(k) = GL(connectedRegions.PixelIdxList{k}(1)); 
end

% one enclosing ellipse per group
ids = unique(gid(keep)); % avoid duplicates
GE = struct('c',{}, 'a',{}, 'b',{}, 'ang',{}, 'members',{}, 'area',{}); % struct containing information about each group
% assigning ellipse parameters for each group 
for i = 1:numel(ids)
    m = find(gid == ids(i)); % member ellipses of this group
    grouparea = sum([individuals(m).Area]) * pixarea; % group flux area in arcsec^2
    if grouparea < mingrouparea % if group is too small get rid of it
        continue;
    end 
    % centers = vertcat(individuals(m).Centroid);
    % weights = [individuals(m).Area].';
    % c = sum(centers .* weights, 1) / sum(weights);
    c = mean(vertcat(individuals(m).Centroid), 1); % average of individual ellipse centroids
    P = []; % array of boundary points to give orientation and size of the new group ellipse
    for k = m
        P = [P; ellipsepoints(individuals(k).Centroid, individuals(k).MajorAxisLength/2, ...
            individuals(k).MinorAxisLength/2, deg2rad(-individuals(k).Orientation), theta).'];
    end
    Q = P - c; % distance of points from centroid of the ellipse
    [V, D] = eig(Q.'*Q/size(Q,1)); % points position along the ellipse's semimajor axis
    [~, o] = sort(diag(D), 'descend');
    u = Q*V(:,o(1));  w = Q*V(:,o(2));
    a0 = max(abs(u));  b0 = max(max(abs(w)), 1); % best semi major/minor axis for the ellipse
    s = max(sqrt((u/a0).^2 + (w/b0).^2)); % scale so every member fits inside
    GE(end+1) = struct('c',c, 'a',1.05*s*a0, 'b',1.05*s*b0, ...
        'ang',atan2(V(2,o(1)), V(1,o(1))), 'members',m, 'area',grouparea); 
end

% figure(1); imagesc(penumbraImg+umbraImg); colormap; colorbar; set(gca, "YDir", "normal");
% figure(1); imagesc(penumbraMask); colormap; colorbar; set(gca, "YDir", "normal");
% figure(2); imagesc(umbraMask); colormap; colorbar; set(gca, "YDir", "normal");
% figure(3); imagesc(img); colormap; colorbar; set(gca, "YDir", "normal"); 
% hold on;

% for k = 1:numel(individuals)
%     cx = individuals(k).Centroid(1); % centroid x coord
%     cy = individuals(k).Centroid(2); % centroid y coord
%     a = individuals(k).MajorAxisLength/2;
%     b = individuals(k).MinorAxisLength/2;
%     angle = deg2rad(-individuals(k).Orientation);
%     % making a parametric ellipse:
%     x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
%     y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
%     plot(x, y, 'r-', 'LineWidth', 1);
% end
% for i = 1:numel(GE)
%     for k = GE(i).members
%         p = ellipsepoints(individuals(k).Centroid, individuals(k).MajorAxisLength/2, ...
%             individuals(k).MinorAxisLength/2, ...
%             deg2rad(-individuals(k).Orientation), theta);
%         plot(p(1,:), p(2,:), 'r', 'LineWidth', 1);
%     end
%     p = ellipsepoints(GE(i).c, GE(i).a, GE(i).b, GE(i).ang, theta);
%     plot(p(1,:), p(2,:), 'b', 'LineWidth', 2);
% end
% hold off;

end