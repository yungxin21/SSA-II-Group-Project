continuumimg = fitsread("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");
% note: image processing toolbox is required
img = fitsread("hmi.M_720s.20190413_210000_TAI.3.magnetogram.fits", "image");
maginfo = fitsinfo("hmi.M_720s.20190413_210000_TAI.3.magnetogram.fits");
% continuumimg = fitsread("hmi.Ic_720s.20130419_160000_TAI.1.continuum.fits", "image");
continfo = fitsinfo("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");

% raw image
figure(1); imagesc(img); colormap("bone"); colorbar; % will probably be white because of problems with colour scaling

% 1. solar disc
% assuming 
info = fitsinfo("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");
keywords = info.Image.Keywords; % retrieve info from HDU file.
valIdx1 = strcmp(keywords(:,1), 'CDELT1');
pixelScale = keywords{valIdx1, 2}
center_x = 2048;
center_y = 2048;
% find the first and last index with positive values larger than 1000
radius = 1900; % empirically chosen value which trims limb within an acceptable amount

% 2. making background = NaN so there is differentiation
% between background and solar disc
[cols, rows] = meshgrid(1:size(img,2), 1:size(img,1));
maskimg = sqrt((cols - center_x).^2 + (rows - center_y).^2) <= radius;
% median filtering
medianimg = medfilt2(img, [3 3]);
% apply the mask to the median filtered image to isolate the solar disc
medianimg(~maskimg) = NaN; 
figure(2); imagesc(medianimg); colormap("bone"); colorbar;

% 3. proceeding with median filtering (clearest filtering):
fluxthreshold = 120; % to be adjusted depending on avg magnetic field strength (TBD)
% maybe set a flux threshold depending on the images time of
% observation within the Hale's Law cycle. 
minimumarea = 90; % num of pixels to say its a region
% apply the threshold to identify field line regions
positivemask = medianimg > fluxthreshold; % white regions of strong lines out of disc
negativemask = medianimg < -fluxthreshold; % black regions of strong lines into disc 
% remove areas of fewer than minimumarea pixels of "flux"
positivemask = bwareaopen(positivemask, minimumarea);
negativemask = bwareaopen(negativemask, minimumarea);
% positive mask & negative mask put together (all flux regions)
maskedimg = positivemask | negativemask;

% 4. fitting an ellipse
% find connected regions (bipolar regions)
cc = bwconncomp(maskedimg);
% gives array of potential centroids of bipolar regions to fit elipses at
S = regionprops(cc, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', ...
    'Area');
% comparing with continuum img taken same time
% note: as hmi images only goes back so far, there is a limited range which we can
% compare with magnetograms.
figure(3); imagesc(continuumimg); colormap("bone"); colorbar; hold on;
theta = linspace(0, 2*pi, 200);
for k = 1:numel(S)
    cx = S(k).Centroid(1); % centroid x coord
    cy = S(k).Centroid(2); % centroid y coord 
    a = S(k).MajorAxisLength/2;
    b = S(k).MinorAxisLength/2;
    angle = deg2rad(-S(k).Orientation);
    % making a parametric elipse:
    x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
    y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
    plot(x, y, 'r', 'LineWidth', 1);
end

hold off;
figure(4); imagesc(medianimg); colormap("bone"); colorbar; hold on;
for k = 1:numel(S)
    cx = S(k).Centroid(1); % centroid x coord
    cy = S(k).Centroid(2); % centroid y coord 
    a = S(k).MajorAxisLength/2; % major axis
    b = S(k).MinorAxisLength/2; % minor axis
    angle = deg2rad(-S(k).Orientation); % orientation of ellipse
    % making a parametric ellipse:
    x = cx + a*cos(theta)*cos(angle) - b*sin(theta)*sin(angle);
    y = cy + a*cos(theta)*sin(angle) + b*sin(theta)*cos(angle);
    plot(x, y, 'r', 'LineWidth', 1);
end
hold off;

% 5. sunspot group/individual spot identification using ellipses & also
% cleanup bad spots

cdelt1 = pixelScale;
minsize_arcsec = 2; % ellipse must be at least this many arcseconds big to be counted
groupseperation_arcsec = 80; % regions closer than this are linked into one group (emperical)
groupseparation_px = groupseperation_arcsec / cdelt1; % max ellipse separation in pixels
mingrouparea = 150; % minimum group flux area (arcsec^2) (emperical and might change)

% the higher this is, more noise and limb darkening effects might be kept
spotfrac = 0.95; % obviously limb darkening isnt an issue when the other code is applied

pixarea = cdelt1 * cdelt2; % arcsec^2 per pixel

% span of an ellipse (semi-axes a,b, angle an) along x and y, in arcsec
span = @(a,b,an) [2*sqrt((a*cos(an))^2 + (b*sin(an))^2)*abs(cdelt1), ...
    2*sqrt((a*sin(an))^2 + (b*cos(an))^2)*abs(cdelt2)];
Iqs = median(continuumimg(isfinite(continuumimg) & continuumimg > 0)); % median brightness of whole image
spotmask = isfinite(continuumimg) & continuumimg > 0 & continuumimg < spotfrac*Iqs;
keep = find(cellfun(@(p) any(spotmask(p)), cc.PixelIdxList)); % keep flux regions that contain a sunspot (removes anomalies)

% drop ellipses with size less than minimum from the max of major/minor axis or orientation
bigenough = arrayfun(@(s) max(span(s.MajorAxisLength/2, s.MinorAxisLength/2, ...
    deg2rad(-s.Orientation))) >= minsize_arcsec, S); % logical array which checks if every ellipse is large enough
keep = keep(bigenough(keep)); % new array of ellipses without the ones that were too small unlike cc array
keptmask = false(size(maskedimg));

% pixels that belong to an ellipse
keptmask(vertcat(cc.PixelIdxList{keep})) = true;

% link the individual ellipses whose edges are within groupseparation_px of each other
GL = labelmatrix(bwconncomp(bwdist(keptmask) <= groupseparation_px/2));
gid = zeros(1, cc.NumObjects);
for k = keep
    gid(k) = GL(cc.PixelIdxList{k}(1)); 
end

% one enclosing ellipse per group
ids = unique(gid(keep)); % avoid duplicates
GE = struct('c',{}, 'a',{}, 'b',{}, 'ang',{}, 'members',{}, 'area',{}); % struct containing information about each group
% assigning ellipse parameters for each group 
for i = 1:numel(ids)
    m = find(gid == ids(i)); % member ellipses of this group
    grouparea = sum([S(m).Area]) * pixarea; % group flux area in arcsec^2
    if grouparea < mingrouparea % if group is too small get rid of it
        continue;
    end 
    c = mean(vertcat(S(m).Centroid), 1); % average of individual ellipse centroids
    P = []; % array of boundary points to give orientation and size of the new group ellipse
    for k = m
        P = [P; ellipsepoints(S(k).Centroid, S(k).MajorAxisLength/2, ...
            S(k).MinorAxisLength/2, deg2rad(-S(k).Orientation), theta).'];
    end
    Q = P - c; % distance of points from centroid of the ellipse
    [V, D] = eig(Q.'*Q/size(Q,1)); % points position along the ellipse's semimajor axis
    [~, o] = sort(diag(D), 'descend');
    u = Q*V(:,o(1));  w = Q*V(:,o(2));
    a0 = max(abs(u));  b0 = max(max(abs(w)), 1); % best semi major/minor axis for the ellipse (max means points on edge)
    s = max(sqrt((u/a0).^2 + (w/b0).^2)); % scale so every member fits inside
    GE(end+1) = struct('c',c, 'a',1.05*s*a0, 'b',1.05*s*b0, ...
        'ang',atan2(V(2,o(1)), V(1,o(1))), 'members',m, 'area',grouparea); 
end
% red ellipses (individual) + blue ellipses surrounding it (group)
figure(5); clf; imagesc(continuumimg); colormap("bone"); colorbar;
axis image; hold on
for i = 1:numel(GE)
    for k = GE(i).members
        p = ellipsepoints(S(k).Centroid, S(k).MajorAxisLength/2, S(k).MinorAxisLength/2, ...
            deg2rad(-S(k).Orientation), theta);
        plot(p(1,:), p(2,:), 'r', 'LineWidth', 1);
    end
    p = ellipsepoints(GE(i).c, GE(i).a, GE(i).b, GE(i).ang, theta);
    plot(p(1,:), p(2,:), 'b', 'LineWidth', 2);
end
hold off
figure(4); imagesc(medianimg); colormap("bone"); colorbar;


% group number:
groupnumber = numel(GE);
fprintf("Sunspot groups: %d\n", groupnumber);

% this handles the parameters for every group ellipse
function p = ellipsepoints(c, a, b, an, theta)
p = [c(1) + a*cos(theta)*cos(an) - b*sin(theta)*sin(an); ...
    c(2) + a*cos(theta)*sin(an) + b*sin(theta)*cos(an)];
end
