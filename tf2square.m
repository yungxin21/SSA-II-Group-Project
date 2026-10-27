% function [data] = coordCorrection (filename)
% coordCorrection %This function takes the pixel coordinates of the sunspot images and corrects
% %for the spherical nature of the shape
% 
% data = fitsread(filename , "image"); 
% 
% 
% end
function [map_norm] = tf2square(data,solarRadius,u_center,v_center, phi_steps, theta_steps)

% info = fitsinfo("hmi.Ic_720s.20190413_210000_TAI.3.continuum.fits");
% keywords = info.Image.Keywords; % retrieve info from HDU file.
% valIdx1 = strcmp(keywords(:,1), 'CDELT1');
% valIdx2 = strcmp(keywords(:,1), 'CRPIX1');
% valIdx3 = strcmp(keywords(:,1), 'CRPIX2');
% pixelScale = keywords{valIdx1, 2}; % arcsec/pixel redundant?
% % x and y centers of the Sun. 
% u_center = round(keywords{valIdx2, 2});
% v_center = round(keywords{valIdx3, 2});

% This takes the array of data from the fits file and composes two seperate
% arrays called "row" and "col". These are effectively 2 arrays with the
% row numbers and column numbers of every pixel co-ordinate of the original
% fits file. 
[row, col]=(ndgrid(1:size(data,1),1:size(data,2)));

%This takes all the intensity values from the data and organises them all
%into a single column that will line the corresponding values up with its
%corresponding row and column number.

% x - 2048.5, 2048.5 - y is done to centre the data properly.

u_centered = col - u_center;
v_centered = row - v_center;

%This defines the x and y transformed co-ordinate arrays as an array of
%zeroes so that they can be iterated on in the transform for loop.
phi_trans = NaN(size(v_centered,1), size(v_centered,2));
theta_trans = NaN(size(u_centered,1), size(u_centered,2));
% z_trans = zeros(size(u_centered,1), size(u_centered,2));

z_heeq = zeros(size(v_centered,1), size(v_centered,2)); 
% Defines heliocentric co-ordinates for the Z axis.

% phi_steps = 4096;
% theta_steps = 4096;
% map = zeros(phi_steps, theta_steps); % comment out if not using nested loop
% map_count = zeros(phi_steps, theta_steps); % comment out if not using nested loop
% map_norm = zeros(phi_steps, theta_steps);
% Defines the map array used to imagesc the transformed image of the sun.
% phi_steps and theta_steps alters the resolution of the image.

% phi_array = [-pi/2: pi/phi_steps : pi/2];
% theta_array = [-pi/2 : pi/theta_steps : pi/2];
% phi_grad = pi; % found by calculating it manually outside of code, just y=mx+c
% theta_grad = pi;
 
% try and use masks here
z_sq = solarRadius^2-(v_centered.^2+u_centered.^2);
% Restrict the transform to pixels projected onto the visible solar disk.
validMask = z_sq >= 0 & isfinite(data);
z_heeq(validMask) = sqrt(z_sq(validMask));
% z(~validMask) = NaN;
% Convert valid projected coordinates to spherical longitude and latitude.
phi_trans(validMask) = atan(v_centered(validMask)./sqrt(z_heeq(validMask).^2 + u_centered(validMask).^2));
theta_trans(validMask) = atan2(u_centered(validMask), z_heeq(validMask));
phi_map = round(1 + (phi_trans(validMask) + pi/2)./pi .*(phi_steps - 1)); % check
theta_map = round(1 + (theta_trans(validMask) + pi/2)./pi .*(theta_steps - 1)); % check

phi_map = min(max(phi_map, 1), phi_steps);
theta_map = min(max(theta_map, 1), theta_steps);

linearIdx = sub2ind([phi_steps, theta_steps], phi_map, theta_map);

map = accumarray(linearIdx(:), data(validMask), [phi_steps * theta_steps, 1], @sum, 0);
map_count = accumarray(linearIdx(:), 1, [phi_steps * theta_steps, 1], @sum, 0);

map = reshape(map, phi_steps, theta_steps);
map_count = reshape(map_count, phi_steps, theta_steps);


% map(phi_map, theta_map) = map(phi_map, theta_map) + data(validMask);
% map_count(phi_map, theta_map) = map_count(phi_map, theta_map)+1;

% for k = 1:size(u_centered,1)
% 
%     for n = 1:size(u_centered,2) 
% 
%         z_sq = solarRadius^2-((v_centered(k,n))^2+(u_centered(k,n))^2); 
%         if z_sq < 0
%             z_heeq(k,n)= NaN;
%             phi_trans(k,n) = NaN; 
%             theta_trans(k,n) = NaN;
%         else
%             z = sqrt(z_sq);
%             %z_heeq(k,n) = z;
%             phi_trans(k,n) = atan(v_centered(k,n)/sqrt(((z)^2+(u_centered(k,n))^2)));
%             theta_trans(k,n) = atan2(u_centered(k,n),z);
% 
%             %phi_map = round((pi/2+phi_trans(k,n))/phi_grad);
% 
%             phi_map = round(1 + (phi_trans(k,n) + pi/2)/pi *(phi_steps - 1));
%             %The array value each transformed co-ord falls into on the map array
% 
% 
%             %theta_map = (pi/2+theta_trans(k,n))/theta_grad;
% 
%             theta_map = round(1 + (theta_trans(k,n) + pi/2)/pi *(theta_steps - 1));
% 
% 
%             map(phi_map, theta_map) = map(phi_map, theta_map) + data(k,n);
%             map_count(phi_map, theta_map) = map_count(phi_map, theta_map)+1;
%         end
%     end 
% end

% for i = 1:phi_steps
%     for j = 1:theta_steps
%         map_norm(i,j)= map(i,j)/(max(1,map_count(i,j)));
%     end
% end
map_norm = map./max(1,map_count);


% Need to add a thing that maps as it goes and normalises the amount of
% plots per space as this will reduce plotting lag. 

% max(1,map_count) when doing the map normalisation as this avoids dividing
% by zero errors.
mask = (map_norm == 0);
map_norm = regionfill(map_norm, mask);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!! maybe add option here to display image for
% debugging
% figure(1); imagesc(map_norm); colorbar; set(gca, "YDir", "normal");


%scatter(phi_trans(:) , theta_trans(:), [], intensity);
%z_bins = 0:5e14:5e17;
% hist(z_trans(:),z_bins)


end