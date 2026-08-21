close all
% In later versions, the file name will have to be adjusted to read multiple files.
img = fitsread("hmi.Ld_720s.20190413_210000_TAI.3.linedepth.fits", "image");
newimg = img;
% figure(2); imagesc(newimg); colormap; colormap("bone"); colorbar;

% 1. Identify the solar disc. 

% automate this by finding the first and last value that is not -1000, 
% and that is not pitch black, in the centermost row/column.
% the intensity of the background space given by the FITS file roughly has values
% under 1000, which is why each value is divided by 1000.

umin=find(img(2048,:)/1000>=1,1,"first");
umax=find(img(2048,:)/1000>=1,1,"last");
diau=umax-umin+1;
ru=diau/2;

vmin=find(img(:,2048)/1000>=1,1,"first");
vmax=find(img(:,2048)/1000>=1,1,"last");
diav=vmax-vmin+1;
rv=diav/2;

dia=floor(mean([diau,diav]));
R = floor(mean([ru,rv])); % the radius of the sun in pixels

% 2. Correct for limb darkening.
% SDO HMI continuum image wavelength: 617.3 nm (this is important for the
% correction factor!). this is because the limb darkening is dependent on
% the wavelength the image was captured at.

center_x = 2048;
center_y = 2048;
% grid of pixel coordinates
[X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
% calculate the radial distance (r) from the center for every pixel
r = sqrt((X - center_x).^2 + (Y - center_y).^2);

% mu = cos(theta), where theta is the heliocentric angle
% also we clamp the values to avoid complex numbers or division by zero outside the limb
mu = sqrt(max(0, 1 - (r ./ R).^2)); 
% limb darkening coefficients (see references):
u1 = 0.836; 
u2 = -0.204;

% equals ~1 at the around the centre (i.e, no correction) and near 0 around
% the edges (a lot of correction)
limb_fit = 1 - u1 - u2 + u1 * mu + u2 * (mu).^2;
 
solar_disk_mask = (r <= R); % basically just 0 where the sun isnt
% and 1 where the sun is, so correction is only applied on sun and not background
% this is kind of redundant because of newimg which removes the background
% but whatever

corrected_img = img;
% apply the corrected intensity values to the array of 0s and 1s
corrected_img(solar_disk_mask) = img(solar_disk_mask) ./ limb_fit(solar_disk_mask);
    
figure(1); imagesc(corrected_img); colormap; colormap("bone"); colorbar;


% view intensity over r/r_solar plot.
% this one is horizontal, from 0 to the left outer edge of the Sun.
figure(2); plot((1:1:ru)/ru,corrected_img(2048, 2041:1:umax)) 
title('Intensity over Radius Plot')
xlabel('r/r_☉') 
ylabel('intensity')
xlim([0 1.5])
ylim([0 40000])

% could maybe figure out a way to do it based on where the cursor is
% hovering over? kinda like in DS9.

