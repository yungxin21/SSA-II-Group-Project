clear all
% In later versions, the file name will have to be adjusted to read mutliple files.
img = fitsread("hmi.Ld_720s.20190413_210000_TAI.3.linedepth.fits", "image");

% identify the solar disc. We can identify the radius of this disc in terms
% of pixel size by using the image data.

% automate this by finding the first and last value that is not -1000, 
% and that is not pitch black, in the centermost row/column.
% the intensity of the background space given by the FITS file has values
% under 1000, which is why each value is divided by 1000.

xmin=find(img(2048,:)/1000>=1,1,"first");
xmax=find(img(2048,:)/1000>=1,1,"last");
diax=xmax-xmin+1;
rx=diax/2;

ymin=find(img(:,2048)/1000>=1,1,"first");
ymax=find(img(:,2048)/1000>=1,1,"last");
diay=ymax-ymin+1;
ry=diay/2;

dia=floor(mean([diax,diay]));
r=floor(mean([rx,ry]));

% tests conclude that rows 2047-2049 will present the same values of r.

% now to correct for limb darkening.

% we'll start simple. Assume that the Sun is a perfect sphere, and we can
% relate two variables: the difference between the value of the surface we
% want to look at, and the surface at the center (set D=0), as well as the
% difference between the intensity value at the specified point, and the
% averaged values of the center.

% an issue with this is it is taking values outside of the disc (since it
% includes the corners of xmin,ymin, and xmax,ymax).
% might be necessary to change to polar coordinates, or figure out another
% method.

% need to figure out which one takes the divisor. (whether diax or diay
% is larger, and so one of the divs will be 1, while the other 1.00001
% (extra row/col of pixels)
% temp fix is to manually adjust entries to create a perfect circle.

x = xmin:1:xmax-1;
y = ymin:1:ymax;

[X,Y] = meshgrid(x,y);
Zsq = r^2-(X-r).^2-(Y-r).^2;
% remove negative numbers to prevent complex array
Zsq(Zsq<0) = NaN;
Z = r-sqrt(Zsq);

% so the above works, but has values up to 1.8k? So somehow we need to
% factor it down to 1 and above.

% alternate way is to measure the change in the drop of intensity over time
% maybe as a derivative or as a gradient (since it is a surface)
% 