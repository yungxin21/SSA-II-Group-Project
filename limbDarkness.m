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
rx=xmax-xmin+1;

ymin=find(img(:,2048)/1000>=1,1,"first");
ymax=find(img(:,2048)/1000>=1,1,"last");
ry=ymax-ymin+1;

r=floor(mean([rx,ry]));

% tests conclude that rows 2047-2049 will present the same values of r.

% now to correct for limb darkening.
% we need to take the averages for the values in the middle, to get an idea
% of how much correcting we need to do towards the edges.

% we'll start simple. Assume that the Sun is a perfect sphere, and we can
% relate two variables: the difference between the value of the surface we
% want to look at, and the surface at the center (set D=0), as well as the
% difference between the intensity value at the specified point, and the
% averaged values of the center.

% we can map out a hemisphere with the function f(x,y)=1-sqrt(1-x^2-y^2), 
% setting the radius of the sun to average of rx and ry.

% an issue with this is it is taking values outside of the disc (since it
% includes the corners of xmin,ymin, and xmax,ymax).
% might be necessary to change to polar coordinates, or figure out another
% method. The values below does not work.

x=xmin:1:xmax-1;
y=ymin:1:ymax;
Z=r-sqrt(r^2-x.^2-y.^2);