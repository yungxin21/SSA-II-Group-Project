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
rx=xmax-xmin;

ymin=find(img(:,2048)/1000>=1,1,"first");
ymax=find(img(:,2048)/1000>=1,1,"first");
ry=ymax-ymin;

% tests conclude that rows 2047-2049 will present the same values of r.

% now to correct for limb darkening.


