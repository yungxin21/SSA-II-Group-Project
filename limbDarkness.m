function limbDarkness(fits_input,fits_output,plot_spacing)
    copyfile(fits_input,"corrected.fits");
    img = fitsread(fits_input, "image");

    % Magnetogram option (in separate branch)

    % AIA EUV option
    % euv_lev1 = fitsread(fits_input, "image");
    % figure(3); imagesc(euv_lev1); colormap("bone"); colorbar;

    % 1. Identify the solar discs

    info = fitsinfo(fits_input);
    keywords = info.Image.Keywords;
    valIdx = strcmp(keywords(:,1), 'CRPIX1');
    valIdx2 = strcmp(keywords(:,1), 'CRPIX2');
    center_x = round(keywords{valIdx, 2});
    center_y = round(keywords{valIdx2, 2});
    
    % find the first and last index with positive values larger than 1000
    xmin=find(img(center_x,:)/1000>=1,1,"first"); 
    xmax=find(img(center_x,:)/1000>=1,1,"last");
    diameter_x=xmax-xmin+1;
    radius_x=diameter_x/2;
    
    ymin=find(img(:,center_y)/1000>=1,1,"first");
    ymax=find(img(:,center_y)/1000>=1,1,"last");
    diameter_y=ymax-ymin+1;
    radius_y=diameter_y/2;
    
    radius = floor(mean([radius_x,radius_y])); % the radius of the sun in pixels
    
    % 2. Correct for limb darkening.
    % SDO HMI continuum image wavelength: 617.3 nm (limb darkness is
    % dependent on captured wavelength)
    
    % grid of pixel coordinates
    [X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
    % calculate the radial distance (r) from the center for every pixel
    r = sqrt((X - center_x).^2 + (Y - center_y).^2);
    
    % here mu = cos(theta), where theta is the heliocentric angle
    mu = sqrt(max(0, 1 - (r ./ radius).^2)); 
    % limb darkening coefficients (see references):
    u1 = 0.836; 
    u2 = -0.204;
    
    % equals ~1 at the around the centre (i.e, no correction) and near 0 around
    % the edges (a lot of correction)
    limb_fit = 1 - u1 - u2 + u1 * mu + u2 * (mu).^2;
    
    % need to change this to work with sunspot.m
    solar_disk_mask = (r <= radius); % basically just 0 where the sun isnt
    % and 1 where the sun is, so correction is only applied on sun and not background
    % solar_disk_mask(solar_disk_mask<0) = NaN;
    
    corrected_img = img;
    % apply the corrected intensity values to the array of 0s and 1s
    corrected_img(solar_disk_mask) = img(solar_disk_mask) ./ limb_fit(solar_disk_mask);        
    figure(1); imagesc(corrected_img); colormap; colormap("bone"); colorbar;
    corrected_img(corrected_img<0) = NaN;
    % add the same metadata from img to corrected_img.
    fitswrite(corrected_img,fits_output);


    % adjustable intensity plot with normalized radii.
    if (plot_spacing > 0)
        correctedModes = [];
        originalModes = [];
        % view intensity over r/r_solar plot.
        for i=1:plot_spacing:radius
            ring = (round(r)==round(i));
            correctedModes = [correctedModes, mode(corrected_img(ring))]; 
            originalModes = [originalModes, mode(img(ring))];
        end
        figure(2); plot((1:plot_spacing:radius)/radius,correctedModes); hold on;
        plot((1:plot_spacing:radius)/radius,originalModes); hold off;
        title('Intensity over Radius Plot')
        xlabel('r/r_☉')
        ylabel('intensity')
        xlim([0 1.25])
        ylim([0 55000])
    end
end