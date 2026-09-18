function [maskedImg, mask] = sigmaClip(data, nSigma, maxI)
    % Make all negative values to be NaN for ease of finding 0 in data
    maskedImg = logical(data);
    maskedImg(:) = 0;
    data(data<0) = NaN;
    % Dedicated map of where NaN's are
    nanMask = isnan(data);
    % Start mask (true = valid, false = clipped)
    mask = true(size(data));
    for i=1:maxI
        % Make current unclipped values as 1D vector
        validData = data(mask & ~nanMask);
        if isempty(validData), break; end
        % Running statistics on valid data and ommiting clipped values
        centreVal = mode(validData);
        stdVal = std(validData);
        if stdVal == 0, break; end
        % Defining boundaries
        lBound = centreVal - nSigma * stdVal;
        % Use bounds to generate new mask 
        newMask = mask & (data >= lBound);
        if isequal(newMask,mask), break; end
        mask = newMask;
    end
    % Replacing clipped values to keep same matrix dimensions
    maskedImg(~mask & ~nanMask) = 1;
end