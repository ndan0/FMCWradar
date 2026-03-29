% Description: detect pixels in a 2D RDM using an improved 2D CFAR algorithm
% Inputs:
%               image       Z(kv, kr) magnitude of 2D RDM
%                           rows: kv
%                           cols: kr
%               Pfa         Prob of false alarm
%               window      rectangular window size, ie. [3, 1] (rows, col)
%
% Outputs:
%               detects     1's based indices of image of detections
%                           col 1: kv
%                           col 2: kr

function detects = cfar_detection(image, Pfa, window)
    if (or(mod(window(1), 2) == 1, mod(window(2), 2) == 1) )
        error('window size must be odd');
    end

    threshold = prod(window) *(Pfa^(-1/(prod(window))) -1);
    
    filt = ones(size_window)/(prod(window)-1);
    filt((window(1)+1)/2, (window(2)+1)/2) = 0;

    if (0)
        U = conv2(image, filt, "same");

        detects = find(threshold*U > image);
    else
        % Improved CA CFAR detector
        % FIXME: should we be replicating instead of zero padding?
        %        or maybe, the filter shouldn't sum values in the pad regions.
        %padded_image = zeros(size(image,1) + window(1)-1, size(image,2) + window(2)-1);
        r_offset = (window(1)+1)/2 - 1;
        c_offset = (window(2)+1)/2 - 1;
        padded_image = padarray(image, [r_offset, c_offset], 'replicate');
        detects = [];

        for r = 1:size(image,1)
            for c = 1:size(image,2)
                U = sum(sum(padded_image(r:r+window(1)-1, c:c+window(2)-1) .* filt));
                if (image(r,c) > threshold*U)
                    detects = [detects; r, c];
                    padded_image(r+r_offset, c+c_offset) = U; % "condB" if pixel is a detect, use it's noise estimate to calculate U for other pixels
                end
            end
        end
    end

end