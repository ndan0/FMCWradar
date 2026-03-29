% Description: detect pixels in a 2D RDM using an improved 2D CFAR algorithm
% Inputs:
%               image       Z(kv, kr) magnitude of 2D RDM
%                           rows: kv
%                           cols: kr
%               Pfa         Prob of false alarm
%               window      rectangular window size, ie. [3, 1] (cols, rows)
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
    
    %padded_image = zeros(size(image,1) + window(1)-1, size(iamge,2) + window(2)-1);
    filt = ones(size_window)/(prod(window)-1);
    filt((window(1)+1)/2, (window(2)+1)/2) = 0;

    % FIXME: This is just a standard square law CFAR algorithm.
    %       How do we make this the "improved" version in the paper
    %       that has conditioning on if (x1,y1) is a target?
    U = conv2(image, filt, "same");

    detects = find(threshold*U > image);

end