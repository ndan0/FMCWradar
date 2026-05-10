function [detects, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window)
%CFAR_DETECTION_PAPER
% Paper-style improved 2D CFAR detector for fused RDM Z.
%
% Inputs:
%   Z       : fused 2D range-Doppler map
%             rows = range bins, columns = Doppler bins
%   Tscale  : threshold scaling factor, S = Tscale * U
%   window  : odd-sized rectangular CFAR window [rows, cols]
%
% Outputs:
%   detects : [Ndet x 2], each row = [range_bin, doppler_bin]
%   detMask : logical detection mask
%   Umap    : local noise estimate map
%   Smap    : threshold map

    if numel(window) ~= 2
        error('window must be [rows, cols].');
    end

    if mod(window(1),2) == 0 || mod(window(2),2) == 0
        error('window size must be odd.');
    end

    r_offset = (window(1)-1)/2;
    c_offset = (window(2)-1)/2;

    % Reference-cell averaging filter excluding CUT
    filt = ones(window);
    filt(r_offset+1, c_offset+1) = 0;
    filt = filt / sum(filt(:));

    % Replicate padding for boundary handling
    padded_Z = padarray(Z, [r_offset, c_offset], 'replicate');

    detMask = false(size(Z));
    Umap    = zeros(size(Z));
    Smap    = zeros(size(Z));
    detects = [];

    for r = 1:size(Z,1)
        for c = 1:size(Z,2)

            local_win = padded_Z(r:r+2*r_offset, c:c+2*c_offset);

            U = sum(local_win .* filt, "all");
            S = Tscale * U;

            Umap(r,c) = U;
            Smap(r,c) = S;

            cutVal = padded_Z(r+r_offset, c+c_offset);

            if cutVal > S
                detMask(r,c) = true;
                detects = [detects; r, c]; %#ok<AGROW>

                % Improved CFAR behavior from paper:
                % replace detected target cell by local noise estimate
                % before continuing the scan.
                padded_Z(r+r_offset, c+c_offset) = U;
            end
        end
    end
end