%% DRAWPOINTS method
function [ I ] = drawPoints(p_array, I, sz_draw, c_rgb)
% draws squares around point of the specified size and color into the image
% provided

    % check inputs
    narginMin = 2;
    narginMax = 4;
    narginchk(narginMin, narginMax);
    if nargin < 3
        sz_draw = 8;
    end
    if nargin < 4
        c_rgb = white(1) * intmax('int16');
    end
    
    sz_I = size(I);
    I_px = -ones(sz_I, 'int16');
    % indeces must be rounded doubles / integers %TODO: might also want to
    % check values to make sure none of them are 0, nor max width and
    % height of I
    p_array = round(p_array);
    p_size = size(p_array);
    p_n = prod(p_size);
    p_h = p_size(1);
    p_w = p_size(2);
    p_frames = p_size(3);
    
    % build arrays of x coords, y coords, channel numbers, and frame
    % numbers starting with x so we can save out some properties that
    % should remain constant throughout
    x_array = reshape(p_array(1,:,:), [p_n/p_h,1,1]);
    sz_x_array = size(x_array);
    % remaining arrays
    y_array = reshape(p_array(2,:,:), sz_x_array)       ;
    r_array = ones(sz_x_array, 'double')                ;
    g_array = ones(sz_x_array, 'double')     *2         ; 
    b_array = ones(sz_x_array, 'double')     *3         ; 
    frame_array = double(reshape(repmat(1:p_frames, [p_w,1,1]), sz_x_array));

    % TODO: add additional points to array to account for "sz_draw"
    
    % clamp x and y arrays from 1 to height & width of image
    x_array = max(min(x_array,size(I, 1)),1);
    y_array = max(min(y_array,size(I, 2)),1);
     
%     % TODO: 
%     %   - create a duplicate function for this so that I can clock them
%     %   against each other
%     %   - concatenate all arays so sub2in() only has to be called once
%     %   - create a c_rgb_array that repmats c_rgb correctly so pixel
%     %   coloring call two steps below only has to be called once
%     xxx_array = repmat(x_array, [3, 1]);
%     yyy_array = repmat(y_array, [3, 1]);
%     rgb_array = [r_array; g_array; b_array];
%     frame_rgb_array = repmat(frame_array, [3, 1]);
%     c_rgb_array = reshape(repmat(c_rgb, size(x_array)), rgb_array);
% 
%     point_indexes = sub2ind(sz_I, xxx_array, yyy_array, rgb_array, frame_rgb_array);
%     
%     I_px(point_indexes) = c_rgb_array(point_indexes);
    
    % create one array for each channel that contains the indexes of the
    % pixels to be colored
    point_indexes_r = sub2ind(sz_I, x_array, y_array, r_array, frame_array);
    point_indexes_g = sub2ind(sz_I, x_array, y_array, g_array, frame_array);
    point_indexes_b = sub2ind(sz_I, x_array, y_array, b_array, frame_array);
    
    % color the pixels
    I_px(point_indexes_r) = c_rgb(1);
    I_px(point_indexes_g) = c_rgb(2);
    I_px(point_indexes_b) = c_rgb(3);

    % TODO: replace this (likely the slowest part) by adding points to
    % points array to account for point size before finding indexes
    % start drawing at the top left of the point
    sz_draw_half = round(sz_draw/2);
    I_px = circshift(I_px, [-sz_draw_half,-sz_draw_half]);
    for pos_x = 1:sz_draw
        for pos_y = 1:sz_draw
            I_tmp = circshift(I_px, [pos_x, pos_y]);
            inds = find(I_tmp>-1);
            I(inds) = I_tmp(inds);
        end
    end
end