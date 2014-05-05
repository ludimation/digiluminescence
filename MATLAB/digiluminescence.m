function [ ...
    output_C_all, ...
    output_D_all, ...
    output_cleanPlate, ...
    output_uMasks_all, ...
    output_j_features, ...
    output_denseCorr_all, ...
    output_denseCorr_masked_all, ...
    output_denseCorr_multiframe_all, ...
    output_digiLum_all, ...
    output_grid_all ...
] ...
    = digiluminescence( ...
        data_C_all, ...
        data_D_all, ...
        data_joint_positions_all, ...
        data_timestamps, ...
        data_mask_thresh, ...
        data_calcDenseCorr ... 
    )
% function [ output_C_all, output_cleanPlate, output_uMasks_all, output_j_features, output_denseCorr_all, output_digiLum_all, output_grid_all, output_denseCorr_masked_all, output_denseCorr_multiframe_all ] = ...
%     digiluminescence(data_C_all, data_D_all, data_joint_positions_all, data_timestamps, data_mask_thresh, data_calcDenseCorr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% TODO: Later 
%     - Generalize this into a GUI that allows user to select data? 
%     - Check formatting of all inputs to make sure this executes properly
%     - Check number of inputs / outputs [nargoutchk/narginchk(minargs,
%     maxargs)]

% TODO: Now
%   - implement with Vivian's data or re-record her data
%       - save out from openframeworks
%       - load into MATLAB
%       - re-run effect rendering on new data
%   - add points that remain still at the corners
%   - composite video to display at end of semester show
%       - compose audio track
%   - post / burn to CD
%   - complete final write-up
%   - update write up with images from vivian's data?

% Start timer
fprintf('====\n');
fprintf('Digiluminescence :: Executing\n');

%% Handle default arguments
tic
fprintf('====\n');
fprintf('Handling default arguments \n');
for i = 1 % For loop is for code collapsing only (so I don't have to look at these)

    % set default value for inputs
    if( nargin < 5 )
        data_mask_thresh = 256;
    end
    if( nargin < 6 )
        data_calcDenseCorr = false;
    end
   
end

% clean up
clear i

% print time
toc

%% Initialize variables
tic
fprintf('----\n');
fprintf('Initializing variables \n');

ui8_max = double(intmax('uint8'));
ui8_hlf = double(round(intmax('uint8')/2));
i16_max = double(intmax('int16'));
i16_2_ui8 = double(2^7);

[ n_width, n_height, ...
  n_chans, n_frames ]   = size(     data_C_all                      );
n_joints                = size(     data_joint_positions_all    , 1 );
% n_frames                = length(   data_timestamps                 );

% print time
toc

%% Preallocate output values
tic
fprintf('----\n');
fprintf('Preallocating output values \n');

output_C_all                    = zeros(    size(data_C_all)                    , 'uint8'   );
output_cleanPlate               = zeros(    size(data_D_all(:,:,1))             , 'int16'   );
output_uMasks_all               = zeros(    size(data_D_all)                    , 'int16'   );
output_denseCorr_all            = ones(     size(data_C_all)                    , 'int16'   ) * double(ui8_hlf);
output_denseCorr_masked_all     =           output_denseCorr_all                            ;
output_denseCorr_multiframe_all = zeros(    size(output_denseCorr_all)          , 'int16'   );
output_digiLum_all              = zeros(    size(data_C_all)                    , 'uint8'   );
output_grid_all                 = zeros(    size(data_C_all)                    , 'uint8'   );

% print time
toc

%% Create a clean plate of the environment
tic
fprintf('----\n');
fprintf('Creating clean plate for depth data \n');

% out_D_cPlate : an image that represents the background of the scene,
%   which is essentially something that approximates the maximum depth
%   values found in all frames of the depth data... while sticking close to
%   the median and mean to minimize noise when calculating user masks later
%
% To write a new cleanplates:
%   1) trash the clean plate image "test_02_Depth_cPlate_all.png" and,
%   2) run this section  by itself, and it will calculate a new one and
%       save it out in the project path

tmp_recalc_cleanPlate = false;

if (...
        exist('test_02_Depth_cPlate_all.png'    , 'file') == 2 ...
    &&  exist('test_02_Depth_cPlate_max.png'    , 'file') == 2 ...
    &&  exist('test_02_Depth_cPlate_mean.png'   , 'file') == 2 ...
    &&  exist('test_02_Depth_cPlate_median.png' , 'file') == 2 ...
    )
    % load it
    output_cleanPlate = int16(imread('test_02_Depth_cPlate_all.png')) * i16_2_ui8;
    
else % otherwise, create a new one
    tmp_recalc_cleanPlate = true;
end

% calculate max
tic; fprintf(' - tmp_output_cleanPlate_max - ');
    tmp_output_cleanPlate_max = max(data_D_all,[],3); % fast (close second)
toc

% use max to create clean_Data_D_all
tic; fprintf(' - cleaning output_D_all - ');
    % create clean depth images by putting cPlate in areas that have no value
    tmp_inds_positive               = data_D_all > -8;
    output_D_all                    = repmat(tmp_output_cleanPlate_max, [1,1,size(data_D_all,3)]);
    output_D_all(tmp_inds_positive) = data_D_all(tmp_inds_positive);
toc

if tmp_recalc_cleanPlate
    % calculate means and medians from clean depth data
    tic; fprintf(' - tmp_output_cleanPlate_mean - ');
        tmp_output_cleanPlate_mean = mean(output_D_all,3); % fastest
    toc    
    tic; fprintf(' - tmp_output_cleanPlate_median - ');
        tmp_output_cleanPlate_median = median(output_D_all,3); % much slower
    toc
    % use their average (mean) for the mask to mimize noise in the user
    % mask calcs later on
    tmp_output_cleanPlate_all = int16([...
                                tmp_output_cleanPlate_max   , ...
                                tmp_output_cleanPlate_mean  , ...
                                tmp_output_cleanPlate_median  ...
                            ]);
    tmp_output_cleanPlate_all = reshape(tmp_output_cleanPlate_all, [n_width, n_height, 3]);
    %     output_cleanPlate = mean(tmp_output_cleanPlate_all, 3);
    output_cleanPlate = tmp_output_cleanPlate_max   ;
    %     output_cleanPlate = tmp_output_cleanPlate_mean  ;
    %     output_cleanPlate = tmp_output_cleanPlate_median;
    % save out all files for comparison
    imwrite( uint8( output_cleanPlate               / i16_2_ui8 ), 'test_02_Depth_cPlate_all.png'    );
    imwrite( uint8( tmp_output_cleanPlate_max       / i16_2_ui8 ), 'test_02_Depth_cPlate_max.png'    );
    imwrite( uint8( tmp_output_cleanPlate_mean      / i16_2_ui8 ), 'test_02_Depth_cPlate_mean.png'   );
    imwrite( uint8( tmp_output_cleanPlate_median    / i16_2_ui8 ), 'test_02_Depth_cPlate_median.png' );
end


% clean up
clear tmp_*

% print time
toc

%% Create user masks using difference matting against environment clean plate
tic
fprintf('----\n');
fprintf('Creating user masks \n');

% out_uMasks_all holds one mask for each frame that represents the part of
% each depth image that varies from the previously calculated cleanplate,
% out_D_cPlate

% calculate masks for areas with large differences between the clean Depth
% image and the cleanPlate
tmp_D_all_diff                  = abs(output_D_all - repmat(output_cleanPlate, [1,1,size(data_D_all,3)]));
tmp_inds_BG                     = abs(tmp_D_all_diff) < data_mask_thresh;
output_uMasks_all               = output_D_all;
output_uMasks_all(tmp_inds_BG)  = i16_max;
% invert the mask so it can be used as a scalar later
output_uMasks_all = i16_max - output_uMasks_all;

% clean up
clear tmp_*

% print time
toc

%% Calculate projective joint positions for all frames
tic
fprintf('----\n');
fprintf('Calculating projective joint positions for all frames \n');

% Re-shape joint_positions_all and split it up so it can be fed into
% projectiveTransformation (x, y, z, to)
%  TODO: 
%    - would be nice if the projectiveTransformation function only took two
%    args (p, to) and did the reshaping on its own
j_pos_all_reshaped = reshape(permute(data_joint_positions_all, [2 1 3]), 3, n_frames*n_joints);
j_pos_all_x = j_pos_all_reshaped(1,:);
j_pos_all_y = j_pos_all_reshaped(2,:);
j_pos_all_z = j_pos_all_reshaped(3,:);

% feed individual x, y, and z lists into projectiveTransformation (x, y, z, to)
[ j_pos_all_x_prjctd, j_pos_all_y_prjctd, j_pos_all_z_prjctd ] = ...
    projectiveTransformation(j_pos_all_x, j_pos_all_y, j_pos_all_z, 'projective');

% build j_pos_all_projective from reshaped return vaules
j_pos_all_reshaped_projective = [ j_pos_all_x_prjctd; j_pos_all_y_prjctd; j_pos_all_z_prjctd ];
j_pos_all_projective = permute(reshape(j_pos_all_reshaped_projective, 3,n_joints,n_frames), [2 1 3]);

% Create second array for easy grabbing of one set of joints from the
% previous frame (frame 1 is compared against frame n which doesn't quite
% make sense, but I'd just like to see what happens)
j_pos_all_projective_start = circshift(j_pos_all_projective, [0,0,1]);

% Create list of feature correspondences in a shape that thin plate dense
% corespondence function requires
%    TODO:
%    - create a 3D thin plate function and update this to include z
j_startingFeatures  = permute(j_pos_all_projective_start(   :, 1:2, :)  , [2,1,3]);
j_endingFeatures    = permute(j_pos_all_projective(         :, 1:2, :)  , [2,1,3]);
% switch x, y positions
j_startingFeatures  = circshift(j_startingFeatures  , 1);
j_endingFeatures    = circshift(j_endingFeatures    , 1);
% flip x values and translate them the width of the screen
output_j_features = [ j_endingFeatures; j_startingFeatures];

% TODO: 
%    - add additional feature points along the lines of larger key limbs \
%    torso
%    - might also need a few additional features outside bounds of figure
%    positions to "tack" down the surrounding frame of the image so it
%    doesn't warp around so much

% cleanup
clear j_pos_all_reshaped j_pos_all_reshaped_projective
clear j_pos_all_x j_pos_all_y j_pos_all_z
clear j_pos_all_x_prjctd j_pos_all_y_prjctd j_pos_all_z_prjctd 
clear j_pos_all_projective j_pos_all_projective_start 

% print time
toc

%% Draw grids
tic
fprintf('----\n');
fprintf('Drawing grids \n');

% drawGrid (I, spcGrid, spcPoints, color)
output_grid_all     = zeros(    size(data_C_all)                    , 'uint8'   );
output_grid_all     = drawGrid(output_grid_all     , 32, 1, [ui8_max, ui8_max, 0]); % yellow

% print time
toc

%% Draw joints and limbs
tic
fprintf('----\n');
fprintf('Drawing debug joints \n');

% drawPoints(points, I, size, color)
output_grid_all = drawPoints(j_startingFeatures, output_grid_all, 16, [ui8_max, ui8_max, 0] ); % yellow
%TODO: draw limbs in white
% drawLimbs(j_startingFeatures, out_grid_all, 4, [ui8_max, ui8_max, 0] ); % yellow

% print time
toc


%% Calculate dense correspondence fields frame by frame
tic
fprintf('----\n');
fprintf('Calculating dense correspondence fields frame by frame \n');

% create dense correspondence fields one frame at a time
if data_calcDenseCorr
    for i_frame = 1:n_frames
        tic
        fprintf([' - frame ' num2str(i_frame) ' of ' num2str(n_frames) ' - ']);
            [ denseCorr, grid ] = thin_plate_denseCorrespondence(...
                output_j_features(:, :, i_frame), ...
                output_grid_all(:,:,:, i_frame) ...
            );
            output_denseCorr_all(:,:,:, i_frame) = denseCorr;
            output_grid_all(:,:,:, i_frame) = grid;
        % print time
        toc
    end
else
    % print time
    toc
end

%% Draw debug joints and reference grid over warped grid
tic
fprintf('----\n');
fprintf('Drawing debug joints and reference grid \n');

% draw grid again sparser and blue to show where it started : drawGrid (I, spcGrid, spcPoints, color)
output_grid_all = drawGrid(output_grid_all, 32, 4, [0 , ui8_max , ui8_max]); % cyan
% draw ending joints in magenta (small radius -- warped joints will appear as large yellow dots)
output_grid_all = drawPoints(j_endingFeatures      , output_grid_all, 8, [ui8_max , 0         , ui8_max] ); % magenta
% draw starting joints in blue (smaller radius)
output_grid_all = drawPoints(j_startingFeatures    , output_grid_all, 6, [0       , ui8_max   , ui8_max] ); % cyan

% cleanup
clear j_pos_all_projective2
clear iterator denseCorr grid

% print time
toc

%% Create digiluminescence effect from dense correspondence fields
tic
fprintf('----\n');
fprintf('Creating digiluminescence effect frame by frame \n');

% Mask dense correspondence field with movement mask
output_denseCorr_masked_all = int16( ...
            double(output_denseCorr_all) ...
        .*  double(repmat(permute(output_uMasks_all, [1,2,4,3]), [1,1,3,1])) ...
        /   double(i16_max) ...
    );
toc

% create multiframe field from faded, circshifted versions of masked denseCorr
output_denseCorr_multiframe_all = zeros(    size(output_denseCorr_all)          , 'int16'   );
i_tmp_max = min(16, n_frames);
for i_tmp = i_tmp_max:-1:1
    tic
    fprintf(['- Circshift ' num2str(i_tmp_max-i_tmp) ' of ' num2str(i_tmp_max) '\n']);

    tmp_inds = repmat(permute(output_uMasks_all, [1,2,4,3]), [1,1,3,1]) > 2^4;
    
    % circshift iteration-1 so it starts with current frame
    output_denseCorr_multiframe_all(tmp_inds) = ...
        circshift(output_denseCorr_masked_all(tmp_inds), [1,1,1,i_tmp - 1] ) ...
        * ((-((i_tmp-1)/i_tmp_max)+1)^(1/2)) ... y = sqrt(-x+1) easing
        .../iteration_max ... scaled based on maximum number of iterations
        ;
    %print time
    toc
end
% for i_tmp = 1:i_tmp_max
%     % circshift iteration-1 so it starts with current frame
%     output_denseCorr_multiframe_all = output_denseCorr_multiframe_all ...
%         + circshift(output_denseCorr_masked_all, [1,1,1,i_tmp - 1] ) ...
%         * ((-((i_tmp-1)/i_tmp_max)+1)^(1/2)) ... y = sqrt(-x+1) easing
%         /i_tmp_max ... scaled based on maximum number of iterations
%         ;
% end

imshow(uint8(output_denseCorr_multiframe_all(:,:,:,1)));

% clean up
clear i_* tmp_*

% print time
toc

%% Plot lines along dense corresponence field
tic
fprintf('- Plotting lines along dense corresponence field - \n');

% TODO: draw lines from 'output_denseCorr_multiframe_all' into 'output_digiLum_all'
output_digiLum_all = zeros( size(data_C_all), 'uint8');

for i_frame = 1:n_frames
    tic
    fprintf([' - frame ' num2str(i_frame) ' of ' num2str(n_frames) ' - \n   - ']);

    % grab x, y, and u, v for each frame
    tmp_img_source = output_denseCorr_multiframe_all(:,:,:,i_frame);
    [tmp_p1_x, tmp_p1_y] = meshgrid(1:size(tmp_img_source,1), 1:size(tmp_img_source, 2));
    tmp_p1_x = double(reshape(tmp_p1_x, [numel(tmp_p1_x), 1]));
    tmp_p1_y = double(reshape(tmp_p1_y, size(tmp_p1_x)));
    tmp_p_dx = double(reshape(tmp_img_source(:,:, 2)', size(tmp_p1_x))) * 2^2;
    tmp_p_dy = double(reshape(tmp_img_source(:,:, 1)', size(tmp_p1_x))) * 2^2;
    tmp_p2_x = tmp_p1_x + tmp_p_dx;
    tmp_p2_y = tmp_p1_y + tmp_p_dy;
    
%     p1_x = [120, 240, 360]'
%     p1_y = [160, 320, 480]'
%     p_dx = [-20, 10, -40]'
%     p_dy = [-20, -40, 10]'
%     p2_x = p1_x + p_dx;
%     p2_y = p1_y + p_dy;

    % create logical array of lines whose magnitude is very small
    tmp_lines_to_skip = double(tmp_p_dx.^2 + tmp_p_dy.^2).^(1/2) < 2^4;
    
    % capture instance of the image to be drawn into
    tmp_img_target = output_digiLum_all(:,:,:,i_frame);

%     % debug
%     imshow(tmp_img_target)
    
    for i_line = 1:numel(tmp_p1_x)
        if tmp_lines_to_skip(i_line)
            continue
        end
        fprintf(['l ' num2str(i_line) ' of ' num2str(numel(tmp_p1_x)) ' - ']);
        
        %// Draw lines from p1 to p2 on matrix tmp_ing_target
        % x = p1(1):p2(1)
        tmp_line_x =    min(tmp_p1_x(i_line), tmp_p2_x(i_line))...
                    : ...
                    max(tmp_p1_x(i_line), tmp_p2_x(i_line));
        tmp_line_x = reshape(tmp_line_x, [numel(tmp_line_x), 1]);

        % round((x - p1(1)) * p_dy / p_dx + p1(2));
        tmp_line_y = round(...
                              (tmp_line_x - tmp_p1_x(i_line)) ...
                            * tmp_p_dy(i_line) ...
                            / tmp_p_dx(i_line) ...
                            + tmp_p1_y(i_line) ...
                    ); ...
                  
        tmp_line_xxx = repmat(tmp_line_x, [3, 1]);
        tmp_line_yyy = repmat(tmp_line_y, [3, 1]);
        tmp_chan = reshape(...
                            repmat(...
                                [1,2,3], ...
                                [numel(tmp_line_x), 1]), ...
                            size(tmp_line_xxx) ...
                        );
        tmp_rgb = reshape(...
                            repmat(...
                                    [0, round(ui8_hlf/2), ui8_max], ... TODO: make this color user-settable
                                    size(tmp_line_x) ...
                                ), ...
                            size(tmp_line_xxx)...
                        ) ;
                    
        % clamp tmp_line_xxx and tmp_line_yyy to indexes withing image bounds
        tmp_line_xxx = max(min(tmp_line_xxx,size(tmp_img_target, 1)),1);
        tmp_line_yyy = max(min(tmp_line_yyy,size(tmp_img_target, 2)),1);
        
        % re-cast as doubles (otherwise, sub2ind() will throw a data type
        % error... not sure why)
        tmp_line_xxx    = double(tmp_line_xxx);
        tmp_line_yyy    = double(tmp_line_yyy);
        tmp_chan        = double(tmp_chan);
        
        % draw color into the indexes
        % m(sub2ind(size(m), y, x, channel)) = 1;
        tmp_inds = sub2ind(size(tmp_img_target), tmp_line_xxx, tmp_line_yyy, tmp_chan);
        tmp_img_target(tmp_inds) = tmp_rgb;
    end
    
    % Draw lines along vectors of the field into the digiluminescence
    % effect output
    output_digiLum_all(:,:,:,i_frame) = tmp_img_target;
    
    % debug
    imshow(tmp_img_target);
    
    % print time
    fprintf('\n')
    toc
end

% apply digilum effect to C_all
output_C_all = data_C_all + output_digiLum_all;

% clean up
clear i_* tmp_*

% print time
toc

%% Save out test files
fprintf('----\n');
fprintf('Saving out some test files \n');

%% output_C_all
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_C_all']);
    % IMG data must be of one of the following classes: double, single, uint8
    tmp_output_C_all                    = uint8(output_C_all                            )               ;
% print time
toc

tic
fprintf(['     > saving image - ']);
    imwrite(tmp_output_C_all(:,:,:,1)                   ,[ 'test_07_Color_ouput.png'            ]);
% print time 
toc

tic
fprintf([' - videos - output_C_all - ']);
    writerObj = VideoWriter(['test_07_Color_ouput.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_C_all)
    close(writerObj);
% print time
toc

%% output_D_all

tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_D_all']);
    tmp_output_D_all = uint8(output_D_all / i16_2_ui8 );
    % must have a [w,h,bitDepth, frames] array for video file writing
    tmp_output_D_all = permute(tmp_output_D_all, [1,2,4,3]);
% print time
toc
tic
fprintf(['     > saving images - ']);
    % save out images
    imwrite(tmp_output_D_all(:,:,:,1), [ 'test_02_Depth_ouput.png' ]);
% print time
toc

% output_D_all
tic
fprintf(['     > saving video - output_D_all - ']);
    writerObj = VideoWriter(['test_02_Depth_ouput.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_D_all)
    close(writerObj);
% print time
toc


%% output_cleanPlate
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_cleanPlate']);
    tmp_output_cleanPlate = uint8(output_cleanPlate / i16_2_ui8 );
% print time
toc

tic
fprintf(['     > saving image - ']);
    imwrite(tmp_output_cleanPlate, [ 'test_02_Depth_cPlate.png' ]);
% print time 
toc

%% output_uMasks_all
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_uMasks_all']);
    % IMG data must be of one of the following classes: double, single, uint8
    tmp_output_uMasks_all = uint8(output_uMasks_all / i16_2_ui8 );
    % must have a [w,h,bitDepth, frames] array for video file writing
    tmp_output_uMasks_all = permute(tmp_output_uMasks_all, [1,2,4,3]);
% print time
toc
tic
fprintf(['     > saving image - ']);
    imwrite(tmp_output_uMasks_all(:,:,:,1)              ,[ 'test_03_uMask.png'                  ]);
% print time 
toc


tic
fprintf(['     > saving video - output_uMasks_all - ']);
    writerObj = VideoWriter(['test_03_uMask.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_uMasks_all)
    close(writerObj);
% print time
toc

%% output_denseCorr_all
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_denseCorr_all, output_denseCorr_masked_all, output_denseCorr_multiframe_all']);
    % set up scale and translate variables used later on for dense correspondence visibility
    if data_calcDenseCorr
        % turn up coefficient to increase visualcontrast in dense
        % correspondence data
        dc_scale = double(2^3); 
    else
        % if dense correspondence data is not being generated, it is already
        % maximal by default, so you don't need to scale it
        dc_scale = double(2^0);
    end
    dc_offset = double(0); % ui8_hlf; % 
    % IMG data must be of one of the following classes: double, single, uint8
    tmp_output_denseCorr_all            = uint8(output_denseCorr_all                    )   + dc_offset ;
    tmp_output_denseCorr_masked_all     = uint8(output_denseCorr_masked_all             )   + dc_offset ;
    tmp_output_denseCorr_multiframe_all = uint8(output_denseCorr_multiframe_all         )   + dc_offset ;
    % scale and translate dense correspondence for visibility
    tmp_output_denseCorr_all            = uint8((double(tmp_output_denseCorr_all            ) - dc_offset) * dc_scale + dc_offset);
    tmp_output_denseCorr_masked_all     = uint8((double(tmp_output_denseCorr_masked_all     ) - dc_offset) * dc_scale + dc_offset);
    tmp_output_denseCorr_multiframe_all = uint8((double(tmp_output_denseCorr_multiframe_all ) - dc_offset) * dc_scale + dc_offset);
    %clean up
    clear dc_*   
% print time
toc

tic
fprintf(['     > saving images - ']);

    imwrite(tmp_output_denseCorr_all(:,:,:,1)           ,[ 'test_04_denseCorr.png'              ]);
    imwrite(tmp_output_denseCorr_masked_all(:,:,:,1)    ,[ 'test_06_denseCorr_masked.png'       ]);
    imwrite(tmp_output_denseCorr_multiframe_all(:,:,:,1),[ 'test_06_denseCorr_multiframe.png'   ]);

% print time
toc

tic
fprintf(['     > saving video - output_denseCorr_all - ']);
    writerObj = VideoWriter(['test_04_denseCorr.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_all)
    close(writerObj);
% print time
toc

% output_denseCorr_masked_all
tic
fprintf(['     > saving video - output_denseCorr_masked_all - ']);
    writerObj = VideoWriter(['test_06_denseCorr_masked.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_masked_all)
    close(writerObj);
% print time
toc

% output_denseCorr_multiframe_all
tic
fprintf(['     > saving video - output_denseCorr_multiframe_all - ']);
    writerObj = VideoWriter(['test_06_denseCorr_multiframe.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_multiframe_all)
    close(writerObj);
% print time
toc

%% output_grid_all
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_grid_all']);
    tmp_output_grid_all = uint8(output_grid_all );
% print time
toc

tic
fprintf(['     > saving image - ']);
    imwrite(tmp_output_grid_all(:,:,:,1), [ 'test_05_grid_warped.png' ]);
% print time 
toc

tic
fprintf(['     > saving video - output_grid_all - ']);
    writerObj = VideoWriter(['test_05_grid_warped.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_grid_all)
    close(writerObj);
% print time
toc

%% output_digiLum_all
tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - output_digiLum_all']);
    tmp_output_digiLum_all = uint8(output_digiLum_all );
% print time
toc

tic
fprintf(['     > saving image - ']);
    imwrite(tmp_output_digiLum_all(:,:,:,1), [ 'test_06_digiLum.png' ]);
% print time 
toc

tic
fprintf(['     > saving video - output_digiLum_all - ']);
    writerObj = VideoWriter(['test_06_digiLum.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_digiLum_all)
    close(writerObj);
% print time
toc

%% data_C_all, data_D_all

tic
fprintf(['     ---- \n'])
fprintf(['     > reformatting data - data_C_all, data_D_all']);
    tmp_data_C_all = uint8(data_C_all );
    tmp_output_D_all = uint8(data_D_all / i16_2_ui8 );
    % must have a [w,h,bitDepth, frames] array for video file writing
    tmp_output_D_all = permute(tmp_output_D_all, [1,2,4,3]);
% print time
toc
tic
fprintf(['     > saving images - ']);
    % save out images
    imwrite(tmp_data_C_all(:,:,:,1), [ 'test_01_Color_data.png' ]);
    imwrite(tmp_output_D_all(:,:,:,1), [ 'test_02_Depth_data.png' ]);
% print time
toc

% data_C_all
tic
fprintf(['     > saving video - data_C_all - ']);
    writerObj = VideoWriter(['test_01_Color_data.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_data_C_all)
    close(writerObj);
% print time
toc

% data_D_all
tic
fprintf(['     > saving video - data_D_all - ']);
    writerObj = VideoWriter(['test_02_Depth_data.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_D_all)
    close(writerObj);
% print time
toc

%% clean up memory
clear writerObj tmp_*


%% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n'); 
%TODO: figure out a way to print all elapsed time for this function
fprintf('====\n');
end