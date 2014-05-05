function [ ...
    output_C_all, ...
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

n_joints            = size(     data_joint_positions_all         , 1         );
n_frames            = length(   data_timestamps                              );

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
%   which is essentially the maximum depth values found in all frames of
%   the depth data
% To write a new cleanplate:
%   1) trash the clean plate image "test_02_Depth_cPlate_all.png" and,
%   2) run this section and it will calculate a new one and save it out in
%       the project path

if (exist('test_02_Depth_cPlate_all.png', 'file') == 2)
    % load it
    output_cleanPlate = int16(imread('test_02_Depth_cPlate_all.png')) * i16_2_ui8;
else
    % create a new one
    output_cleanPlate = max(data_D_all,[],3);
    % and save it out
    tmp_output_cleanPlate = uint8(output_cleanPlate / i16_2_ui8 );
    imwrite(tmp_output_cleanPlate, 'test_02_Depth_cPlate_all.png');
end

% print time
toc

%% Create user masks using difference matting against environment clean plate
tic
fprintf('----\n');
fprintf('Creating user masks \n');

% out_uMasks_all holds one mask for each frame that represents the part of
% each depth image that varies from the previously calculated cleanplate,
% out_D_cPlate

% clean up depth images by putting cPlate in areas that have no value
inds_positive               = data_D_all > -8; 
D_all_cPlate                = repmat(output_cleanPlate, [1,1,size(data_D_all,3)]);
D_all_clean                 = D_all_cPlate;
D_all_clean(inds_positive)  = data_D_all(inds_positive);
% calculate masks for areas with large differences between the clean Depth
% image and the cleanPlate
D_all_diff                  = abs(D_all_clean - D_all_cPlate);
inds_BG                     = abs(D_all_diff) < data_mask_thresh;
output_uMasks_all           = D_all_clean;
output_uMasks_all(inds_BG)  = i16_max;
% invert the mask so it can be used as a scalar later
output_uMasks_all = i16_max - output_uMasks_all;

% clean up
clear inds_positive inds_BG
clear D_all_cPlate D_all_clean D_all_diff

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
% drawLimbs(j_pos_all_projective, out_grid_all, 4, [ui8_max, ui8_max, 0] ); % yellow

% print time
toc


%% Calculate dense correspondence fields frame by frame
tic
fprintf('----\n');
fprintf('Calculating dense correspondence fields frame by frame \n');

% create dense correspondence fields one frame at a time
if data_calcDenseCorr
    for iterator = 1:n_frames
        tic
        fprintf([' - frame ' num2str(iterator) ' of ' num2str(n_frames) ' - ']);
            [ denseCorr, grid ] = thin_plate_denseCorrespondence(...
                output_j_features(:, :, iterator), ...
                output_grid_all(:,:,:, iterator) ...
            );
            output_denseCorr_all(:,:,:, iterator) = denseCorr;
            output_grid_all(:,:,:, iterator) = grid;
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
clear j_pos_all_projective j_pos_all_projective2
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

% create multiframe field from faded, circshifted versions of masked denseCorr
output_denseCorr_multiframe_all = zeros(    size(output_denseCorr_all)          , 'int16'   );
iteration_max = min(16, n_frames);
for iteration = 1:iteration_max
    % circshift iteration-1 so it starts with current frame
    output_denseCorr_multiframe_all = output_denseCorr_multiframe_all ...
        + circshift(output_denseCorr_masked_all, [1,1,1,iteration - 1] ) ...
        * ((-((iteration-1)/iteration_max)+1)^(1/2))...
        .../iteration_max ... y = sqrt(-x+1) easing
        ;
end

% print time
toc

%% Plot lines along dense corresponence field
tic
fprintf('- Plotting lines along dense corresponence field - \n');

% TODO: draw lines from 'output_denseCorr_multiframe_all' into 'output_digiLum_all'
output_digiLum_all = zeros( size(data_C_all), 'uint8');

for iterator = 1:n_frames
    tic
    fprintf([' - frame ' num2str(iterator) ' of ' num2str(n_frames) ' - \n   - ']);

    % grab x, y, and u, v for each frame
    tmp_img_source = output_denseCorr_multiframe_all(:,:,:,iterator);
    [tmp1_x, tmp1_y] = meshgrid(1:size(tmp_img_source,1), 1:size(tmp_img_source, 2));
    tmp1_x = reshape(tmp1_x, [numel(tmp1_x), 1]);
    tmp1_y = reshape(tmp1_y, size(tmp1_x));
    tmp2_x = reshape(tmp_img_source(:,:, 1), size(tmp1_x));
    tmp2_y = reshape(tmp_img_source(:,:, 2), size(tmp1_x));
    
    lines_to_skip = double(tmp2_x.^2 + tmp2_y.^2).^(1/2) < 8;
    
%     tmp1_x_size = size(tmp1_x)
%     tmp1_y_size = size(tmp1_y)
%     tmp2_x_size = size(tmp2_x)
%     tmp2_y_size = size(tmp2_y)
    
    % capture instance of the image to be drawn over
    tmp_img_target = output_digiLum_all(:,:,:,iterator);

    % debug
    imshow(tmp_img_target)
    
    for line_i = 1:numel(tmp1_x)
        if lines_to_skip(line_i)
            continue
        end
        
        fprintf(['l ' num2str(line_i) ' of ' num2str(numel(tmp1_x)) ' - ']);
        
        %// Draw lines from p1 to p2 on matrix tmp_ing_target
        % x = p1(1):p2(1)
        tmp_line_x = min(tmp1_x(line_i),tmp2_x(line_i)):max(tmp1_x(line_i),tmp2_x(line_i));
        tmp_line_x = reshape(tmp_line_x, [numel(tmp_line_x), 1]);

        % round((x - p1(1)) * (p2(2) - p1(2)) / (p2(1) - p1(1)) + p1(2));
        tmp_line_y = reshape(...
                                round(...
                                          (tmp_line_x - tmp1_x(line_i)) ...
                                        * (tmp2_y(line_i) - tmp1_y(line_i)) ...
                                        / (tmp2_x(line_i) - tmp1_x(line_i)) ...
                                        + tmp1_y(line_i) ...
                                    )...
                                , size(tmp_line_x)...
                      ); ...
                      % ;
        tmp_line_xxx = repmat(tmp_line_x, [3, 1]);
        tmp_line_yyy = repmat(tmp_line_y, [3, 1]);
        tmp_chan = reshape(...
                            repmat(...
                                [1,2,3], ...
                                [numel(tmp_line_x), 1]), ...
                            size(tmp_line_xxx) ...
                        );
                    
        % clamp tmp_line_xxx and tmp_line_yyy to indexes withing image bounds
        tmp_line_xxx = max(min(tmp_line_xxx,size(tmp_img_target, 1)),1);
        tmp_line_yyy = max(min(tmp_line_yyy,size(tmp_img_target, 2)),1);
        
        % re-cast as doubles (otherwise, sub2ind() will throw a data type
        % error... not sure why)
        tmp_line_xxx    = double(tmp_line_xxx);
        tmp_line_yyy    = double(tmp_line_yyy);
        tmp_chan        = double(tmp_chan);
        
        % draw white into the indexes
        % m(sub2ind(size(m), y, x, channel)) = 1;
        inds = sub2ind(size(tmp_img_target), tmp_line_xxx, tmp_line_yyy, tmp_chan);
        tmp_img_target(inds) = ui8_max;
    end
    
    % Draw colored lines along vectors of the field into the digiluminescence
    % effect output
    output_digiLum_all(:,:,:,iterator) = tmp_img_target;
    
    % debug
    imshow(tmp_img_target);
    
    % print time
    fprintf('\n')
    toc
end
% %test with quiver plot
% quiver(tmp1_x, tmp1_y, tmp2_x, tmp2_y, ...
%     1, ... scale
%     '-','filled', ... linespec, markerfilling
%     'Linewidth',1, 'MarkerSize',1); % other properties

% %% Assign plots to digilum effect
% % TODO: replace these hack values for digilum effect
% output_digiLum_all = uint8(output_denseCorr_multiframe_all + ui8_hlf);
% output_digiLum_all(:,:,3,:) = ui8_max;

% clean up
clear iteration*

% print time
toc

%% Save out test files
fprintf('----\n');
fprintf('Saving out some test files \n');

% reformat data
tic
fprintf([' - reformatting data - ']);
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
%     tmp_data_C_all                      = uint8(data_C_all                              )               ;
%     tmp_data_D_all                      = uint8(data_D_all              / i16_2_ui8     )               ;
    tmp_output_cleanPlate               = uint8(output_cleanPlate       / i16_2_ui8     )               ;
    tmp_output_uMasks_all               = uint8(output_uMasks_all       / i16_2_ui8     )               ;
    tmp_output_denseCorr_all            = uint8(output_denseCorr_all                    )   + dc_offset ;
    tmp_output_grid_all                 = uint8(output_grid_all                         )               ;
    tmp_output_denseCorr_masked_all     = uint8(output_denseCorr_masked_all             )   + dc_offset ;
    tmp_output_denseCorr_multiframe_all = uint8(output_denseCorr_multiframe_all         )   + dc_offset ;
    tmp_output_digiLum_all              = uint8(output_digiLum_all                      )               ;
    % must have a [w,h,bitDepth, frames] array for video file writing
%     tmp_data_D_all = permute(tmp_data_D_all, [1,2,4,3]);
    tmp_output_uMasks_all = permute(tmp_output_uMasks_all, [1,2,4,3]);
    % scale and translate dense correspondence for visibility
    tmp_output_denseCorr_all            = uint8((double(tmp_output_denseCorr_all            ) - dc_offset) * dc_scale + dc_offset);
    tmp_output_denseCorr_masked_all     = uint8((double(tmp_output_denseCorr_masked_all     ) - dc_offset) * dc_scale + dc_offset);
%     tmp_output_denseCorr_multiframe_all = uint8((double(tmp_output_denseCorr_multiframe_all ) - dc_offset) * dc_scale + dc_offset);
    %clean up
    clear dc_*
% print time
toc

%%%%%%%%%%%%
% images
%%%%%%%%%%%%
tic
fprintf([' - images - ']);
    % save out images
%     imwrite(tmp_data_C_all(:,:,:,1)                     ,[ 'test_01_Color.png'                  ]);
%     imwrite(tmp_data_D_all(:,:,:,1)                     ,[ 'test_02_Depth.png'                  ]);
    imwrite(tmp_output_cleanPlate                       ,[ 'test_02_Depth_cPlate.png'           ]);
    imwrite(tmp_output_uMasks_all(:,:,:,1)              ,[ 'test_03_uMask.png'                  ]);
    imwrite(tmp_output_denseCorr_all(:,:,:,1)           ,[ 'test_04_denseCorr.png'              ]);
    imwrite(tmp_output_grid_all(:,:,:,1)                ,[ 'test_05_grid_warped.png'            ]);
    imwrite(tmp_output_denseCorr_masked_all(:,:,:,1)    ,[ 'test_06_denseCorr_masked.png'       ]);
    imwrite(tmp_output_denseCorr_multiframe_all(:,:,:,1),[ 'test_06_denseCorr_multiframe.png'   ]);
    imwrite(tmp_output_digiLum_all(:,:,:,1)             ,[ 'test_06_digiLum.png'                ]);
% print time 
toc

%%%%%%%%%%%%
% videos
%%%%%%%%%%%%

% % data_C_all
% tic
% fprintf([' - videos - data_C_all - ']);
%     writerObj = VideoWriter(['test_01_Color.mp4'], 'MPEG-4');
%     open(writerObj);
%     writeVideo(writerObj,tmp_data_C_all)
%     close(writerObj);
% % print time
% toc
% 
% % data_D_all
% tic
% fprintf([' - videos - data_D_all - ']);
%     writerObj = VideoWriter(['test_02_Depth.mp4'], 'MPEG-4');
%     open(writerObj);
%     writeVideo(writerObj,tmp_data_D_all)
%     close(writerObj);
% % print time
% toc

% output_uMasks_all
tic
fprintf([' - videos - output_uMasks_all - ']);
    writerObj = VideoWriter(['test_03_uMask.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_uMasks_all)
    close(writerObj);
% print time
toc

% output_denseCorr_all
tic
fprintf([' - videos - output_denseCorr_all - ']);
    writerObj = VideoWriter(['test_04_denseCorr.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_all)
    close(writerObj);
% print time
toc

% output_grid_all
tic
fprintf([' - videos - output_grid_all - ']);
    writerObj = VideoWriter(['test_05_grid_warped.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_grid_all)
    close(writerObj);
% print time
toc

% output_denseCorr_masked_all
tic
fprintf([' - videos - output_denseCorr_masked_all - ']);
    writerObj = VideoWriter(['test_06_denseCorr_masked.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_masked_all)
    close(writerObj);
% print time
toc

% output_denseCorr_multiframe_all
tic
fprintf([' - videos - output_denseCorr_multiframe_all - ']);
    writerObj = VideoWriter(['test_06_denseCorr_multiframe.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_denseCorr_multiframe_all)
    close(writerObj);
% print time
toc

% output_digiLum_all
tic
fprintf([' - videos - output_digiLum_all - ']);
    writerObj = VideoWriter(['test_06_digiLum.mp4'], 'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,tmp_output_digiLum_all)
    close(writerObj);
% print time
toc

% clean up memory
clear writerObj tmp_*


%% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n'); 
%TODO: figure out a way to print all elapsed time for this function
fprintf('====\n');
end