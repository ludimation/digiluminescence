function [ output_C_all, output_cleanPlate, output_uMasks_all, output_j_features, output_denseCorr_all, output_digiLum_all, output_grid_all, output_masked_denseCorr_all ] = ...
    digiluminescence(data_C_all, data_D_all, data_joint_positions_all, data_timestamps, data_mask_thresh, data_calcDenseCorr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% TODO: 
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

ui8_max = intmax('uint8');
ui8_hlf = round(intmax('uint8')/2);
i16_max = intmax('int16');
i16_2_ui8 = 2^7;

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
output_masked_denseCorr_all     =           output_denseCorr_all                            ;
output_digiLum_all              = ones(     size(data_C_all)                    , 'uint8'   ) * double(ui8_hlf);
output_grid_all                 = zeros(    size(data_C_all)                    , 'uint8'   );

% print time
toc

%% Draw grids
tic
fprintf('----\n');
fprintf('Drawing grids \n');

% drawGrid (I, spcGrid, spcPoints, color)
output_grid_all     = drawGrid(output_grid_all     , 32, 1, [ui8_max, ui8_max, 0]); % yellow

% print time
toc

%% Create a clean plate of the environment
tic
fprintf('----\n');
fprintf('Creating clean plate for depth data \n');

% out_D_cPlate is an image that represents the background of the scene,
% which is essentially the maximum depth values found in all frames of the
% depth data
output_cleanPlate = max(data_D_all,[],3);

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
            fprintf([' - frame ' num2str(iterator) ' - ']);
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

%% Draw debug joints and grid over warped grid
tic
fprintf('----\n');
fprintf('Drawing debug joints and grid \n');

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
output_masked_denseCorr_all = uint8( ...
            double(output_denseCorr_all) ...
        .*  double(repmat(permute(output_uMasks_all, [1,2,4,3]), [1,1,3,1])) ...
        /   double(i16_max) ...
    );

% create digiLum field from faded, circshifted versions of masked denseCorr
iteration_max = 16;
for iteration = 1:iteration_max
    % circshift iteration-1 so it starts with current frame
    output_digiLum_all = output_digiLum_all ...
        + circshift(output_masked_denseCorr_all, [1,1,1,iteration - 1] ) ...
        * (-((iteration-1)/iteration_max)+1)^(1/2) ... y = sqrt(-x+1) easing
        ;
end

% Draw colored lines along vectors of the field into the digiluminescence
% effect output

% print time
toc

%% Save out test files
fprintf('----\n');
fprintf('Saving out some test files \n');

%%%%%%%%%%%%
% images
%%%%%%%%%%%%
tic
fprintf([' - images - ']);
% scale dense correspondence
if data_calcDenseCorr
    % turn up coefficient to increase visualcontrast in dense
    % correspondence data
    dc_scale = double(2^3); 
else
    % if dense correspondence data is not being generated, it is already
    % maximal by default, so you don't need to scale it
    dc_scale = double(2^0);
end
dc_offset = double(0); % ui8_hlf;
tmp_output_denseCorr_all = (double(output_denseCorr_all) - dc_offset) * dc_scale + dc_offset;
tmp_output_masked_denseCorr_all = (double(output_masked_denseCorr_all) - dc_offset) * dc_scale + dc_offset;
% in grid images
imwrite( data_C_all(:,:,:,1)                                            ,[ 'test_01_Color.png'              ]);
imwrite(uint8( data_D_all(:,:,1)                / i16_2_ui8 )           ,[ 'test_02_Depth.png'              ]);
imwrite(uint8( output_cleanPlate                / i16_2_ui8 )           ,[ 'test_02_Depth_cPlate.png'       ]);
imwrite(uint8( output_uMasks_all(:,:,1)         / i16_2_ui8 )           ,[ 'test_03_uMask.png'              ]);
imwrite(uint8( tmp_output_denseCorr_all(:,:,:,1) )  + ui8_hlf           ,[ 'test_04_denseCorr.png'          ]);
imwrite(uint8( output_grid_all(:,:,:,1) )                               ,[ 'test_05_grid_warped.png'        ]);
imwrite(uint8( tmp_output_masked_denseCorr_all(:,:,:,1) )  + ui8_hlf    ,[ 'test_06_denseCorr_masked.png'   ]);
imwrite(uint8( output_digiLum_all(:,:,:,1) )                            ,[ 'test_06_digiLum.png'            ]);
% print time 
toc

%%%%%%%%%%%%
% videos
%%%%%%%%%%%%
tic
fprintf([' - videos - reformatting data - ']);
% must have a [w,h,bitDepth, frames] array for video file writing
tmp_data_D_all = permute(data_D_all, [1,2,4,3]);
tmp_output_uMasks_all = permute(output_uMasks_all, [1,2,4,3]);
% IMG must be of one of the following classes: double, single, uint8
tmp_data_C_all                  = uint8(data_C_all                              )               ;
tmp_data_D_all                  = uint8(tmp_data_D_all          / i16_2_ui8     )               ;
tmp_output_uMasks_all           = uint8(tmp_output_uMasks_all   / i16_2_ui8     )               ;
tmp_output_denseCorr_all        = uint8(tmp_output_denseCorr_all                )   + ui8_hlf   ;
tmp_output_grid_all             = uint8(output_grid_all                         )               ;
tmp_output_masked_denseCorr_all = uint8(tmp_output_masked_denseCorr_all         )   + ui8_hlf   ;
tmp_output_digiLum_all          = uint8(output_digiLum_all                      )               ;
% print time
toc

% data_C_all
tic
fprintf([' - videos - data_C_all - ']);
writerObj = VideoWriter(['test_01_Color.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,tmp_data_C_all)
close(writerObj);
% print time
toc

% data_D_all
tic
fprintf([' - videos - data_D_all - ']);
writerObj = VideoWriter(['test_02_Depth.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,tmp_data_D_all)
close(writerObj);
% print time
toc

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

% output_masked_denseCorr_all
tic
fprintf([' - videos - output_masked_denseCorr_all - ']);
writerObj = VideoWriter(['test_06_denseCorr_masked.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,tmp_output_masked_denseCorr_all)
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
clear C_all D_all joint_positions_all timestamps writerObj
clear writerObj tmp_*


%% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n'); 
%TODO: figure out a way to print all elapsed time for this function
fprintf('====\n');
end