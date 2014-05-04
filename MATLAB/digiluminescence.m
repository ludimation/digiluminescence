function [ out_dl_all, out_D_cPlate, out_uMasks_all, out_j_features, out_denseCorr_all, out_grid_all ] = ...
    digiluminescence(C_all, D_all, joint_positions_all, timestamps, calcDenseCorr)
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

    % set default value for first argument
    if( nargin < 5 )
        calcDenseCorr = false;
    end
   
end

% print time
toc

%% Initialize variables
tic
fprintf('----\n');
fprintf('Initializing variables \n');

out_D_cPlate       	= zeros(    size(D_all(:,:,1))          , 'int16'   );
out_uMasks_all      = zeros(    size(D_all)                 , 'int16'   );
out_denseCorr_all   = zeros(    size(C_all)                 , 'uint8'   );
out_grid_all        = zeros(    size(C_all)                 , 'uint8'   );
grid_template       = zeros(    size(out_grid_all(:,:,:,1)) , 'uint8'   );
out_dl_all          = zeros(    size(C_all)                 , 'uint8'   );

n_joints            = size(     joint_positions_all         , 1         );
n_frames            = length(   timestamps                              );

ui8_max = intmax('uint8');
ui8_hlf = round(intmax('uint8')/2);
i16_max = intmax('int16');
u16_2_ui8 = 2^7;

% Draw grids : drawGrid (I, spcGrid, spcPoints, color)
grid_template   = drawGrid(grid_template    , 32, 2, [ui8_max, ui8_max, 0]); % yellow
out_grid_all    = drawGrid(out_grid_all     , 32, 1, [ui8_max, ui8_max, 0]); % yellow

% print time
toc

%% Create a clean plate of the environment
tic
fprintf('----\n');
fprintf('Creating clean plate for depth data \n');

% out_D_cPlate is an image that represents the background of the scene,
% which is essentially the maximum depth values found in all frames of the
% depth data
out_D_cPlate = max(D_all,[],3);

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
inds_positive           = find(D_all > -8); 
D_all_cPlate            = repmat(out_D_cPlate, [1,1,size(D_all,3)]);
D_all_clean             = D_all_cPlate;
D_all_clean(inds_positive) = D_all(inds_positive);
% calculate masks for areas with large differences between the clean Depth
% image and the cleanPlate
D_all_diff              = abs(D_all_clean - D_all_cPlate);
inds_BG                 = find(abs(D_all_clean - D_all_cPlate) < 1024); %TODO: Make this threshold user-specifiable
out_uMasks_all          = D_all_clean;
out_uMasks_all(inds_BG) = -8;

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
j_pos_all_reshaped = reshape(permute(joint_positions_all, [2 1 3]), 3, n_frames*n_joints);
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
out_j_features = [ j_endingFeatures; j_startingFeatures];

% TODO: 
%    - add additional feature points along the lines of larger key limbs \
%    torso
%    - might also need a few additional features outside bounds of figure
%    positions to "tack" down the surrounding frame of the image so it
%    doesn't warp around so much

% draw joints in yellow: drawPoints(points, I, size, color)
out_grid_all = drawPoints(j_startingFeatures, out_grid_all, 16, [ui8_max, ui8_max, 0] ); % yellow
%TODO: draw limbs in white
% drawLimbs(j_pos_all_projective, out_grid_all, 4, [ui8_max, ui8_max, 0] ); % yellow

% cleanup
clear j_pos_all_reshaped j_pos_all_reshaped_projective
clear j_pos_all_x j_pos_all_y j_pos_all_z
clear j_pos_all_x_prjctd j_pos_all_y_prjctd j_pos_all_z_prjctd 

% print time
toc

%% Calculate dense correspondence fields frame by frame
% tic
fprintf('----\n');
fprintf('Calculating dense correspondence fields frame by frame \n');

% create dense correspondence fields one frame at a time
if calcDenseCorr
    for iterator = 1:n_frames
        tic
            fprintf([' - frame ' num2str(iterator) ' - ']);
            [ denseCorr, grid ] = thin_plate_denseCorrespondence(...
                out_j_features(:, :, iterator), ...
                out_grid_all(:,:,:, iterator) ...
            );
            out_denseCorr_all(:,:,:, iterator) = denseCorr;
            out_grid_all(:,:,:, iterator) = grid;
        % print time
        toc
    end
end

% draw grid again sparser and blue to show where it started : drawGrid (I, spcGrid, spcPoints, color)
out_grid_all = drawGrid(out_grid_all, 32, 4, [0 , ui8_max , ui8_max]); % cyan
% draw ending joints in magenta (small radius -- warped joints will appear as large yellow dots)
out_grid_all = drawPoints(j_endingFeatures      , out_grid_all, 8, [ui8_max , 0         , ui8_max] ); % magenta
% draw starting joints in blue (smaller radius)
out_grid_all = drawPoints(j_startingFeatures    , out_grid_all, 6, [0       , ui8_max   , ui8_max] ); % cyan

% cleanup
clear j_pos_all_projective j_pos_all_projective2
clear iterator denseCorr grid

% print time
% toc

%% Create digiluminescence effect from dense correspondence fields
tic
fprintf('----\n');
fprintf('Creating digiluminescence effect frame by frame \n');

% Mask dense correspondence field with movement mask

% Fade effect from previous frame

% Draw colored lines along vectors of the field

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
% TODO: include drawings of old, new, and warped positions of joints/limbs
% in grid images
imwrite( C_all(:,:,:,1)                                     ,[ 'test_01_Color.png'          ]);
imwrite(uint8( D_all(:,:,1)                 / u16_2_ui8 )   ,[ 'test_02_Depth.png'          ]);
imwrite(uint8( out_D_cPlate                 / u16_2_ui8 )   ,[ 'test_02_Depth_cPlate.png'   ]);
imwrite(uint8( out_uMasks_all(:,:,1)        / u16_2_ui8 )   ,[ 'test_03_uMask.png'          ]);
imwrite(uint8( out_denseCorr_all(:,:,:,1) )                 ,[ 'test_04_denseCorr.png'      ]);
imwrite(uint8( grid_template(:,:,:) )                       ,[ 'test_05_grid_template.png'  ]);
imwrite(uint8( out_grid_all(:,:,:,1) )                      ,[ 'test_05_grid_warped.png'    ]);
% print time
toc

%%%%%%%%%%%%
% videos
%%%%%%%%%%%%
tic
fprintf([' - videos - reformatting data - ']);
% must have a [w,h,bitDepth, frames] array for video file writing
D_all = permute(D_all, [1,2,4,3]);
out_uMasks_all = permute(out_uMasks_all, [1,2,4,3]);
% IMG must be of one of the following classes: double, single, uint8
C_all                   = uint8(C_all                               );
D_all                   = uint8(D_all               / u16_2_ui8     );
out_uMasks_all          = uint8(out_uMasks_all      / u16_2_ui8     );
out_denseCorr_all       = uint8(out_denseCorr_all                   );
out_grid_all            = uint8(out_grid_all                        );
% print time
toc

% C_all
tic
fprintf([' - videos - C_all - ']);
writerObj = VideoWriter(['test_01_Color.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,C_all)
close(writerObj);
% print time
toc

% D_all
tic
fprintf([' - videos - D_all - ']);
writerObj = VideoWriter(['test_02_Depth.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,D_all)
close(writerObj);
% print time
toc

% out_uMasks_all
tic
fprintf([' - videos - out_uMasks_all - ']);
writerObj = VideoWriter(['test_03_uMask.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,out_uMasks_all)
close(writerObj);
% print time
toc

% out_denseCorr_all
tic
fprintf([' - videos - out_denseCorr_all - ']);
writerObj = VideoWriter(['test_04_denseCorr.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,out_denseCorr_all)
close(writerObj);
% print time
toc

% out_grid_all
tic
fprintf([' - videos - out_grid_all - ']);
writerObj = VideoWriter(['test_05_grid_warped.mp4'], 'MPEG-4');
open(writerObj);
writeVideo(writerObj,out_grid_all)
close(writerObj);
% print time
toc

% clean up memory
clear C_all D_all joint_positions_all timestamps writerObj


%% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n'); 
%TODO: figure out a way to print all elapsed time for this function
fprintf('====\n');
end

%% DRAWGRID method
function [ I ] = drawGrid (I, spcGrid, spcP, color)
    if nargin < 2
        spcGrid = 10;
    end
    if nargin < 3
        spcP = 1;
    end
    if nargin < 4
        color = white(1) * ui8_max;
    end
   
    I(1                     , spcP:spcP:end       , 1,:) = color(1); % first row r
    I(1                     , spcP:spcP:end       , 2,:) = color(2); % first row g
    I(1                     , spcP:spcP:end       , 3,:) = color(3); % first row b
    I(spcP:spcP:end         , 1                   , 1,:) = color(1); % first column r
    I(spcP:spcP:end         , 1                   , 2,:) = color(2); % first column g
    I(spcP:spcP:end         , 1                   , 3,:) = color(3); % first column b
    I(spcGrid:spcGrid:end   , spcP:spcP:end       , 1,:) = color(1); % everything in between r
    I(spcGrid:spcGrid:end   , spcP:spcP:end       , 2,:) = color(2); % everything in between g
    I(spcGrid:spcGrid:end   , spcP:spcP:end       , 3,:) = color(3); % everything in between b
    I(spcP:spcP:end         , 10:spcGrid:end      , 1,:) = color(1); 
    I(spcP:spcP:end         , 10:spcGrid:end      , 2,:) = color(2);
    I(spcP:spcP:end         , 10:spcGrid:end      , 3,:) = color(3); % TODO: there HAS to be a more elegant way to do this
end

%% DRAWPOINTS method
function [ I ] = drawPoints(p_array, I, sz_draw, c_rgb)
% draws squares around point of the specified size and color into the image
% provided

    % check inputs
    narginMin = 2;
    narginMax = 4;
    narginchk(narginMin, narginMax);
    if nargin < 3
        p_size = 8;
    end
    if nargin < 4
        c_rgb = white(1) * intmax('int16');
    end
    
    sz_I = size(I);
    I_px = -ones(sz_I, 'int16');
    % indeces must be rounded doubles / integers %TODO: might also want to
    % check values to make sure none of them are 0, nor max width and
    % height of I
%     p_array
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

    % create one array for each channel that contains the indexes of the
    % pixels to be colored
    point_indexes_r = sub2ind(sz_I, x_array, y_array, r_array, frame_array);
    point_indexes_g = sub2ind(sz_I, x_array, y_array, g_array, frame_array);
    point_indexes_b = sub2ind(sz_I, x_array, y_array, b_array, frame_array);
    
    % color the pixels
    I_px(point_indexes_r) = c_rgb(1);
    I_px(point_indexes_g) = c_rgb(2);
    I_px(point_indexes_b) = c_rgb(3);

    % TODO: add points to points array to account for point size 
    % start drawing at the top left of the point
    sz_draw_half = round(sz_draw/2);
    I_px = circshift(I_px, [-sz_draw_half,-sz_draw_half]);
    for pos_x = 1:sz_draw
        for pos_y = 1:sz_draw
            I_tmp = circshift(I_px, [pos_x, pos_y]);
            inds = find(I_tmp>-1);
            I(inds) = I_tmp(inds);
%             I = I + circshift(I_px, [pos_x, pos_y]);
        end
    end
end