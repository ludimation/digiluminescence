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
out_denseCorr_all   = zeros(    size(C_all)                 , 'int8'    );
out_grid_all        = zeros(    size(C_all)                 , 'int8'    );
out_dl_all          = zeros(    size(C_all)                 , 'int8'    );

n_joints            = size(     joint_positions_all         , 1         );
n_frames            = length(   timestamps                              );

% drawGrid (I, spcGrid, spcPoints, color)
grid_template = zeros(size(out_grid_all(:,:,:,1)), 'int8');
grid_template = drawGrid(grid_template, 10, 2, [2^8, 2^8, 0]);
out_grid_all = drawGrid(out_grid_all, 10, 1, [2^8, 2^8, 0]);

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

% Create second array for easy grabbing of one set of joints from the
% previous frame (frame 1 is compared against frame n which doesn't quite
% make sense, but I'd just like to see what happens)
j_pos_all_projective2 = circshift(j_pos_all_projective, 1);

% Create list of feature correspondences in a shape that thin plate dense
% corespondence function needs
%    TODO:
%    - create a 3D thin plate function and update this to include z
out_j_features = [j_pos_all_projective(:, 1:2, :), j_pos_all_projective2(:, 1:2, :)];
out_j_features = permute(out_j_features, [2,1,3]);

% create dense correspondence fields one frame at a time
if calcDenseCorr
    for iterator = 1:n_frames
        tic
            fprintf([' - frame ' num2str(iterator) ' - ']);
            [ denseCorr, grid ] = thin_plate_denseCorrespondence(out_j_features(:, :, iterator), out_grid_all(:,:,:, iterator) );
            out_denseCorr_all(:,:,:, iterator) = denseCorr;
            out_grid_all(:,:,:, iterator) = grid;
        % print time
        toc
    end
end

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
imwrite( C_all(:,:,:,1)                             ,[ 'test_01_Color.png'          ]);
imwrite(uint8( D_all(:,:,1)                 / 2^8 ) ,[ 'test_02_Depth.png'          ]);
imwrite(uint8( out_D_cPlate                 / 2^8 ) ,[ 'test_02_Depth_cPlate.png'   ]);
imwrite(uint8( out_uMasks_all(:,:,1)        / 2^8 ) ,[ 'test_03_uMask.png'          ]);
imwrite(uint8( out_denseCorr_all(:,:,:,1) )         ,[ 'test_04_denseCorr.png'      ]);
imwrite(uint8( grid_template(:,:,:) )               ,[ 'test_05_grid_template.png'  ]);
imwrite(uint8( out_grid_all(:,:,:,1) )              ,[ 'test_05_grid_warped.png'    ]);
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
C_all                   = uint8(C_all                       );
D_all                   = uint8(D_all               / 2^8   );
out_uMasks_all          = uint8(out_uMasks_all      / 2^8   );
out_denseCorr_all       = uint8(out_denseCorr_all           );
out_grid_all            = uint8(out_grid_all                );
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

%% Local methods
function [ I ] = drawGrid (I, spcGrid, spcPoints, color)
    if nargin < 2
        spcGrid = 10;
    end
    if nargin < 3
        spcPoints = 1;
    end
    if nargin < 4
        color = white(1) * 2^8;
    end
    
    I(1                     , 1:spcPoints:end       , 1,:) = color(1); % first row r
    I(1                     , 1:spcPoints:end       , 2,:) = color(2); % first row g
    I(1                     , 1:spcPoints:end       , 3,:) = color(3); % first row b
    I(1:spcPoints:end       , 1                     , 1,:) = color(1); % first column r
    I(1:spcPoints:end       , 1                     , 2,:) = color(2); % first column g
    I(1:spcPoints:end       , 1                     , 3,:) = color(3); % first column b
    I(spcGrid:spcGrid:end   , 1:spcPoints:end       , 1,:) = color(1); % everything in between r
    I(spcGrid:spcGrid:end   , 1:spcPoints:end       , 2,:) = color(2); % everything in between g
    I(spcGrid:spcGrid:end   , 1:spcPoints:end       , 3,:) = color(3); % everything in between b
    I(1:spcPoints:end       , 10:spcGrid:end        , 1,:) = color(1); 
    I(1:spcPoints:end       , 10:spcGrid:end        , 2,:) = color(2);
    I(1:spcPoints:end       , 10:spcGrid:end        , 3,:) = color(3); % TODO: there HAS to be a more elegant way to do this
end


%% Load images
% % TODO: Generalize this into a GUI that allows you to select images
% 
% % Check dimensions of all images to make sure they match, as well as bit
% % depth of S0 & T0 to make sure they match
% if size(S0,1) ~= w | size(S0,2) ~= h | size(S0,3) ~= chan | size(m,1) ~= w | size(m,2) ~= h 
%     % Throw error and return from function
%     fprintf('====\n'); 
%     fprintf('Error - Poisson_dallen2ndAttempt:\n'); 
%     fprintf('--\n'); 
%     fprintf('source, target, and mask image dimensions must match,\n'); 
%     fprintf('and source, target must have same number of channels\n'); 
%     fprintf(' - Target is %d x %d x %d\n', size(T0,1), size(T0,2), size(T0,3)); 
%     fprintf(' - Source is %d x %d x %d\n', size(S0,1), size(S0,2), size(S0,3)); 
%     fprintf(' - Mask is %d x %d x %d\n', size(m,1), size(m,2), size(m,3));
%     fprintf('====\n'); 
%     return;
% end
% 
% %% Save
% % Debugging images
% if false % These never execute
%     imwrite(uint8( M0*255   ),[ iFileName '-M0.png'         ]);
%     imwrite(uint8( Lt       ),[ iFileName '-Lt.png'         ]);
%     imwrite(uint8( Ls       ),[ iFileName '-Ls.png'         ]);
%     imwrite(uint8( LtAmp    ),[ iFileName '-LtAmp.png'      ]);
%     imwrite(uint8( LsAmp    ),[ iFileName '-LsAmp.png'      ]);
%     imwrite(uint8( LgVF*255 ),[ iFileName '-LgVF.png'       ]);
%     imwrite(uint8( Lgm*255  ),[ iFileName '-Lgm.png'        ]);
% end
% 
% % Move into debugging if statement as they begin working as expected
% imwrite(uint8( L0       ),[ iFileName '-L0.png'         ]);
% 
% % Save final image
% imwrite(uint8( I0        ),[ iFileName '.png'        ]);
% 
% 
% %% Cleanup Variables (some of these could be cleared earlier to save memory)
% clear iFileName nonC iterations dimReturnsThreshold verbose
% clear w h chan 
% clear S0 M0 T0 K Lt Ls LtAmp LsAmp LgVF Lgm Lm L0 I0
% 
% end

