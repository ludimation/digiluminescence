function [ out_dl_all, out_D_cPlate, out_uMasks_all, out_denseCorr_all ] = digiluminescence(C_all, D_all, joint_positions_all, timestamps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Start timer
fprintf('====\n');
fprintf('Digiluminescence :: Executing\n');

%% Handle default arguments
tic
fprintf('====\n');
fprintf('Handling default arguments \n');
for i = 1 % For loop is for code collapsing only (so I don't have to look at these)

    % set default value for first argument
    if( nargin < 1 )
        ;
    end
   
end

% print time
toc

%% Initialize variables
tic
fprintf('----\n');
fprintf('Initializing variables \n');

out_D_cPlate                = zeros(size(D_all(:,:,1))  , 'int16'   );
out_uMasks_all              = zeros(size(D_all)         , 'int16'   );
out_denseCorr_all           = zeros(size(C_all)         , 'int8'    );
out_dl_all                  = zeros(size(C_all)         , 'int8'    );

n_joints            = size(joint_positions_all  , 1         );
n_frames            = length(timestamps                     );

% print time
toc

%% Create a clean plate for the depth data
tic
fprintf('----\n');
fprintf('Creating clean plate for depth data \n');

% out_D_cPlate is an image that represents the background of the scene,
% which is essentially the maximum depth values found in all frames of the
% depth data

% for j = 1:size(D_all, 3) 
%     out_D_cPlate = max(out_D_cPlate, D_all(:,:,j));
% end
out_D_cPlate = max(D_all,[],3); % does same as above, without a for loop

% print time
toc

%% Create user masks
tic
fprintf('----\n');
fprintf('Creating user masks \n');

% clean up depth images by putting cPlate in areas that have no value
inds_positive       = find(D_all > -8); 
D_all_cPlate        = repmat(out_D_cPlate, [1,1,size(D_all,3)]);
D_all_clean         = D_all_cPlate;
D_all_clean(inds_positive) = D_all(inds_positive);
% calculate the difference between the clean Depth image and the cleanPlate
% (BG)
D_all_diff          = abs(D_all_clean - D_all_cPlate);
inds_BG             = find(abs(D_all_clean - D_all_cPlate) < 1024); %TODO: Make this threshold user-specifiable
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

% Re-shape joint_positions_all so it can be fed into projectiveTransformation (p, to)
j_pos_all_reshaped = reshape(permute(joint_positions_all, [2 1 3]), 3, n_frames*n_joints);
j_pos_all_x = j_pos_all_reshaped(1,:);
j_pos_all_y = j_pos_all_reshaped(2,:);
j_pos_all_z = j_pos_all_reshaped(3,:);
% feed individual x, y, and z lists into projectiveTransformation (x, y, z, to)
[ j_pos_all_x_prjctd, j_pos_all_y_prjctd, j_pos_all_z_prjctd ] = ...
    projectiveTransformation(j_pos_all_x, j_pos_all_y, j_pos_all_z, 'projective');
j_pos_all_reshaped_projective = [ j_pos_all_x_prjctd; j_pos_all_y_prjctd; j_pos_all_z_prjctd ];
j_pos_all_projective = permute(reshape(j_pos_all_reshaped_projective, 3,n_joints,n_frames), [2 1 3]);

% cleanup
clear j_pos_all_reshaped j_pos_all_reshaped_projective
clear j_pos_all_x j_pos_all_y j_pos_all_z
clear j_pos_all_x_prjctd j_pos_all_y_prjctd j_pos_all_z_prjctd 

% print time
toc


%% Calculate dense correspondence fields frame by frame
tic
fprintf('----\n');
fprintf('Calculating densedense correspondence fields frame by frame \n');

% create second array for easy grabbing of one set of joints and the
% previous frame (frame 1 is compared against frame n just to see what
% happens)
j_pos_all_projective2 = circshift(j_pos_all_projective, 1);

% create dense correspondence fields
for iterator = 1:size(j_pos_all_projective, 3)
    j_pos_all_projective(:,1:2,iterator)';
    j_pos_all_projective2(:,1:2,iterator)';
end

% print time
toc

%% Create digiluminescence effect from dense correspondence fields
tic
fprintf('----\n');
fprintf('Creating digiluminescence effect frame by frame \n');

% Mask dense correspondence field with movement mask

% Fade effect from previous frame

% Draw colored lines along vectors of the field

% print time
toc

%% Save out some test images
tic
fprintf('----\n');
fprintf('Saving out some test images \n');

imwrite( C_all(:,:,:,1)                     ,[ 'test_01_Color.png'         ]);
imwrite(uint8( D_all(:,:,1) / 256 )         ,[ 'test_02_Depth.png'         ]);
imwrite(uint8( out_D_cPlate / 256 )         ,[ 'test_03_Depth_cPlate.png'  ]);
imwrite(uint8( out_uMasks_all(:,:,1) / 256 ),[ 'test_04_uMask.png'         ]);
% dense correspondence field
% effect

% clean up memory
clear C_all D_all joint_positions_all timestamps

% print time
toc


%% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n'); 
%TODO: figure out a way to print all elapsed time for this function
fprintf('====\n');



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

