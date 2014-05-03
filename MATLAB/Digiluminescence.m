function [cleanPlate, digiluminescence_all] = digiluminescence(C_all, D_all, joint_positions_all, timestamps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Start timer
tic
fprintf('====\n');
fprintf('Digiluminescence :: Executing\n');

%% Handle default arguments
for i = 1 % For loop is for code collapsing only (so I don't have to look at these)

    % set default value for first argument
    if( nargin < 1 )
        ;
    end
   
end

fprintf('====\n');
fprintf('handled default arguments \n');
toc


%% Initialize return values
cleanPlate = zeros(size(D_all(:,:,1)), 'int16');
digiluminescence_all = zeros(size(C_all),'int8');

% cleanPlateSize = size(cleanPlate)
% j_max = size(D_all, 3)

fprintf('----\n');
fprintf('initialized return values \n');
toc

%% Create a clean plate for the depth data
for j = 1:size(D_all, 3) % TODO: ask Il if there's a way to do this without a for loop
    cleanPlate = max(cleanPlate, D_all(:,:,j));
end

imshow(cleanPlate);
% TODO: save this out

fprintf('----\n');
fprintf('created clean plate for depth data \n');
toc

%% Create digiluminescence effect frame by frame 
% TODO: ask Il if there's a way to do this without a for loop (i.e.:
% calculating large matrixes for all frames and doing comparisons at once
% for entire piece

% Calculate projective joint positions of this frame

% Skip the following on the first frame

% Calculate vector field for current frame based on join positions of this
% frame and the previous frame

% Calculate movement mask based on clean slate difference matting with
% current frame

% Mask vector field with movement mask

% Fade effect from previous frame

% Draw colored lines along vectors of the field

fprintf('----\n');
fprintf('create digiluminescence effect frame by frame \n');
toc

% Report timestamp
fprintf('====\n');
fprintf('Digiluminescence :: End\n');
toc
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

