%% Clear previous output data
clear n_framesToProcess ans output_*

%% Load data
tic
fprintf('====\n');
fprintf('Demo :: Loading kinect data \n');

% Only load data if the data has not yet been loaded
if (...
        ~exist('data_C_all'                 , 'var'	) || ...
        ~exist('data_D_all'                 , 'var'	) || ...
        ~exist('data_joint_positions_all'   , 'var' ) || ...
        ~exist('data_timestamps'            , 'var' )...
        )
    load('Images/20140429_data_fromDanKruse/david_kinect_data2.mat');
end

% print time
toc
fprintf('====\n');

%% Run digiluminescence
data_mask_thresh = 256
data_calcDenseCorr = false;
data_calcDenseCorr = true; %comment this line out to skip dense correspondence caculcations
data_n_framesToProcess = 31:32; % 2 frames = smallest test size
% data_n_framesToProcess = 31:35; % 5 frames
% data_n_framesToProcess = 31:80; % 50 frames
% data_n_framesToProcess = 1:300;
% data_n_framesToProcess = 1:length(data_timestamps); % Uncomment this line for full data processing, comment it out for smaller set above

% TODO: 
%    - feed dense correspondences (if any) into this equation
%    - break down effect into more manageable function call chunks
[ output_C_all, output_cleanPlate, output_uMasks_all, output_j_features, output_denseCorr_all, output_digiLum_all, output_grid_all ] = digiluminescence( ...
        data_C_all(:,:,:, data_n_framesToProcess), ...
        data_D_all(:,:, data_n_framesToProcess), ...
        data_joint_positions_all(:,:, data_n_framesToProcess), ...
        data_timestamps(data_n_framesToProcess), ...
        data_mask_thresh, ...
        data_calcDenseCorr ... 
    );

% TODO: time-stamp and date this file
save digiluminescence_outputdata.mat output_*

%% Test display
% % TODO: add imshow to this figure
% 
% figure(1);
% h_p = plot3(0,0,0,'.b');
% axis equal;
% axis([-2 2 -2 2 1 3]);
% 
%%
% % TODO: update image being shown as well
% for k=1:length(timestamps)
%     set(h_p,'XData',joint_positions_all(:,1,k), ...
%             'YData',joint_positions_all(:,2,k), ...
%             'ZData',joint_positions_all(:,3,k));
%     pause(0.01);
%         drawnow;
% end