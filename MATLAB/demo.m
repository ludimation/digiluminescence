%% Clear previous output data
clear n_framesToProcess ans output_*

%% Load data
tic
fprintf('====\n');
fprintf('Demo :: Loading kinect data \n');

% Only load data if the data has not yet been loaded
if (...
        ~exist('data_C_all'                 , 'var'	) | ...
        ~exist('data_D_all'                 , 'var'	) | ...
        ~exist('data_joint_positions_all'   , 'var' ) | ...
        ~exist('data_timestamps'            , 'var' )...
        )
    load('Images/20140429_data_fromDanKruse/david_kinect_data2.mat');
end

% print time
toc
fprintf('====\n');

%% Run digiluminescence
n_framesToProcess = 30:35; % smaller test amount
% n_framesToProcess = 1:300;
% n_framesToProcess = 1:length(data_timestamps); % Uncomment this line for full data processing, comment it out for smaller set above
[ output_digiluminescence_all, output_cleanPlate, output_uMasks_all, output_j_features, output_denseCorr_all, output_grid_all ] = ...
    digiluminescence(...
        data_C_all(:,:,:, n_framesToProcess), ...
        data_D_all(:,:, n_framesToProcess), ...
        data_joint_positions_all(:,:, n_framesToProcess), ...
        data_timestamps(n_framesToProcess)...
        , true(1) ... comment this line out with '...' preceding ', true...' to skip calculate dense correspondence
    );

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