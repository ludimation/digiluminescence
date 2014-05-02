%% 
load('Images/20140429_data_fromDanKruse/david_kinect_data2.mat');

%%
[cleanPlate, digiluminescence_all] = digiluminescence(C_all, D_all, joint_positions_all, timestamps);

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