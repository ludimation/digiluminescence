clear all; close all;

kinect = RobotRaconteur.Connect('tcp://127.0.0.1:4444/KinectService/Kinect');

x0 = [1;0;0]; y0 = [0;1;0]; z0 = [0;0;1]; zed = [0;0;0];

%%

kinect.EnableAllStreams();

%% 

dt = 1/30;
T = 30;
t = 0:dt:T;

C_all = zeros(480,640,3,length(t),'uint8');
D_all = zeros(480,640,length(t),'int16');
joint_positions_all = zeros(20,3,length(t));
timestamps = zeros(1,length(t));

jp = zeros(20,3);

tic;
for k=1:length(t)
    timestamps(k) = toc;

    c = kinect.current_color_image;
    d = kinect.current_depth_image;

    C = derasterBGRaImage(c,[640 480]);
    D = reshape(d,[640 480])';

    if kinect.num_tracked_people > 0
        p = kinect.get_person(int32(0));
        jp = reshape(p.joint_positions,[3 20])';
    end
    
    C_all(:,:,:,k) = C;
    D_all(:,:,k) = D;
    joint_positions_all(:,:,k) = jp;

    
    t1 = toc;
    while t1 - timestamps(k) < dt
        t1 = toc;
    end
end

%% Test display

figure(1);
h_p = plot3(0,0,0,'.b');
axis equal;
axis([-2 2 -2 2 1 3]);
%%

for k=1:length(timestamps)
    set(h_p,'XData',joint_positions_all(:,1,k), ...
            'YData',joint_positions_all(:,2,k), ...
            'ZData',joint_positions_all(:,3,k));
    pause(0.01);
        drawnow;
end