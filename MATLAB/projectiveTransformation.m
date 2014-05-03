function [ x_trns, y_trns, z_trns ] = projectiveTransformation (x, y, z, to)
%Assumes a 640 x 480 projection, and kinect-specific depth -- could be
%     updated to take more args

%% Initialize calculation variables
WIDTH       = 640           ;
HEIGHT      = 480           ;
MAXDEPTH    = 10000         ;
HALFWIDTH   = WIDTH / 2     ;
HALFHEIGHT  = HEIGHT / 2    ;
HFOV        = 1.01447       ; %fXtoZ
VFOV        = 0.789809      ; % fYtoZ
COEFX       = WIDTH / HFOV  ;
COEFY       = HEIGHT / VFOV ;
HFOVTANHALF = tan(HFOV/2)   ;
VFOVTANHALF = tan(VFOV/2)   ;

% x_trns = zeros(size(x),'double'); 
% y_trns = zeros(size(y),'double'); 
% z_trns = zeros(size(z),'double'); 

if to == 'projective'
    %     ofPoint projective;
    % 	  projective.x = COEFX * p.x / p.z + HALFWIDTH;
    %     projective.y = HALFHEIGHT - COEFY * p.y / p.z;
    %     projective.z = p.z;
    % 	return projective;
    % }

    x_trns = COEFX * x ./ z + HALFWIDTH;
    y_trns = HALFHEIGHT - COEFY * y ./ z;
    z_trns = z;
else
    if to == 'world'
        %     ofPoint world;
        %     world.x = (p.x / WIDTH - 0.5) * p.z * HFOV;
        %     world.y = (0.5 - p.y / HEIGHT) * p.z * VFOV;
        %     world.z = p.z;
        % 	return world;
        % }
        x_trns = (x / WIDTH - 0.5) .* z * HFOV;
        y_trns = (0.5 - y / HEIGHT) .* z * VFOV;
        z_trns = z;
    else
        % 'to' argument is not one of the valid types, throw an error
    end
end

end
