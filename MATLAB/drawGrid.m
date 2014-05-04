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
    I(spcGrid:spcGrid:end   , spcP:spcP:end       , 3,:) = color(3); % everything in between b, etc.
    I(spcP:spcP:end         , spcGrid:spcGrid:end , 1,:) = color(1); 
    I(spcP:spcP:end         , spcGrid:spcGrid:end , 2,:) = color(2);
    I(spcP:spcP:end         , spcGrid:spcGrid:end , 3,:) = color(3); 
    
    % TODO: there HAS to be a more elegant way to do this
end