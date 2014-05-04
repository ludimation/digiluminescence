function [ out_denseCorrespondence, out_I_interp, out_w, out_a, out_b ] = thin_plate_denseCorrespondence(in_feature_pairs, in_I)
%THIN_PLATE_DENSECORRESPONDENCE 
% Based on script provided by Il-Young Son as example of 
%     Scattered Data Interpolation
%
% [denseCorrespondence, I_interp, w, a, b ] = thin_plate_denseCorrespondence(feature_pairs, I)
%
% INPUTS:
%
%   feats_pairs : feature correspondences. It should be a 4 by n matrix.  
%   First two rows are (x,y) coordinates of feature points in image 1.  
%   Last two rows are corresponding (x,y) coordinates of feature points in
%   image 2. 
%
%   I is image 1. 
%
% OUTPUTS:
%
%   denseCorrespondence : a matrix of the vectors from which each pixel in
%   image I_interp was calculated from image I
%
%   I_interp : image 1 interpolated based on the vector field created by
%   the feature correspondence pairs.
%

%% Handle default arguments
for i = 1 % For loop is for code collapsing only (so I don't have to look at these)

    no_out_I_interp = false;
    
    % set default value of no_out_I_interp based on whether or not there is
    % an out argument for the interpoloated image
    if( nargout < 2 )
        no_out_I_interp = true;
    end
   
end


%% Setup key variables

% store some sizes for future reference
n_feats = size(in_feature_pairs,2);
im_sz = size(in_I);
im_sz_wh = im_sz(1:2);
% im_sz = [size(in_I,1), size(in_I,2)]; % same as previous two lines?

% reformat I -- why is this necessary?
in_I = double(in_I);

% separate features for image1 and image2
feats_I1 = in_feature_pairs(1:2,:); % [x1, x2, ..., xn; y1, y2, ..., yn] for first set of features
feats_I2 = in_feature_pairs(3:4,:); % [x1, x2, ..., xn; y1, y2, ..., yn] for second set of features

%% set up the linear system

% creat matrix of r_ij values % TODO: still don't fully understand how this is constructed
x_i = repmat(feats_I1,n_feats,1);
x_j = repmat(reshape(feats_I1, numel(feats_I1), 1),1,n_feats);
z_ij = reshape(x_i-x_j,2,numel(x_i-x_j)/2);
r_ij = reshape(sqrt(sum(z_ij.*z_ij)),n_feats,n_feats);

% calculate phi matrix for A from r_ij vaulues
phi_lsys = (r_ij.^2).*(log(r_ij+...
    eye(n_feats)));
% since log of 1 is zero we can add 1s to the diagonal entries of r_ij
% matrix to prevent taking log of zero

% build A and b, then calcuate x
A = [phi_lsys feats_I1' ones(n_feats,1);feats_I1 zeros(2,3);ones(1,n_feats) zeros(1,3)];
b = [feats_I2';zeros(3,2)];
x = A\b;

% store out_w, out_a, and out_b portions of x
out_w = x(1:n_feats,:); % [ feats_n x 2 ] matrix for [w11, w21; w12, w22; ... w1n, w2n]
out_a = x(n_feats+1:n_feats+2,:); % [ 2 x 2 ] matrix for [a11, a21; a12, a22]
out_b = x(end,:); % [ 1 x 2 ] matrix for [b1, b2]

%% construct dense correspondence matrix

% create a [pixCount x 2] matrix of all the x,y pixel indixes possible in
% image I (pix_inputs) to feed into the correspondence formula
% store length of this list for future use
n_pix_inputs = prod(im_sz_wh);
xx = 1:im_sz_wh(1); % list of all x indixes [1,2,3,..., im_sz(1)]
yy = 1:im_sz_wh(2); % list of all y indixes [1,2,3,..., im_sz(2)]
xx = reshape(repmat(xx',1,im_sz_wh(2)),n_pix_inputs,1); % complete list of repeating cycles of x indeces [1;2;3;...;n;1;2;3;...;n;...]
yy = reshape(repmat(yy,im_sz_wh(1),1),n_pix_inputs,1); % complete list of stacked increasing y values [1;1;1;...;2;2;2;...;n;n;n;...]
pix_inputs = [xx yy]'; % [1,1; 2,1; 3,1; ...; n,1; 1,2; 2,2; 3,2; ...; n,2; ...]

% construct r_1 matrix to feed into phi_dc % TODO: still don't fully understand how this is constructed
x_ins = repmat(pix_inputs,n_feats,1);
x_is = repmat(reshape(feats_I1, numel(feats_I1), 1),1,n_pix_inputs);
zz = reshape(x_ins-x_is,2,numel(x_ins-x_is)/2);
r_i = reshape(sqrt(sum(zz.*zz)),n_feats,n_pix_inputs);
% since log of 1 is zero we can add 1s to the zero entries of r_i matrix to
% prevent taking log of zero
r_i(r_i==0) = 1; 

% compute phi_dc for dense correspondence formula
phi_dc = (r_i.^2).*(log(r_i));

% create dense correspondence field
x_interps = out_w'*phi_dc + out_a'*pix_inputs + repmat(out_b',1,n_pix_inputs);

% reshape dense correspondence field into x, y, and (TODO:) z vectors
out_denseCorrespondence = zeros(im_sz,'int16');
i16_max = intmax('int16');
i16_hlf = int16(i16_max / 2);
out_denseCorrespondence(:,:,1) = int16(reshape((x_interps(1,:) - xx'), im_sz_wh(1),im_sz_wh(2)));
out_denseCorrespondence(:,:,2) = int16(reshape((x_interps(2,:) - yy'), im_sz_wh(1),im_sz_wh(2)));
out_denseCorrespondence(:,:,3) = 0; % TODO: should eventually include z in dense correspondence calculations

if ~no_out_I_interp
    vr = reshape(in_I(:,:,1),n_pix_inputs,1);
    vg = reshape(in_I(:,:,2),n_pix_inputs,1);
    vb = reshape(in_I(:,:,3),n_pix_inputs,1);

    Fr = scatteredInterpolant(xx, yy, vr);
    Fg = scatteredInterpolant(xx, yy, vg);
    Fb = scatteredInterpolant(xx, yy, vb);

    out_I_interp(:,:,1) = reshape(Fr(x_interps(1,:)',x_interps(2,:)'),im_sz_wh(1),im_sz_wh(2));
    out_I_interp(:,:,2) = reshape(Fg(x_interps(1,:)',x_interps(2,:)'),im_sz_wh(1),im_sz_wh(2));
    out_I_interp(:,:,3) = reshape(Fb(x_interps(1,:)',x_interps(2,:)'),im_sz_wh(1),im_sz_wh(2));
end

end

