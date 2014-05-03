function [ w, a, b, I_interp ] = thin_plate_interp2(x_feats,I)
%THIN_PLATE_INTERP Thin-plate spline interpolation
%Provided by Il-Young Son as example of 
%     Scattered Data Interpolation
%   [w, a, b, x_interps] = thin_plate_interp2(x_feats, x_inputs)
%   x_feats : feature correspondences. It should be a 4 by n matrix.  First
%   two rows are (x,y) coordinates of feature points in image 1.  Last two
%   rows are corresponding (x,y) coordinates of feature points in image 2.
%


n_feats = size(x_feats,2);
im_sz = size(I);
im_sz = im_sz(1:2);

I = double(I);

x_feats1 = x_feats(1:2,:);
x_feats2 = x_feats(3:4,:);
% creat matrix of r_ij values
x_i = repmat(x_feats1,n_feats,1);
x_j = repmat(reshape(x_feats1, numel(x_feats1), 1),1,n_feats);
z_ij = reshape(x_i-x_j,2,numel(x_i-x_j)/2);
r_ij = reshape(sqrt(sum(z_ij.*z_ij)),n_feats,n_feats);

% since log of 1 is zero we can add 1s to the diagonal entries of r_ij
% matrix to prevent taking log of zero
phi = (r_ij.^2).*(log(r_ij+eye(n_feats)));

% set up the linear system
A = [phi x_feats1' ones(n_feats,1);x_feats1 zeros(2,3);ones(1,n_feats) zeros(1,3)];
y = [x_feats2';zeros(3,2)];

z = A\y;

w = z(1:n_feats,:);
a = z(n_feats+1:n_feats+2,:);
b = z(end,:);

% compute r_i's to feed into phi_in for interpolation
xx = 1:im_sz(1);
yy = 1:im_sz(2);

xx = reshape(repmat(xx',1,im_sz(2)),prod(im_sz),1);
yy = reshape(repmat(yy,im_sz(1),1),prod(im_sz),1);
x_inputs = [xx yy]';
n_inputs = size(x_inputs,2);

x_ins = repmat(x_inputs,n_feats,1);
x_is = repmat(reshape(x_feats1, numel(x_feats1), 1),1,n_inputs);
zz = reshape(x_ins-x_is,2,numel(x_ins-x_is)/2);
r_i = reshape(sqrt(sum(zz.*zz)),n_feats,n_inputs);
r_i(find(r_i==0)) = 1;
% compute phi_in for interpolation
phi_in = (r_i.^2).*(log(r_i));


% interpolate points
x_interps = w'*phi_in + a'*x_inputs + repmat(b',1,n_inputs);

vr = reshape(I(:,:,1),prod(im_sz),1);
vg = reshape(I(:,:,2),prod(im_sz),1);
vb = reshape(I(:,:,3),prod(im_sz),1);

Fr = scatteredInterpolant(xx,yy,vr);
Fg = scatteredInterpolant(xx,yy,vg);
Fb = scatteredInterpolant(xx,yy,vb);

I_interp(:,:,1) = reshape(Fr(x_interps(1,:)',x_interps(2,:)'),im_sz(1),im_sz(2));
I_interp(:,:,2) = reshape(Fg(x_interps(1,:)',x_interps(2,:)'),im_sz(1),im_sz(2));
I_interp(:,:,3) = reshape(Fb(x_interps(1,:)',x_interps(2,:)'),im_sz(1),im_sz(2));
end

