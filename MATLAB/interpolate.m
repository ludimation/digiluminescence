function [im1w,xsamp,ysamp,xout,yout] = interpolate(im1,im2,p,q,method,factor,showfigs);
% Provided by Professor Radke as example of 
%     Scattered Data Interpolation

if nargin<7
    showfigs = 0;
end

% p, q are N x 2;

im1 = double(imresize(im1,1/factor));
im2 = double(imresize(im2,1/factor));

p = p/factor;
q = q/factor;

[h,w,N] = size(im1);
[xin, yin] = meshgrid(1:w,1:h);

corners = [1 1; 1 w; h/2 w; h w; h 1; h w/2; 1 w/2; h/2 1];
corners = corners(:,[2 1]);

p = [p;corners];
q = [q;corners];

% 1st col of p = columns
% 2nd col of p = rows

n = length(p);

close all

sp = 2;  % pixel grid spacing

xout = [];
yout = [];
z = [];
beta = 1e-4;

ivec = 1:sp:h;
jvec = 1:sp:w;

[xsamp, ysamp] = meshgrid(jvec,ivec);

switch method
    
    case 'tps'
        
        phisub = zeros(n);
        for i=1:n
           for j=1:n
                if i ~= j
                    r = norm(p(i,:) - p(j,:));
                    phisub(i,j) = r^2*log(r);
                else
                    phisub(i,j) = 0;
                end
            end
        end
        
        phimat = [phisub, p, ones(n,1); p', zeros(2,3); ones(1,n) zeros(1,3)];
        rhs = [q;zeros(3,2)];
        
        wab = phimat\rhs;

        % i = which row
        % j = which col
        
        for ii=1:length(ivec);
            for jj=1:length(jvec);
                
                i = ivec(ii);
                j = jvec(jj);
                
                for k=1:n
                    r = norm(p(k,:) - [j i]);
                    if r>0
                        b(k) = r^2*log(r);
                    else
                        b(k) = 0;
                    end
                end

                xout(ii,jj) = [b j i 1]*wab(:,1);
                yout(ii,jj) = [b j i 1]*wab(:,2);

                for k=1:3
                    z{k}(ii,jj) = im1(i,j,k);
                end
            end
        end
        
        im1w = delinterp(xout,yout,z,xin,yin);            
        
    case 'arap'  % See Schaefer
        
        for ii=1:length(ivec);
            for jj=1:length(jvec);
                
                i = ivec(ii);
                j = jvec(jj);
               
                % i = which row
                % j = which col
                
                % 1st col of p = which column
                % 2nd col of p = which row
                
                if (jj == length(jvec))
                    disp('');
                end
                
                for k=1:n
                    wt(k) = (beta + norm([j,i] - p(k,:)))^(-2);
                end
                
                wtsum = sum(wt);
                wtmat = wt'*[1 1];

                pstar = sum(wtmat.*p)/wtsum;
                qstar = sum(wtmat.*q)/wtsum;

                phat = p - ones(n,1)*pstar;
                qhat = q - ones(n,1)*qstar;

%                    mu = sum(diag(phat'*diag(wt)*phat));  
                % for similarity transformations

                A = zeros(2);
                mu1 = 0;
                mu2 = 0;

                for k=1:n
                    A = A + wt(k)*[qhat(k,1) qhat(k,2); -qhat(k,2) qhat(k,1)]*[phat(k,1) -phat(k,2); phat(k,2) phat(k,1)];
                    mu1 = mu1 + wt(k)*(qhat(k,1)*phat(k,1) + qhat(k,2)*phat(k,2));
                    mu2 = mu2 + wt(k)*(qhat(k,2)*phat(k,1) - qhat(k,1)*phat(k,2));
                end

                mu = sqrt(mu1^2 + mu2^2);

                A = A/mu;

                t = qstar' - A*pstar';

                if any(isnan([A(:); t]))
                    disp('Nan!');
                end

                % Now q = A*p + t is the forward transformation.

                % i = which row
                % j = which col
                
                % 1st col of p = which column
                % 2nd col of p = which row
                
                fwd = A*[j;i] + t;

                % 1st elt of fwd = which column
                % 2nd elt of fwd = which row     
                
                xout(ii,jj) = fwd(1);
                yout(ii,jj) = fwd(2);

                % i = which row
                % j = which col
                
                for k=1:3
                    z{k}(ii,jj) = im1(i,j,k);
                end
   
            end
        end
        
        im1w = delinterp(xout,yout,z,xin,yin);
       
end

im1w = uint8(im1w);


if showfigs
    showpts(im1,p,1);
    showpts(im2,q,2);
    showpts(im1w,q,3);
    showmesh(im1,xsamp,ysamp,4);
    showmesh(im2,xout,yout,5);
end

function im1w = delinterp(xout,yout,z,xin,yin)

dt = DelaunayTri(xout(:),yout(:));

for k=1:3
    F = triscatteredinterp(dt,z{k}(:),'natural');
    im1w(:,:,k) = F(xin,yin);
end

% 
%         figure(1)
%         clf;
%         imshow(uint8(im1));
%         hold on
% %        plotmesh(xsamp,ysamp);
%         plot(p(:,1),p(:,2),'g.','markersize',25)
% 
%         figure(2)
%         clf;
%         imshow(uint8(im2));
%         hold on
%  %       plotmesh(xout,yout);
%         plot(q(:,1),q(:,2),'g.','markersize',25)
        
%         sm = zeros(length(p)-8,6);

        
%         figure(4)
%         clf
%         hold on
%         axis equal
%         set(gca,'ydir','reverse')
%         for i=1:size(xout,1)
%             for j=1:size(xout,2)
%                 col = [z{1}(i,j) z{2}(i,j) z{3}(i,j)];
%                 col = col/255;
%                 col(col<0) = 0;
%                 col(col>1) =1;
%                 plot(xout(i,j),yout(i,j),'.','markersize',25,'color',col);
%             end
%         end
%         triplot(dt);