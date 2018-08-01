clear all;
A = imread('/Users/thiyaga/Documents/thiy/Codes/Camera_Test/build/Debug/VRightBig/31.bmp');
G = rgb2gray(A);
G = double(G)/255; % Scaling between 0 - 1

[H,W]=size(G);
OW = 1008;
OH = 1008;
fovFactor = 1.0;
params = struct('fx',629.2711604878951, 'fy', 628.5435111583766, 'cx', 500.5283689158744, 'cy', 488.5061520580846, 'k1', -0.4081, 'k2', 0.21765, 'k3',  -0.06604786876627262, 'p1', -0.00071, 'p2', -0.00014);
fx = params.fx/W;
fy = params.fy/H;
cx = params.cx/W;
cy = params.cy/H;
k1 = params.k1;
k2 = params.k2;
k3 = params.k3;
p1 = params.p1;
p2 = params.p2;

K = [fx 0 cx; 0 fy cy; 0 0 1];


I = zeros(OH,OW);

[i j] = find(~isnan(I));

% Xp = the xyz vals of points on the z plane
Xp = inv(K)*[j i ones(length(i),1)]';

% Now we calculate how those points distort i.e forward map them through the distortion

x = Xp(1,:)/OW;
y = Xp(2,:)/OH;
x = (x - 0.5) * fovFactor + 0.5;
y = (y - 0.5) * fovFactor + 0.5;

r2 = Xp(1,:).^2+Xp(2,:).^2;
x = x.*(1+k1*r2 + k2*r2.^2 + k3*r2.^4) + 2*p1.*x.*y + p2*(r2 + 2*x.^2);
y = y.*(1+k1*r2 + k2*r2.^2 + k3*r2.^4) + 2*p2.*x.*y + p1*(r2 + 2*y.^2);

% u and v are now the distorted cooridnates
u = reshape(fx*x + cx,size(I));
v = reshape(fy*y + cy,size(I));


I = interp2(G, u, v);
 subplot(121); imshow(G);%imagesc(G);
 subplot(122); imshow(I);%imagesc(I);

