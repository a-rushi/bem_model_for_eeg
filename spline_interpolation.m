function [A_potential] = spline_interpolation(A_vertices,B_vertices,B_potential)
% Input A_vertices are the vertices where you need to find potential (electrodes). It would have a dimension of 97x1. 
% Input B_vertices are the vertices of known potential points (1500). It would have a dimension of 1500x3.
% Input B_potential is the potential at B_vertices.
% Output A_potential is the interpolated potential at A_vertices.
x1 = B_vertices(:,1);  x2 = A_vertices(:,1); 
y1 = B_vertices(:,2);  y2 = A_vertices(:,2); 
z1 = B_vertices(:,3);  z2 = A_vertices(:,3); 
I = length(x1); 

%% Compute matrix P, Q
dsq = @(x11,y11,z11,x22,y22,z22)(((x11-x22)^2)+((y11-y22)^2)+((z11-z22)^2));
kfunc = @(dsq,w)((dsq^2)*(log10(dsq+w^2)));
w = 1;
for i = 1:I
    for j = 1:I
        dsqq = dsq(x1(i,1),y1(i,1),z1(i,1),x1(j,1),y1(j,1),z1(j,1));
        K(i,j) = kfunc(dsqq,w);
    end
end
for i = 1:I
E(i,:) = [1 x1(i,1) y1(i,1) x1(i,1)^2 x1(i,1)*y1(i,1) y1(i,1)^2 z1(i,1) z1(i,1)*x1(i,1) z1(i,1)*y1(i,1) z1(i,1)^2];
end
Amat = [K E; E' zeros(10,10)];
vmat = [B_potential; zeros(10,1)];
XXX = linsolve(Amat,vmat);
P = XXX(1:I,1);
Q = XXX(I+1:end,1);

%% Interpolation
clear K; clear E;
count = 0;
for i = 1:length(x2(:,1))
    for ii = 1:length(x2(1,:))
        count = count+1;
    for j = 1:I
        dsqq = dsq(x2(i,ii),y2(i,ii),z2(i,ii),x1(j,1),y1(j,1),z1(j,1));
        K(count,j) = kfunc(dsqq,w);
    end
    end
end
count = 0;
for i = 1:length(x2(:,1))
    for ii = 1:length(x2(1,:))
    count = count+1;
 E(count,:) = [1 x2(i,ii) y2(i,ii) x2(i,ii)^2 x2(i,ii)*y2(i,ii) y2(i,ii)^2 z2(i,ii) z2(i,ii)*x2(i,ii) z2(i,ii)*y2(i,ii) z2(i,ii)^2];
end
end
ans1 = K*P + E*Q;
A_potential = reshape(ans1,[length(x2(1,:)),length(x2(:,1))]);
A_potential = A_potential';