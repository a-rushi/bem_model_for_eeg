clear all;
close all;
clc;

%%
load bnd;
f = bnd(1);
bnd(1) = bnd(3);
bnd(3) = f;
% %% Scalp
sgm1 = 0.33;   % conductivity in S/m
sgm_air = 0 ;
figure;
patch('vertices',bnd(1).pos,'faces',bnd(1).tri,'edgecolor',[0 0 0],'facecolor',[0.5 1 1],'facelighting','none');
rotate3d;
axis equal
%% Skull
sgm2 = 0.0041;
% figure;
% patch('vertices',bnd(2).pos,'faces',bnd(2).tri,'edgecolor',[0 0 0],'facecolor',[0.5 1 1],'facelighting','none');
% rotate3d;
    %% Brain
    sgm3 = 0.33; 
% figure;
% patch('vertices',bnd(3).pos,'faces',bnd(3).tri,'edgecolor',[0 0 0],'facecolor',[0.5 1 1],'facelighting','none');
% rotate3d;


%% Electrode placement & Centroid calculation

% Electrode coordinates
load elec
coords = elec.elecpos;
I = length(coords);

% Source Placement
dip.pos = [5;-5;5]; 
dip.mom = dip.pos/(norm(dip.pos));

%%
size1= size(bnd(1).tri,1); 
size2= size(bnd(2).tri,1);
size3= size(bnd(3).tri,1);

centroid1=zeros((size1),3);
centroid2=zeros((size1),3);
centroid3=zeros((size1),3);

for j = 1:size1
    centroid1(j,1)=(bnd(1).pos(bnd(1).tri(j,1),1)+ bnd(1).pos(bnd(1).tri(j,2),1)+ bnd(1).pos(bnd(1).tri(j,3),1))/3;
    centroid1(j,2)=(bnd(1).pos(bnd(1).tri(j,1),2)+ bnd(1).pos(bnd(1).tri(j,2),2)+ bnd(1).pos(bnd(1).tri(j,3),2))/3;
    centroid1(j,3)=(bnd(1).pos(bnd(1).tri(j,1),3)+ bnd(1).pos(bnd(1).tri(j,2),3)+ bnd(1).pos(bnd(1).tri(j,3),3))/3;    
end
for j = 1:size2
    centroid2(j,1)=(bnd(2).pos(bnd(2).tri(j,1),1)+ bnd(2).pos(bnd(2).tri(j,2),1)+ bnd(2).pos(bnd(2).tri(j,3),1))/3;
    centroid2(j,2)=(bnd(2).pos(bnd(2).tri(j,1),2)+ bnd(2).pos(bnd(2).tri(j,2),2)+ bnd(2).pos(bnd(2).tri(j,3),2))/3;
    centroid2(j,3)=(bnd(2).pos(bnd(2).tri(j,1),3)+ bnd(2).pos(bnd(2).tri(j,2),3)+ bnd(2).pos(bnd(2).tri(j,3),3))/3;
end
for j = 1:size3
    centroid3(j,1)=(bnd(3).pos(bnd(3).tri(j,1),1)+ bnd(3).pos(bnd(3).tri(j,2),1)+ bnd(3).pos(bnd(3).tri(j,3),1))/3;
    centroid3(j,2)=(bnd(3).pos(bnd(3).tri(j,1),2)+ bnd(3).pos(bnd(3).tri(j,2),2)+ bnd(3).pos(bnd(3).tri(j,3),2))/3;
    centroid3(j,3)=(bnd(3).pos(bnd(3).tri(j,1),3)+ bnd(3).pos(bnd(3).tri(j,2),3)+ bnd(3).pos(bnd(3).tri(j,3),3))/3;
end
centroids= cat(1,centroid1,centroid2,centroid3);

%% Potential Calulation at electrodes

k1 = (2*sgm3)/(sgm1+sgm_air);
k2 = (2*sgm3)/(sgm2+sgm1);
k3 = (2*sgm3)/(sgm3+sgm2);

d_loc = dip.pos;
sigma = sgm3;   % brain conductivity
M = dip.mom;
%Io = random('normal',1.5,1,[1 100]);   % mean=1.5 and st.dev=1
% G1 = forward(coords1,d_loc,sigma);
% G2 = forward(coords2,d_loc,sigma);
% G3 = forward(coords3,d_loc,sigma);

G1 = forward(centroid1.',d_loc,sigma);
G2 = forward(centroid2.',d_loc,sigma);
G3 = forward(centroid3.',d_loc,sigma);

%% Signal Matrix 
time = (1:250)/250;           % manually create a time axis
S = sin(10*time*2*pi);

Vo1 = k1*(G1*M*S);   % 5000*t (t=1) potential on scalp at t different time instants
Vo2 = k2*(G2*M*S);   
Vo3 = k3*(G3*M*S);  

Vo = cat(1,Vo1, Vo2, Vo3); 
%% Area of triangle

area1=zeros(size1,3);
area2=zeros(size2,3);
area3=zeros(size3,3);
A = zeros(1,3);
B = zeros(1,3);
b1=0;
b2=0;
b3=0;

for j = 1:size1
    A(1,1)=bnd(1).pos(bnd(1).tri(j,1),1)- bnd(1).pos(bnd(1).tri(j,2),1);
    A(1,2)=bnd(1).pos(bnd(1).tri(j,1),2)- bnd(1).pos(bnd(1).tri(j,2),2);
    A(1,3)=bnd(1).pos(bnd(1).tri(j,1),3)- bnd(1).pos(bnd(1).tri(j,2),3);
            
    B(1,1)=bnd(1).pos(bnd(1).tri(j,1),1)- bnd(1).pos(bnd(1).tri(j,3),1);
    B(1,2)=bnd(1).pos(bnd(1).tri(j,1),2)- bnd(1).pos(bnd(1).tri(j,3),2);
    B(1,3)=bnd(1).pos(bnd(1).tri(j,1),3)- bnd(1).pos(bnd(1).tri(j,3),3);
            
    Z= cross(A,B);
    area1(j-b1,:) = Z/2;
%            area1(j,1)=.5*norm(Z);
end
for j = 1:size2
    A(1,1)=bnd(2).pos(bnd(2).tri(j,1),1)- bnd(2).pos(bnd(2).tri(j,2),1);
    A(1,2)=bnd(2).pos(bnd(2).tri(j,1),2)- bnd(2).pos(bnd(2).tri(j,2),2);
    A(1,3)=bnd(2).pos(bnd(2).tri(j,1),3)- bnd(2).pos(bnd(2).tri(j,2),3);
            
    B(1,1)=bnd(2).pos(bnd(2).tri(j,1),1)- bnd(2).pos(bnd(2).tri(j,3),1);
    B(1,2)=bnd(2).pos(bnd(2).tri(j,1),2)- bnd(2).pos(bnd(2).tri(j,3),2);
    B(1,3)=bnd(2).pos(bnd(2).tri(j,1),3)- bnd(2).pos(bnd(2).tri(j,3),3);
            
    Z= cross(A,B);
    area2(j-b2,:) = Z/2;
%   area2(j,1)=.5*norm(Z);
end
for j = 1:size3
    A(1,1)=bnd(3).pos(bnd(3).tri(j,1),1)- bnd(3).pos(bnd(3).tri(j,2),1);
    A(1,2)=bnd(3).pos(bnd(3).tri(j,1),2)- bnd(3).pos(bnd(3).tri(j,2),2);
    A(1,3)=bnd(3).pos(bnd(3).tri(j,1),3)- bnd(3).pos(bnd(3).tri(j,2),3);
            
    B(1,1)=bnd(3).pos(bnd(3).tri(j,1),1)- bnd(3).pos(bnd(3).tri(j,3),1);
    B(1,2)=bnd(3).pos(bnd(3).tri(j,1),2)- bnd(3).pos(bnd(3).tri(j,3),2);
    B(1,3)=bnd(3).pos(bnd(3).tri(j,1),3)- bnd(3).pos(bnd(3).tri(j,3),3);
            
    Z= cross(A,B);
    area3(j-b3,:) = Z/2;

end
area = cat(1,area1,area2,area3);

%% lamda matrix

 l11 = (sgm1 - sgm_air)/(sgm1 +sgm_air);
 l12 = (sgm1 - sgm_air)/(sgm2 +sgm1);
 l13 = (sgm1 - sgm_air)/(sgm3 +sgm2);
 l21 = (sgm2 - sgm1)/(sgm1 +sgm_air);
 l22 = (sgm2 - sgm1)/(sgm2 +sgm1);
 l23 = (sgm2 - sgm1)/(sgm3 +sgm2);
 l31 = (sgm3 - sgm2)/(sgm1 +sgm_air);
 l32 = (sgm3 - sgm2)/(sgm2 +sgm1);
 l33 = (sgm3 - sgm2)/(sgm3 +sgm2);
N = size1;
L11 = l11 * ones(N, N);
L12 = l12 * ones(N, N);
L13 = l13 * ones(N, N);
L21 = l21 * ones(N, N);
L22 = l22 * ones(N, N);
L23 = l23 * ones(N, N);
L31 = l31 * ones(N, N);
L32 = l32 * ones(N, N);
L33 = l33 * ones(N, N);
lamda = (1/2*pi)*[L11 L12 L13;L21 L22 L23;L31 L32 L33];

%% 
BH = B_helper(centroids,3*N, area);

%%
B = BH.*lamda;
V2 = inv(eye(8988)-B)*Vo;  % since there is no singulrity in inv(eye(8988)-B), let's not do deflation as that's the main purpose of doing deflation.
% C = B - 1/(3*N);
% AA= inv(1- C);
% 
% 
%% Calculating potential for only the brain for subtraction 
BH_dash = BH(((2*N) + 1):3*N,((2*N) + 1):3*N);
lamba_dash = (1/2*pi);
B_dash = lamba_dash*BH_dash;
%C_dash = B_dash - 1/(N);
k_dash = (sgm3+sgm2)/(sgm3);
Vo3_dash = k_dash*Vo3;
V_dash = inv(eye(2996) - B_dash)*Vo3_dash;
% 
%% IPA
beta = sgm2/sgm3;
V_ipa = isolated(Vo,V_dash,beta);
% 
%% Potential at electrodes
% V_final = AA*V_ipa;
% V_scalp = V_final(1:N,:);
%V_electrodes = spline_interpolation(coords, centroid1 , V_scalp );
V_final = inv(eye(8988)-B)*V_ipa;
V_scalp = V_final(1:N,:);

%% 
AA = weightedSparse(coords, centroid1, 1);

%%
V_electrodes = AA*V_scalp;


% %% Signal Matrix 
% time = (1:250)/250;           % manually create a time axis
% S = sin(10*time*2*pi);
% 
% V_time = V_electrodes*S;