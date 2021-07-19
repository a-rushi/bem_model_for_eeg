function [G] = forward(coords,d_loc,sigma)
%coords are sensor location
%d_loc is source location
%sigma is head conductivity
P = length(d_loc(1,:));    % length of first row of d_loc (P=1) = no. of dipoles
I = size(coords,2);     % second dim of coords (I=5000)
G = [];
for i = 1:P
    g_p = [];
    for j = 1:I
            b(1,:) = (coords(:,j)-d_loc(:,i))/(4*pi*sigma*((norm(coords(:,j)-d_loc(:,i)))^3));      
        g_p = [g_p; b];
    end   % g_p is 64*3  
    G = [G,g_p];   % G = g_p for 1 dipole
end
end