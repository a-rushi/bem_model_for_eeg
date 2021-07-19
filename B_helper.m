function [G] = B_helper(coords,N,area)
% coords are centroids of all triangles
% coords1 are centroids of triangles of surf1
% area is the matrix of area of all triangles (15000*3)
% N = no. of triangles = N1+N2+N3
G = [];
for i = 1:N
    g_p = [];
    for j = 1:N
            if ~(i==j)
                b(1,:) = (coords(j,:)-coords(i,:))/((norm(coords(j,:)-coords(i,:)))^3);
                p = dot(b,area(j,:));
                g_p = [g_p, p];
            else 
                g_p = [g_p, 0];
            end
    end
    G = [G; g_p];  
end
end