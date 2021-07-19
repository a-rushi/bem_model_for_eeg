function [w] = weight(v1, v2, tau)
%v1 and v2 are supposed to be 3*1 vectors containing the coordinates
d = distance(v1, v2);
w = exp(-(d.^2)/(2*(tau.^2)));