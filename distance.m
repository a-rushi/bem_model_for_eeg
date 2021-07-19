function [d] = distance(v1, v2)
xdi = v1(1) - v2(1);
ydi = v1(2) - v2(2);
zdi = v1(3) - v2(3);

d = sqrt(xdi.^2 + ydi.^2 + zdi.^2);