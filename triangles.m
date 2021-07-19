clear all;
close all;
clc;

load bnd

vertices= bnd(2).pos;
faces= bnd(2).tri;
P = length(faces);
a=[];
for i= 1:P
    if faces(i,1) == faces(i,2)
        a=[a;0];
    end
end
