function [P] = isolated(V,V_k,beta)
S = (2*size(V,1))/3;
g=V*beta;
k = (2*beta)/(beta+1);
for i = (S+1):size(V,1)
    g(i) = g(i) - k*(V_k(i-S));
end
P= g;
end