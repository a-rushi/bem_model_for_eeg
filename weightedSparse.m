function [D] = weightedSparse(coords, centroid, tau)

%this function is to be applied only on the same surface
sizcor = size(coords);
sizcen = size(centroid);
D = [];
for i = 1:sizcor
    v = [coords(i,1) coords(i,2) coords(i,3)];
    w = [];
    for j = 1:sizcen
        w = [w weight(v, centroid(j,:), tau)];
    end
    totalw = sum(w);
    w = w./totalw;
    D = [D; w];
end

