function precOn = MB(data,lambda) 
% implementation of MB
% pgm708 hw#2, problem 1.3

precOn = zeros(size(data,2));
for i=1:size(data,2)
    w = my_lasso(data,data(:,i),lambda);
    onInd = find(abs(w)>0);
    for j=1:length(onInd)
        precOn(i,i) = 1;
        precOn(i,onInd(j)) = 1;
        precOn(onInd(j),i) = 1;
    end
end
