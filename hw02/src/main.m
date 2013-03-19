function main()

% load data
load('Xinput.mat');

% MB algo
lambdaMb = [0,20,30,40];
for i=1:length(lambdaMb)
    precMb{i} = MB(X,lambdaMb(i));
    figure, imagesc(precMb{i}); colormap(gray);
end

% glasso
lambdaGlasso = [0,0.2,0.5,0.8];
precGlasso = {};
for i=1:length(lambdaGlasso)
    precGlasso{i} = glasso(X,lambdaGlasso(i));
    figure, imagesc(precGlasso{i}); colormap(gray);
end

