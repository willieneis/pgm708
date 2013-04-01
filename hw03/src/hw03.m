function hw03()

% parse data
load('hw3data.mat');
alpha = 0.1*ones(1,size(beta_matrix,2));
tic; out = ldavb(data,beta_matrix,alpha); time=toc;
fprintf('Computation took %g seconds.\n',time);
disp('Iterations for convergence (per document):');
disp(mat2str(out{3}));

% figures
figure,imagesc(out{2}{1}); colorbar;
figure,imagesc(out{1}); colorbar;

% last part
alphas = [0.01,1,10];
for i=1:length(alphas)
    out = ldavb(data,beta_matrix,alphas(i));
    figure,imagesc(out{1}); colorbar;
    fprintf('alpha = %g. Mean number of iterations for convergence: %g.\n',alphas(i),mean(out{3}));
end
