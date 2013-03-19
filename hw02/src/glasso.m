function precOn = glasso(data,lambda)
% implementation of glasso
% pgm708 hw#2, problem 1.3

covmat = cov(data);
dim = size(covmat,1); maxIters = 100; eps = 0.001;
lastW = (lambda*eye(dim))+covmat; nextW = lastW;
iter = 1; endLoop=0;
while endLoop==0 && iter<maxIters
    for i = dim:-1:1
        inot = setdiff(1:dim,i);
        [V,D] = eig(nextW(inot,inot));
        d = diag(D);
        X = V*diag(sqrt(d))*V';
        Y = V*diag(1./sqrt(d))*V'*covmat(inot,i);
        w = my_lasso(X,Y,lambda);
        nextW(inot,i) = nextW(inot,inot)*w;
        nextW(i,inot) = nextW(inot,i)';
    end
    if norm(nextW-lastW,1) < eps, endLoop=1; end
    lastW = nextW;
    iter = iter+1;
end
prec = nextW^-1;
precOn = zeros(size(prec)); 
precOn(find(abs(prec)>0.01)) = 1;
