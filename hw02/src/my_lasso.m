function w = my_lasso(X,Y,lambda)
% lasso implementation
% pgm708 hw#2, problem 1.3

maxIters = 100; eps = 0.001; [n,dim] = size(X);
if dim > n, w = zeros(p,1); else w = X \ Y; end
endLoop = 0; iter = 1; w_old = w; 
prod1 = X'*X; prod2 = X'*Y;
while endLoop==0 && iter<maxIters,
    for i = 1:dim,
        inot = setdiff(1:dim,i); S = prod1(i,inot)*w(inot)-prod2(i);
        if S > lambda, w(i) = (lambda-S)/norm(X(:,i),2)^2;
        elseif S < -lambda, w(i) = -(lambda+S)/norm(X(:,i),2)^2;
        else w(i) = 0; end
    end
    if norm(w-w_old,1) < eps, endLoop=1; end;
    iter = iter+1; w_old = w;
end
