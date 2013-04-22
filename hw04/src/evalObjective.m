function objVal = evalObjective(X,y,G,l1,l2,b)
% evaluate sparseGroupLasso objective function

    t1 = 0.5*sum(y-X*b');
    
    t2 = l1*sum(abs(b)); 
    for g=1:size(G,2)
       sqrtGpVec(g) = sqrt(sum(b().^2));
    end

    t3 = l2*sum(sqrtGpVec);

    objVal = t1+t2+t3;
