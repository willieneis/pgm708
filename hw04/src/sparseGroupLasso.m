function b = sparseGroupLasso(X,y,G,l1,l2)
% Carries out sparse group lasso on (n x p) data matrix X,
%   (n x 1) label vector y, group index matrix G (with 
%   format specified in hw4_prob4_input_readme.txt), and
%   lambda regularization params l1 and l2.
   
    X = X-repmat(mean(X),size(X,1),1);
    tol = 0.005;
    [n,p] = size(X);

    % init b
    %b = zeros(1,p);
    b = randn(1,p);

    % sparse group lasso
    for iter=1:1000
        bOld = b;
        for g=1:size(G,2)
            r = compute_r(g); 
            minJ = findMinJ(g,r);
            if minJ<=1
                b(G(1,g):G(2,g)) = 0;
            else
                b(G(1,g):G(2,g)) = findFinalParams(g,r);
            end
            fprintf('finished: iter=%d g=%d lam1=%g lam2=%g\n',iter,g,l1,l2);
        end
        if max(abs(b-bOld))<tol
            break;
        end
        fprintf('change=%g\n',max(abs(b-bOld)));
    end

    % ----------
    function res = compute_r(gInd)
        % compute residual value r
        XNew = X; XNew(:,G(1,gInd):G(2,gInd)) = [];
        bNew = b; bNew(G(1,gInd):G(2,gInd)) = [];
        res = y-(XNew*bNew');
    end

    function minVal = findMinJ(gInd,res)
        % compute minimum of (first) J function
        Z = X(:,G(1,gInd):G(2,gInd));
        a = Z'*res; % k x 1
        tHat = a/l2;
        tHat(find(tHat>1)) = 1;
        tHat(find(tHat<-1)) = -1;
        minVal = (1/l1^2)*sum((a-(l2*tHat)).^2);
    end

    function theta = findFinalParams(gInd,res);
        % coordinate descent to find final set of params
        theta = b(G(1,gInd):G(2,gInd));
        Z = X(:,G(1,gInd):G(2,gInd));
        for iter2=1:1000
            theta_old = theta;
            for j=1:length(theta)
                w_j = res-(Z*theta')+(Z(:,j)*theta(j)); % n x 1
                ZwNoAbs = dot(Z(:,j),w_j);
                Zw = abs(ZwNoAbs);
                if Zw<l2
                    theta(j) = 0;
                else
                    %fHandle = @(x)findJ2Val(j,theta,Z,res,x);
                    %opt = optimset('fminsearch');
                    %opt.MaxFunEvals = 225;
                    %theta(j) = fminsearch(fHandle,0,opt);
                    theta(j) = wthresh(ZwNoAbs,'s',l2)/(dot(Z(:,j),Z(:,j))+(l2/norm(theta)));
                end
                %fprintf(['\tinner loop finished: iter=%d, j=%d, theta=[',repmat('%g,',1,length(theta)),']\n'],iter2,j,theta);
            end
            if max(abs(theta-theta_old))<tol
                break;
            end
            %fprintf('\tchange=%g\n',max(abs(theta-theta_old)))
        end
    end

    function J2Val = findJ2Val(jInd,th,Zmat,resVec,th_j_val)
        % return J2 function value (where th(jInd) replaced with th_j_val)
        th(jInd) = th_j_val;
        diffVals = resVec-Zmat*th';
        t1 = 0.5*sum(diffVals.^2);
        t2 = l1*norm(th); 
        t3 = l2*sum(abs(th));
        J2Val = t1+t2+t3;
    end

end
