% PGM708, hw01, Willie Neiswanger
% sum product algorithm (Baum-Welch)

function marginalProb = sumProduct()
    % observations
    O_str = ['R','F','F','H','F','H','H','H','H','H', ...
            'H','H','H','R','H','H','H','R','H','H'];
    O = [2,3,3,1,3,1,1,1,1,1,...
        1,1,1,2,1,1,1,2,1,1];
    % transition matrix (col/row order S,M,A,W)
    T = [0.8,0.17,0.02,0.01;...
        0.1,0.7,0.19,0.01;...
        0.02,0.05,0.7,0.23;...
        0.2,0.01,0.04,0.7];
    % initial vector
    I = [0.15,0.6,0.2,0.05];
    % emission matrix (col order H,R,F)
    E = [0.4,0.3,0.3;...
        0.5,0.45,0.05;...
        0.3,0.4,0.3;...
        0.0001,0.2499,0.75];
    % do sumProduct
    marginalProb = sumProd();


    function M = sumProd()
        alpha = getAlphaMat();
        beta = getBetaMat();
        M = getMarginals(alpha,beta);
    end

    function A = getAlphaMat()
        A = I'.*E(:,O(1));
        for n = 2:length(O)
            A = nextAlpha(A);
        end
    end

    function matA = nextAlpha(matA)
        t = size(matA,2)+1;
        matA(:,t) = 0;
        for i = 1:size(matA,1)
            tmp = dot(matA(:,t-1),T(:,i));
            matA(i,t) = tmp*E(i,O(t));
        end
        matA(:,t) = matA(:,t) / sum(matA(:,t));
    end

    function B = getBetaMat()
        B = ones(length(I),1);
        for n = length(O)-1:-1:1
            B = nextBeta(B);
        end
    end

    function matB = nextBeta(matB)
        t = length(O) - size(matB,2);
        matB = [zeros(size(matB,1),1),matB];
        for i = 1:size(matB,1)
            matB(i,1) = sum(matB(:,2).*E(:,O(t+1)).*T(:,i));
        end
        matB(:,1) = matB(:,1) / sum(matB(:,1));
    end

    function marg = getMarginals(A,B)
        for i = 1:size(A,2)
            marg(:,i) = A(:,i).*B(:,i);
            marg(:,i) = marg(:,i)/sum(marg(:,i));
        end
    end

end
