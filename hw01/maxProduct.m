% PGM708, hw01, Willie Neiswanger
% max product algorithm (Viterbi)

function mostLikelySeq = sumProduct()
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
    % do maxProduct
    mostLikelySeq = maxProd();


    function ind = maxProd()
        ind = [];
        W = log(I) + log(E(:,O(1))');
        for n = 2:length(O)
            [W,ind] = nextW(W,ind);
        end
        [todel,ind(:,end+1)] = max(W);
    end

    function [mat,indMat] = nextW(mat,indMat)
        t = size(indMat,2)+2;
        K = size(mat,2);
        indMat(1:K,end+1) = zeros(K,1);
        for i = 1:K
            [maxVal,indMat(i,end)] = max(log(T(:,i)) + mat');
            mat(i) = log(E(i,O(t))) + maxVal;
        end
    end

    function whichInd = getMax(wVec,t)
        whichInd = max(log() + wVec); 
    end 

end
