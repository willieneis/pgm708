function [bCell,objVals] = hw04()
% Carries out question 4 for hw04

    load('hw4_prob4_input.mat');
    lamVec = [0.1,0.01,0.001,0.0001];
    bCell = cell(1,length(lamVec));
    for i=1:length(lamVec)
        bCell{i} = sparseGroupLasso(X,y,G,lamVec(i),lamVec(i));
        objVals(i) = evalObjective(X,y,G,lamVec(i),lamVec(i),bCell{i});
        save('hw04_results');
    end
