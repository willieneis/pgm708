Please load input data in MATLAB.

>> load('hw4_prob4_input.mat');

Then you will see input matrix X, output vector y, and group indices G.

X: N by P matirx, where N is sample size and P is the number of features.
y: N by 1 matirx, where N is sample.
G: Group indices. Here number of columns is the number of groups.
   Each column indicates the start and end position of each group.
   (E.g.) G(:,i) = [j,k]' means that i-th group starts from j-th feature to
k-th feature.


