% The total raw stress.

function f=MDS_training_cost_total(X,Dist)
% X: matrix of MDS codes
% Dist: distance matrix
% f: raw stress

N=size(Dist,1);
f=0;
for i=1:N-1
    x=X(i,:);
    XX=X(i+1:N,:);
    V_dist=Dist(i+1:N,i);
    
    XX=bsxfun(@minus,XX,x);
    V=sum(XX.*XX,2);
    V=sqrt(V)-V_dist;
    f = f + sum(V.*V);
end
