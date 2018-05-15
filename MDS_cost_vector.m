
% The target function to be minimized at each step. 

function f=MDS_cost_vector(x,V_dist,X)

X=bsxfun(@minus,X,x);
V=sum(X.*X,2);
f=sqrt(V)-V_dist;