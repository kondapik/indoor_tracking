function [ X ] = nmds( proximities )

% extracting number of points from proximities matrix
N = size(proximities,1);

% creting identity matrix of size N x N
identity = eye(N);

% this is an equally sized matrix of 1's
one = ones(N);

% calculating centering matrix, commonly referred to as "J"
centering_matrix = identity - (1/N) * one;
J = centering_matrix;

% Applying double centering, i.e multiplying square of proximity matrix with J on each side and with a coefficient of -.5 
B = -.5*J*(proximities).*(proximities)*J;

% Extracting M eigen values and vectors, where M number of dimension to project
M = 2; % so we can plot in 2D
[eigvec1,eigvalue, ~] = svd(B);

% calculating largest M eigen values and vectors
eigvec = eigvec1(:,1:M); % Note that eigenvectors are in columns
eigval = eigvalue(1:M,1:M);
A = eigval;
A = A^(1/2);

% If we multiply eigenvectors by eigenvalues, we get coordinates, X:
eigvec = transpose(eigvec);
X = A*eigvec;
X = transpose(X);

%Adjestment Stage - iteratively reducing stess using LMA
epsilon = 10^-6;
Dist = sqrt(proximities);
iter = 5; %maximum number of iterations
% Setting options for lsqnonlin function to use MLA and optimization parameters
options = optimset('Display', 'off', ...
    'Algorithm', {'levenberg-marquardt',0.01},'MaxIter',5,'MaxFunEvals',100); 
total_cost=zeros(iter+1,1); % initializing cost function
for t=1:iter
    total_cost(t)=MDS_training_cost_total(X,Dist);
    % checking if error is less than a threshold
    if t>1 && abs(total_cost(t-1)-total_cost(t))/(total_cost(t)+eps)<epsilon 
        total_cost(t+1:end)=total_cost(t);
        return; % retrning values of cost is lessthan threshold
    end
    % setting order for adjesting coordinates of each point
    rand_order = randperm(N);
    count = 0;
    for k = rand_order
%    for k = 5:1
        count = count+1;
        x = X(k,:); % copying position of coordinate to be adjested
        nodes = [1:k-1 k+1:N]; 
        X_ = X(nodes,:); % array of remaining coordinates
        V_dist = Dist(nodes,k); % vector of distances
        x = lsqnonlin(@(x)MDS_cost_vector(x,V_dist,X_),x,[],[],options);
        X(k,:)=x; % replacing adjested position in original array
    end
end
total_cost(iter+1)=MDS_training_cost_total(X,Dist); % calculating cost
end

