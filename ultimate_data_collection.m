clc
clear all
close all;
%%
iter = 5; %maximum number of iterations
angles = zeros(100,101,iter+1,10); %initializing angles, values in angles (data_point,loop,iterations,angles or reliability)
for poi = 1:100
    % loading proximites matrix and actual angles
    load (sprintf('set2_%i.mat',poi));
    i = int16(1);
    %angles = zeros(100,11,6);
    % running algorithm 100 times on the data point 
    while i < 101
    % code from NMDS function and MDS   
    X = zeros(2,5);
    N = size(proximities,1);
    identity = eye(N);
    one = ones(N);
    centering_matrix = identity - (1/N) * one;
    J = centering_matrix;
    B = -.5*J*(proximities).*(proximities)*J;
    M = 2;
    [eigvec1,eigvalue, ~] = svd(B);
    eigvec = eigvec1(:,1:M);
    eigval = eigvalue(1:M,1:M);
    A = eigval;
    A = A^(1/2);
    eigvec = transpose(eigvec);
    X = A*eigvec;
    X = transpose(X);
    angles(poi,i,1,1) = rad2deg(atan2(abs(det([X(5,:)-X(1,:);X(2,:)-X(1,:)])),dot(X(5,:)-X(1,:),X(2,:)-X(1,:))));
    angles(poi,i,1,2) = rad2deg(atan2(abs(det([X(5,:)-X(1,:);X(3,:)-X(1,:)])),dot(X(5,:)-X(1,:),X(3,:)-X(1,:))));
    angles(poi,i,1,3) = rad2deg(atan2(abs(det([X(5,:)-X(1,:);X(4,:)-X(1,:)])),dot(X(5,:)-X(1,:),X(4,:)-X(1,:))));
    angles(poi,i,1,4) = rad2deg(atan2(abs(det([X(2,:)-X(1,:);X(4,:)-X(1,:)])),dot(X(2,:)-X(1,:),X(4,:)-X(1,:))));
    angles(poi,i,1,5) = rad2deg(atan2(abs(det([X(2,:)-X(1,:);X(3,:)-X(1,:)])),dot(X(2,:)-X(1,:),X(3,:)-X(1,:))));
    angles(poi,i,1,6) = rad2deg(atan2(abs(det([X(3,:)-X(1,:);X(4,:)-X(1,:)])),dot(X(3,:)-X(1,:),X(4,:)-X(1,:))));
    epsilon = 10^-6;
    Dist = sqrt(proximities);    
    options = optimset('Display', 'off', ...
        'Algorithm', {'levenberg-marquardt',0.01},'MaxIter',5,'MaxFunEvals',100);
    total_cost=zeros(iter+1,1);
    for t=1:iter
        total_cost(t)=MDS_training_cost_total(X,Dist);
        rand_order = randperm(N);
        count = 0;
        for k = rand_order
            count = count+1;
            x = X(k,:);
            nodes = [1:k-1 k+1:N];
            X_ = X(nodes,:);
            V_dist = Dist(nodes,k);
            x = lsqnonlin(@(x)MDS_cost_vector(x,V_dist,X_),x,[],[],options);
            X(k,:)=x;
        end
        angles(poi,i,1+t,1) = rad2deg(atan2(det([X(5,:)-X(1,:);X(2,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(2,:)-X(1,:))));
        angles(poi,i,1+t,2) = rad2deg(atan2(det([X(5,:)-X(1,:);X(3,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(3,:)-X(1,:))));
        angles(poi,i,1+t,3) = rad2deg(atan2(det([X(5,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(4,:)-X(1,:))));
        angles(poi,i,1+t,4) = rad2deg(atan2(det([X(2,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(2,:)-X(1,:),X(4,:)-X(1,:))));
        angles(poi,i,1+t,5) = rad2deg(atan2(det([X(2,:)-X(1,:);X(3,:)-X(1,:)]),dot(X(2,:)-X(1,:),X(3,:)-X(1,:))));
        angles(poi,i,1+t,6) = rad2deg(atan2(det([X(3,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(3,:)-X(1,:),X(4,:)-X(1,:))));
    end
    %     if (abs(min(angles(i,iter+1,:))>180 || abs(max(angles(i,iter+1,:)))>180))
    %     X(:,1) = -X(:,1);
    %     end
    %     angles(i,iter+1,1) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(2,1),X(2,2));
    %     angles(i,iter+1,2) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(3,1),X(3,2));
    %     angles(i,iter+1,3) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(4,1),X(4,2));
    %     angles(i,iter+1,4) = angle(X(1,1),X(1,2),X(2,1),X(2,2),X(4,1),X(4,2));
    %     angles(i,iter+1,5) = angle(X(1,1),X(1,2),X(2,1),X(2,2),X(3,1),X(3,2));
    %     angles(i,iter+1,6) = angle(X(1,1),X(1,2),X(3,1),X(3,2),X(4,1),X(4,2));

    if angles(poi,i,iter+1,4) < 0.0
          angles(poi,i,iter+1,:) = -1*angles(poi,i,iter+1,:);
    end

    % Filtering with and condition, i.e checking if all angles are correct if only beacons are correct w.r.t actual angles
    if (abs(angles(poi,i,1+iter,4) - act_angles(4)) < 5 && abs(angles(poi,i,1+iter,5) - act_angles(5)) < 5 && abs(angles(poi,i,1+iter,6) - act_angles(6)) < 5)
          angles(poi,i,1+iter,7) = 1;
    end
    if (abs(angles(poi,i,1+iter,1) - act_angles(1)) < 5 && abs(angles(poi,i,1+iter,2) - act_angles(2)) < 5 && abs(angles(poi,i,1+iter,3) - act_angles(3)) < 5 && abs(angles(poi,i,1+iter,4) - act_angles(4)) < 5 && abs(angles(poi,i,1+iter,5) - act_angles(5)) < 5 && abs(angles(poi,i,1+iter,6) - act_angles(6)) < 5)
          angles(poi,i,1+iter,8) = 1;
    end
   
    % Filtering with or condition, i.e checking if all angles are correct if only beacons are correct w.r.t actual angles
    if (abs(angles(poi,i,1+iter,4) - act_angles(4)) < 5 || abs(angles(poi,i,1+iter,5) - act_angles(5)) < 5 || abs(angles(poi,i,1+iter,6) - act_angles(6)) < 5)
          angles(poi,i,1+iter,9) = 1;
    end
    if (abs(angles(poi,i,1+iter,1) - act_angles(1)) < 5 || abs(angles(poi,i,1+iter,2) - act_angles(2)) < 5 || abs(angles(poi,i,1+iter,3) - act_angles(3)) < 5 || abs(angles(poi,i,1+iter,4) - act_angles(4)) < 5 || abs(angles(poi,i,1+iter,5) - act_angles(5)) < 5 || abs(angles(poi,i,1+iter,6) - act_angles(6)) < 5)
          angles(poi,i,1+iter,10) = 1;
    end
        i = i+1;
    end
    
    % calulating preformance of the code w.r.t acutal angles
    i = int16(1);
    angles(101,:,:) = 0;
    while i < 101
        k = int16(1);
        while k < iter + 2
            % incrementing the value in 101th row if angle value is correct
            if abs(angles(poi,i,k,1) - act_angles(1)) < 5
                angles(poi,101,k,1) = 1 + angles(poi,101,k,1);
            end
            if abs(angles(poi,i,k,2) - act_angles(2)) < 5
                angles(poi,101,k,2) = 1 + angles(poi,101,k,2);
            end
            if abs(angles(poi,i,k,3) - act_angles(3)) < 5
                angles(poi,101,k,3) = 1 + angles(poi,101,k,3);
            end
            if abs(angles(poi,i,k,4) - act_angles(4)) < 5
                angles(poi,101,k,4) = 1 + angles(poi,101,k,4);
            end
            if abs(angles(poi,i,k,5) - act_angles(5)) < 5
                angles(poi,101,k,5) = 1 + angles(poi,101,k,5);
            end
            if abs(angles(poi,i,k,6) - act_angles(6)) < 5
                angles(poi,101,k,6) = 1 + angles(poi,101,k,6);
            end
            % incrementing the value in 101th row if reliabiltity value is 1
            if angles(poi,i,k,7) == 1
                angles(poi,101,k,7) = 1 + angles(poi,101,k,7);
            end
            if angles(poi,i,k,8) == 1
                angles(poi,101,k,8) = 1 + angles(poi,101,k,8);
            end
            if angles(poi,i,k,9) == 1
                angles(poi,101,k,9) = 1 + angles(poi,101,k,9);
            end
            if angles(poi,i,k,10) == 1
                angles(poi,101,k,10) = 1 + angles(poi,101,k,10);
            end
            k = k + 1;
        end
        i = i + 1;
    end
    % %%
    % i = int16(1);
    % angles(101,:,:) = 0;
    % while i < 101
    %     k = int16(1);
    %     while k < 12
    %         if abs(angles(i,k,1)- 45) < 5
    %             angles(101,k,1) = 1 + angles(101,k,1);
    %         end
    %         if abs(angles(i,k,2)- 0) < 5
    %             angles(101,k,2) = 1 + angles(101,k,2);
    %         end
    %         if abs(angles(i,k,3)- 45) < 5
    %             angles(101,k,3) = 1 + angles(101,k,3);
    %         end
    %         if abs(angles(i,k,4)- 90) < 5
    %             angles(101,k,4) = 1 + angles(101,k,4);
    %         end
    %         if abs(angles(i,k,5)- 45) < 5
    %             angles(101,k,5) = 1 + angles(101,k,5);
    %         end
    %         if abs(angles(i,k,6)- 45) < 5
    %             angles(101,k,6) = 1 + angles(101,k,6);
    %         end
    %         k = k + 1;
    %     end
    %     i = i + 1;
    % end
    %%
    % Plotting
    % subplot(2,3,1);
    % mesh(angles(:,:,1));
    % subplot(2,3,2);
    % mesh(angles(:,:,2));
    % subplot(2,3,3);
    % mesh(angles(:,:,3));
    % subplot(2,3,4);
    % mesh(angles(:,:,4));
    % subplot(2,3,5);
    % mesh(angles(:,:,5));
    % subplot(2,3,6);
    % mesh(angles(:,:,6));
    %%
    % %%
    % % Step 3: calculating angles
    % % finding angles from 1:2 to 1:5 and 1:3 to 1:5 and 1:4 to 1:5
    % angles(1) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(2,1),X(2,2));
    % angles(2) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(3,1),X(3,2));
    % angles(3) = angle(X(1,1),X(1,2),X(5,1),X(5,2),X(4,1),X(4,2));
    % angles(4) = angle(X(1,1),X(1,2),X(2,1),X(2,2),X(4,1),X(4,2));
    % if angles(4) < 0.0
    %     angles = -1*angles;
    % end
    % %%
    % % Plotting the points 
    % plot(X(:,1),X(:,2),'o')
    % title('Result of Non Metric MDS using ILMA');
    % % labeling
    % text(X(:,1), X(:,2), label, 'VerticalAlignment','bottom', ...
    %                              'HorizontalAlignment','right')
    % %%
    % % Plotting stress
    % figure;
    % plot(total_cost,'-b^','LineWidth',2);
    % title('raw stress in each training iteration');
    % xlabel('iteration');
    % ylabel('raw stress');
end
