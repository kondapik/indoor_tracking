clc
clear all
close all;
%%
% Step 1: loading proximity data and lables
load set2_1.mat;
% labels: 1 - tripod, 2&3&4 - beacons, 5 - tag
%%
% Step 2: Non Metric Multidimentional Scaling
[X] = nmds(proximities); % calling NMDS function
reliability = 0; % initializing the reliability parameter for NMDS output 
%%
% Step 3: calculating angles
% finding angles from 1:2 to 1:5 , 1:3 to 1:5 , 1:4 to 1:5 , 1:4 to 1:2 , 1:3 to 1:2 , 1:4 to 1:3 
angles(1) = rad2deg(atan2(det([X(5,:)-X(1,:);X(2,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(2,:)-X(1,:))));
angles(2) = rad2deg(atan2(det([X(5,:)-X(1,:);X(3,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(3,:)-X(1,:))));
angles(3) = rad2deg(atan2(det([X(5,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(5,:)-X(1,:),X(4,:)-X(1,:))));
angles(4) = rad2deg(atan2(det([X(2,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(2,:)-X(1,:),X(4,:)-X(1,:))));
angles(5) = rad2deg(atan2(det([X(2,:)-X(1,:);X(3,:)-X(1,:)]),dot(X(2,:)-X(1,:),X(3,:)-X(1,:))));
angles(6) = rad2deg(atan2(det([X(3,:)-X(1,:);X(4,:)-X(1,:)]),dot(X(3,:)-X(1,:),X(4,:)-X(1,:))));
% if the predicted positions are flipped, all calculated angles will be in reverse direction w.r.t the actual angles 
% Accordng to prefixed configuration, angles(4) is positive, hence it is used as benchmark to detect flipping 
if angles(4) < 0.0 % detecting flip
    angles(:) = -1*angles(:); % reversing angles back to original direction
end
%%
% Setting reliability value as discussed in NMDS section
if (abs(angles(4) - act_angles(4)) < 5 || abs(angles(5) - act_angles(5)) < 5 || abs(angles(6) - act_angles(6)) < 5)
     reliability = 1;
end
%%
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
