function [X_out,P_out] = ekf_6dof2(X, P, noise_gyro, noise_accel,...
                       gyroX,gyroY,gyroZ,...
                       accX,accY,accZ,...
                       dt)
% function y = ekf_6dof2(gyroX,gyroY,gyroZ,accX,accY,accZ,...
%                         P,Q,R...
%                          dt ,TX)
%     X=[q0;q1;q2;q3]
%     P 4X4:
%     Q 4X3:
%     R:
%     meas 3X3: gyroX,gyroY,gyroZ,...        %������ԭʼֵ
%               accX,accY,accZ,...         %���ٶ�ԭʼֵ m/s^2
%               magX,magY,magZ,...         %�ش�ԭʼֵ
%     dt:              %����ʱ��
 X_p = zeros(size(X));    %һ��Ԥ��״̬��X_p, ����ռ�
 P_p = zeros(size(P));    %һ��Ԥ��Э�������P_p ������ռ�
 F  = zeros(size(P));
 gravity = 9.81007;          % Gravity magnitude, m/s^2
%״̬ת�ƾ���..��
F= [1  , -0.5*dt*gyroX   ,-0.5*dt*gyroY  , -0.5*dt*gyroZ;...
    0.5*dt*gyroX , 1  ,  0.5*dt*gyroZ  ,  -0.5*dt*gyroY;...
    0.5*dt*gyroY , -0.5*dt*gyroZ , 1  ,   0.5*dt*gyroX;...
    0.5*dt*gyroZ ,   0.5*dt*gyroY ,   -0.5*dt*gyroX,  1];
%����������
G = [-X(2)  -X(3)  -X(4); ...
      X(1)  -X(4)   X(3); ...
      X(4)   X(1)  -X(2); ...
     -X(3)   X(2)   X(1)] / 2;
%% һ��Ԥ��
X_p = F*X;
X_p = X_p / norm(X_p);    %��Ԫ����һ��
%Q = (noise_gyro * dt)^2 * (G * G');
Q = G *(noise_gyro* dt )^2 * ( G');
P_p = F * P * F'+ Q;
%%
%������
a =[accX; accY; accZ];
% [ -2*q2,  2*q3, -2*q0, 2*q1]
% [  2*q1,  2*q0,  2*q3, 2*q2]
% [  2*q0, -2*q1, -2*q2, 2*q3]
H =[ -2*X_p(3), 2*X_p(4), -2*X_p(1), 2*X_p(2);...
      2*X_p(2), 2*X_p(1),  2*X_p(4), 2*X_p(3);...
      2*X_p(1), -2*X_p(2), -2*X_p(3),2*X_p(4)];
  
%��Ϣ
innvo = a/norm(a) - H*X_p;

%R
R_internal = (noise_accel / norm(a))^2 * eye(3);
R_external = (1-gravity/norm(a))^2 * eye(3);
R = R_internal + R_external;

%����������kg
Kg = P_p * H'/(H*P_p*H'+R); 

%�������X��P
X_out = X_p + Kg* innvo;
%P_P_est = P_P_prd -P_P_prd * H;    eye(4)
P_out = (eye(4) - Kg * H) * P_p;



% 
% %״̬ת�ƾ���..��
% A= [1  , -0.5*dt*gyroX   ,-0.5*dt*gyroY  , -0.5*dt*gyroZ;...
%     0.5*dt*gyroX , 1  ,  0.5*dt*gyroZ  ,  -0.5*dt*gyroY;...
%     0.5*dt*gyroY , -0.5*dt*gyroZ , 1  ,   0.5*dt*gyroX;...
%     0.5*dt*gyroZ ,   0.5*dt*gyroY ,   -0.5*dt*gyroX,  1];
% %Q
% %R
% %--- Dataset parameters
% %deltat = 1/200;             % Sampling period
% %noise_gyro = 2.4e-2;        % Gyroscope noise(discrete), rad/s
% %noise_accel = 2.83e-2;      % Accelerometer noise, m/s^2
% 
% noise_gyro = 2.4e-3;        % Gyroscope noise(discrete), rad/s
% noise_accel = 2.83e-2; 
% gravity = 9.81007;          % Gravity magnitude, m/s^2
% 
% persistent states_est P_P_est                % Initial state conditions
% if isempty(states_est)
%     %states_est = [1;0;0;0];             % x_est=[x,y,Vx,Vy,Ax,Ay]'
%     states_est = TX; 
%     P_P_est = 1e-10 * eye(4);
% end
% 
% % ״̬Ԥ�⣬Э����Ԥ��
% states_prd = A * states_est;
% states_prd = states_prd / norm(states_prd);    %��Ԫ����һ��
% G = [-states_est(2)  -states_est(3)  -states_est(4); ...
%           states_est(1)  -states_est(4)   states_est(3); ...
%           states_est(4)   states_est(1)  -states_est(2); ...
%          -states_est(3)   states_est(2)   states_est(1)] / 2;
% %Q = (noise_gyro * dt)^2 * (G * G');
% Q = G *(noise_gyro* dt )^2 * ( G');
% P_P_prd = A * P_P_est * A' + Q;
% 
% %������
% a =[accX; accY; accZ];
% % [ -2*q2,  2*q3, -2*q0, 2*q1]
% % [  2*q1,  2*q0,  2*q3, 2*q2]
% % [  2*q0, -2*q1, -2*q2, 2*q3]
% H =[ -2*states_prd(3), 2*states_prd(4), -2*states_prd(1), 2*states_prd(2);...
%       2*states_prd(2), 2*states_prd(1),  2*states_prd(4), 2*states_prd(3);...
%       2*states_prd(1), -2*states_prd(2), -2*states_prd(3),2*states_prd(4)];
% % H =[ (98*states(3))/5, -(98*states(4))/5,  (98*states(1))/5, -(98*states(2))/5;...
% %     -(98*states(2))/5, -(98*states(1))/5, -(98*states(4))/5, -(98*states(3))/5;...
% %     -(98*states(1))/5,  (98*states(2))/5,  (98*states(3))/5, -(98*states(4))/5];
% innvo = a/norm(a) - H*states_prd;
% %  Measurement noise R.
% R_internal = (noise_accel / norm(a))^2 * eye(3);
% R_external = (1-gravity/norm(a))^2 * eye(3);
% R = R_internal + R_external;
% 
% Kg = P_P_prd * H'/(H*P_P_prd*H'+R);                   %Kg= P(k|k-1)HT/(HP(k|k-1)HT+R)
% states_est = states_prd + Kg* innvo;
% %P_P_est = P_P_prd -P_P_prd * H;    eye(4)
% P_P_est = (eye(4) - Kg * H) * P_P_prd;
% y = states_est;
end 
