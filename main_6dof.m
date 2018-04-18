close all;                          % close all figures

clear;                              % clear all variables

clc; 
%--- Dataset parameters
deltat = 1/200;             % Sampling period

%--- Load data
% Each row of IMU data : 
% (timestamp, wx, wy, wz, ax, ay, az)
% Each row of Ground truth data : 
% (time, position, quaternion, velocity, gyroscope bias, accelerometer bias)
data = load('data/attitude_data.mat');
imu_data = data.imu_data;   % IMU readings
grt_data = data.grt_data;   % Ground truth (GT)
grt_q = grt_data(:, 5:8);   % GT quaternions

bias_w = grt_data(1, 12:14);    % gyroscope bias
bias_a = grt_data(1, 15:17);    % accelerometer bias

%--- Container of the results
N = length(imu_data);
allX = zeros(N, 4);

%--- Initialization
x = grt_q(1,:)';            % Initial state (quaternion)
% P = 1e-10 * eye(4);         % Initial covariance
allX(1,:) = x';
%%
% ahrs = QEAHRSInit(deltat);
% att = zeros(N, 3);
% att(1,:) = [0;0;0];
%%
for k = 2 : N
    % Gyroscope measurements
    w = (imu_data(k-1, 2:4) + imu_data(k, 2:4))/2;
    w = w - bias_w;
    %Acc measurements
    a = imu_data(k, 5:7);
    a = a - bias_a;
     y = ekf_6dof(w(1),w(2),w(3),a(1),a(2),a(3),deltat ,grt_q(1,:)');
%     ahrs = QEAHRSUpdate(ahrs, [w(1);w(2);w(3)], [a(1);a(2);a(3)], [0;0;0], deltat);
%     avp(k,:) = [m2att(ahrs.Cnb); ahrs.kf.xk(5:7); diag(ahrs.kf.Pxk); ahrs.tk];
%       att(k,:)= m2att(ahrs.Cnb);
     allX(k,:) = y'; 
end

%--- Compare the results with ground truth
q_Ws0 = quatinv(grt_q(1,:));
for i=1:N
    grt_q(i,:) = quatmultiply(q_Ws0, grt_q(i,:)); 
    allX(i,:) = quatmultiply(q_Ws0, allX(i,:));
end
[psi, theta, phi] = quat2angle(allX);
[grt_psi, grt_theta, grt_phi] = quat2angle(grt_q);

figure, hold on
plot(1:N, psi,   'r-.', 1:N, grt_psi,   'r');
legend('psi预测值' , 'psi真实值' );
plot(1:N, theta, 'g-.', 1:N, grt_theta, 'g');
legend('theta预测值' , 'theta真实值' );
plot(1:N, phi,   'b-.', 1:N, grt_phi,   'b');
legend('phi预测值' , 'phi真实值' );

% %%
% figure, hold on
% plot(1:N, att(:,1),   'r-.', 1:N, grt_psi,   'r');
% legend('psi预测值' , 'psi真实值' );
% plot(1:N, att(:,2), 'g-.', 1:N, grt_theta, 'g');
% legend('theta预测值' , 'theta真实值' );
% plot(1:N, att(:,3),   'b-.', 1:N, grt_phi,   'b');
% legend('phi预测值' , 'phi真实值' );
