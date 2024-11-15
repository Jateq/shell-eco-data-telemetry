% Model Initialization
% ADT MPC Controller for Lusail Circuit Based on US Highway Model

%% Add Images to the path
addpath(genpath('Images'));

%% Load and parse the UTM coordinates for Lusail track
% Load the UTM coordinates from the new CSV file
track_data = readtable('Shifted_UTM_Coordinates__Relative_to_Origin_.csv', 'VariableNamingRule', 'preserve');
track_data.Properties.VariableNames = {'Easting', 'Northing'};

% Define main track centerline from UTM coordinates
xRef = track_data.Easting;
yRef = track_data.Northing;

% Define Track Width and Boundaries
track_width = 12;  % Total width in meters
half_width = track_width / 2;

% Calculate track boundaries based on the centerline
dx = gradient(xRef);
dy = gradient(yRef);
norm_factor = sqrt(dx.^2 + dy.^2);  % Normalize for perpendicular direction

% Left and Right boundaries (offset by half of the track width)
xLeft = xRef + half_width * (-dy ./ norm_factor);
yLeft = yRef + half_width * (dx ./ norm_factor);
xRight = xRef - half_width * (-dy ./ norm_factor);
yRight = yRef - half_width * (dx ./ norm_factor);

%% Optional: Plot the track with boundaries for visualization





%% Original Vehicle Parameters and Simulation Settings (from US Highway)
%% Define data for velocity lookup table
lookUpt = readmatrix('velocityDistributionHighway.xlsx');
xlt = lookUpt(2:42,1);
ylt = lookUpt(1,2:31);
vel = lookUpt(2:42,2:31)*4/5;

%% Specify simulation stop time
Ts = 1000*5/4;

% Define vehicle parameters used in the models
L = 10; % Bicycle length
ld = 100; % Lookahead distance
X_o = xRef(1); % Initial vehicle x position
Y_o = yRef(1); % Initial vehicle y position 
dx_initial = xRef(2) - xRef(1);
dy_initial = yRef(2) - yRef(1);
psi_o = atan2(dy_initial, dx_initial);
vehicle_x_position = [];
vehicle_y_position = [];


%% Simulation Parameters
dt = 0.1; % Time step
Ts = 100; % Total simulation time
num_steps = floor(Ts / dt);
vehicle_x_position = zeros(num_steps, 1);
vehicle_y_position = zeros(num_steps, 1);
vehicle_yaw = zeros(num_steps, 1);

% Initialize vehicle state
vehicle_x_position(1) = X_o;
vehicle_y_position(1) = Y_o;
vehicle_yaw(1) = psi_o;

% Simulation loop
for k = 2:num_steps
    % Compute the desired position from the reference path
    % You may need to interpolate if your time steps don't match
    idx = min(k, length(xRef2s));
    desired_x = xRef2s(idx);
    desired_y = yRef2s(idx);

    % Compute control inputs (steering angle, acceleration)
    % For example, using a simple proportional controller
    error_x = desired_x - vehicle_x_position(k-1);
    error_y = desired_y - vehicle_y_position(k-1);
    steering_angle = atan2(error_y, error_x) - vehicle_yaw(k-1);
    acceleration = 0; % Assuming constant speed for simplicity

    % Update vehicle state
    vehicle_speed = Vx; % Assuming constant speed
    vehicle_x_position(k) = vehicle_x_position(k-1) + vehicle_speed * cos(vehicle_yaw(k-1)) * dt;
    vehicle_y_position(k) = vehicle_y_position(k-1) + vehicle_speed * sin(vehicle_yaw(k-1)) * dt;
    vehicle_yaw(k) = vehicle_yaw(k-1) + (vehicle_speed / L) * tan(steering_angle) * dt;
end

% ... [plotting code] ...

figure;
plot(xRef, yRef, 'k-', 'LineWidth', 1.5); % Main track
hold on;
plot(xLeft, yLeft, 'r--', 'LineWidth', 1); % Left boundary
plot(xRight, yRight, 'b--', 'LineWidth', 1); % Right boundary
plot(X_o, Y_o, 'go', 'MarkerFaceColor', 'g'); % Starting point
plot(vehicle_x_position, vehicle_y_position, 'go-', 'LineWidth', 1); % Vehicle path
title('Vehicle Path Along Lusail Circuit Track');
xlabel('Easting (m)');
ylabel('Northing (m)');
legend('Center Line', 'Left Boundary', 'Right Boundary', 'Starting Position', 'Vehicle Path');
grid on;
axis equal;
hold off;


plot(vehicle_x_position, vehicle_y_position, 'go-', 'LineWidth', 1); % Vehicle path


% Calculate distance vector for the Lusail track (replaces the US Highway data)
refPose = [xRef, yRef];
distancematrix = squareform(pdist(refPose));
distancesteps = zeros(length(refPose)-1,1);
for i = 2:length(refPose)
    distancesteps(i-1,1) = distancematrix(i,i-1);
end
totalDistance = sum(distancesteps); % Total traveled distance
distbp = cumsum([0; distancesteps]);
gradbp = linspace(0, totalDistance, 100); % Linearize distance

% Interpolate x and y reference points based on distance
xRef2 = interp1(distbp, xRef, gradbp, 'pchip');
yRef2 = interp1(distbp, yRef, gradbp, 'pchip');
xRef2s = smooth(gradbp, xRef2);
yRef2s = smooth(gradbp, yRef2);

% Calculate curvature vector
curvature = getCurvature(xRef2, yRef2);

%% Define Reference Time for Plotting
Ts = 45 * 5/4;
tRef = linspace(0, Ts, length(gradbp));

% Original MPC Model Parameters (same as US Highway)

tau = 0.5;
Vx = 10;
m = 57; % Vehicle mass (updated)
Iz = 128; % Yaw moment of inertia for updated mass
Lf = 1.4;
Lr = 1.6;
Cf = 12e3;
Cr = 11e3;

% Longitudinal and lateral model matrices for MPC
A1 = [-1/tau 0; 1 0];
B1 = [1/tau; 0];
C1 = [0 1];
D1 = 0;

A2 = [-2*(Cf+Cr)/m/Vx -Vx-2*(Cf*Lf-Cr*Lr)/m/Vx; -2*(Cf*Lf-Cr*Lr)/Iz/Vx -2*(Cf*Lf^2+Cr*Lr^2)/Iz/Vx];
B2 = 2*Cf * [1/m; Lf/Iz];
C2 = [1 0; 0 1];
D2 = [0; 0];

% Combined state-space model for MPC
A = [A1 zeros(2,2); zeros(2,2) A2];
B = [B1 zeros(2,1); zeros(2,1) B2];
C = [C1 zeros(1,2); zeros(2,2) C2];
D = [D1 zeros(1,1); zeros(2,1) D2];

%% MPC Pedal Map (Optional)
rho = 1.21;
Cd = 0.3;
Af = 2;
tire_r = 0.309;

% Define acceleration and velocity ranges for torque calculation
accel_vec = (-4:0.5:4)'; % Acceleration range
vel_vec = 0:2:20; % Velocity range
torque_map = zeros(length(accel_vec), length(vel_vec));

for i = 1:length(accel_vec)
    for j = 1:length(vel_vec)
        torque_map(i,j) = tire_r * ((m * accel_vec(i)) + (0.5 * rho * Cd * Af * vel_vec(j)^2) + 160);
    end
end

% Convert torque to pedal map
pedal_map = torque_map;
max_prop_torque = 425 * 9.5;
pedal_map(pedal_map > 0) = pedal_map(pedal_map > 0) / max_prop_torque;
pressure_conv = (0.2 * 7.5e6 * pi * 0.05 * 0.05 * .177 * 2 / 4) * 4 * 1.01;
pedal_map(pedal_map < 0) = pedal_map(pedal_map < 0) / pressure_conv;

%% Curvature Calculation Function
function curvature = getCurvature(xRef, yRef)
    DX = gradient(xRef);
    D2X = gradient(DX);
    DY = gradient(yRef);
    D2Y = gradient(DY);
    curvature = (DX .* D2Y - DY .* D2X) ./ (DX.^2 + DY.^2).^(3/2);
end


