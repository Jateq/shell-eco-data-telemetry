%% Add Images to the path
addpath(genpath('Images'));

%% Load the scene data files
% Load data from Driving Scenario Designer
load('USHighway.mat');
refPose = data.ActorSpecifications(1,46).Waypoints;

%% Load the new track data
% Load CSV data from the new track file
track_data = readmatrix('Shifted_UTM_Coordinates__Relative_to_Origin_.csv');


%% Define data for velocity lookup table
lookUpt = readmatrix('velocityDistributionHighway.xlsx');
xlt = lookUpt(2:42,1);
ylt = lookUpt(1,2:31);
vel = lookUpt(2:42,2:31)*4/5;

%% Specify simulation stop time
Ts = 1000*5/4;

%% Define vehicle parameters used in the models
L = 3; % bicycle length
ld = 4; % lookahead distance
X_o = 0; % initial vehicle position
Y_o = 0; % initial vehicle position 
psi_o = 88.5*(pi/180); % initialize yaw angle

%% Calculate reference pose vectors
% Calculate distance vector based on new track data
distancematrix = squareform(pdist([xRef, yRef]));
distancesteps = zeros(length(xRef)-1,1);
for i = 2:length(xRef)
    distancesteps(i-1,1) = distancematrix(i,i-1);
end
totalDistance = sum(distancesteps); % Total traveled distance
distbp = cumsum([0; distancesteps]); % Distance for each waypoint
gradbp = linspace(0, totalDistance, 100); % Linearize distance

% Interpolate X and Y vectors based on distance
xRef2 = interp1(distbp, xRef, gradbp, 'pchip');
yRef2 = interp1(distbp, yRef, gradbp, 'pchip');
yRef2s = smooth(gradbp, yRef2);
xRef2s = smooth(gradbp, xRef2);

% Calculate curvature vector
curvature = getCurvature(xRef2, yRef2);

%% Define reference time for plotting 
s = size(xRef);
tRef = linspace(0, Ts, length(gradbp)); % Used in 2D visualization for plotting reference points

%% Calculate A, B, C matrices used in MPC Model
% Vehicle parameters
tau = 0.5;
Vx = 10;
m = 57;
Iz = 128;
Lf = 1.4;
Lr = 1.6;
Cf = 12e3;
Cr = 11e3;

% Longitudinal model
A1 = [-1/tau 0; 1 0];
B1 = [1/tau; 0];
C1 = [0 1];
D1 = 0;

% Lateral model
A2 = [-2*(Cf+Cr)/m/Vx -Vx-2*(Cf*Lf-Cr*Lr)/m/Vx;...
      -2*(Cf*Lf-Cr*Lr)/Iz/Vx -2*(Cf*Lf^2+Cr*Lr^2)/Iz/Vx];
B2 = 2*Cf*[1/m; Lf/Iz];
C2 = [1 0; 0 1];
D2 = [0; 0];

% Combined model
A = [A1 zeros(2,2); zeros(2,2) A2];
B = [B1 zeros(2,1); zeros(2,1) B2];
C = [C1 zeros(1,2); zeros(2,2) C2];
D = [D1 zeros(1,1); zeros(2,1) D2];

%% MPC Pedal Map
% Additional vehicle parameters
rho = 1.21;
Cd = 0.3;
Af = 2;
tire_r = 0.309;

% Bounds for 2-D lookup table
accel_vec = (-4:0.5:4)'; % Acceleration range
vel_vec = 0:2:20; % Velocity range
torque_map = zeros(length(accel_vec), length(vel_vec));

% Calculate required torque
for i = 1:length(accel_vec)
    for j = 1:length(vel_vec)
        % Calculate torque based on tractive and resistive forces
        torque_map(i,j) = tire_r*((m*accel_vec(i))+(0.5*rho*Cd*Af*vel_vec(j)^2)+160);
    end
end

% Convert torque to pedal based on powertrain parameters
pedal_map = torque_map;
max_prop_torque = 425*9.5;
pedal_map(pedal_map>0) = pedal_map(pedal_map>0)/max_prop_torque;

% Calculate the conversion from torque to maximum pressure
pressure_conv = (0.2*7.5e6*pi*0.05*0.05*.177*2/4)*4*1.01; % Calibrated value
pedal_map(pedal_map<0) = pedal_map(pedal_map<0)/pressure_conv;

sim('mpcHighway.slx')

%% Curvature Function
function curvature = getCurvature(xRef, yRef)
    % Calculate gradient by the gradient of the X and Y vectors
    DX = gradient(xRef);
    D2X = gradient(DX);
    DY = gradient(yRef);
    D2Y = gradient(DY);
    curvature = (DX.*D2Y - DY.*D2X) ./(DX.^2 + DY.^2).^(3/2);
end 
