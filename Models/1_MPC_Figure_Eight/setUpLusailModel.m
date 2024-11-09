% Model Initialization
% MPC Controller for Lusail International Circuit

%% Add Images to the path (optional if needed)
addpath(genpath('Images'));

%% Load and parse the GeoJSON file for Lusail track
geojsonData = jsondecode(fileread('lusail_track.geojson'));
lusailCoordinates = geojsonData.features(1).geometry.coordinates;

% Extract x and y reference points from coordinates
xRef = lusailCoordinates(:, 1); % Longitudes as x-coordinates
yRef = lusailCoordinates(:, 2); % Latitudes as y-coordinates

% Define reference pose (waypoints)
refPose = [xRef, yRef];

%% Define vehicle parameters used in the models
L = 3; % bicycle model length
ld = 4; % lookahead distance
X_o = xRef(1); % initial vehicle x position
Y_o = yRef(1); % initial vehicle y position 
psi_o = 0; % initialize yaw angle (adjust as necessary)

%% Calculating reference pose vectors for Lusail circuit
% Calculate distance vector
distancematrix = squareform(pdist(refPose));
distancesteps = zeros(length(refPose)-1,1);
for i = 2:length(refPose)
    distancesteps(i-1,1) = distancematrix(i,i-1);
end
totalDistance = sum(distancesteps);
distbp = cumsum([0; distancesteps]);
gradbp = linspace(0, totalDistance, 2500); % Linearize distance

% Interpolate x and y reference points based on distance
xRef2 = interp1(distbp, xRef, gradbp, 'pchip');
yRef2 = interp1(distbp, yRef, gradbp, 'pchip');
yRef2s = smooth(gradbp, yRef2);
xRef2s = smooth(gradbp, xRef2);

% Calculate curvature vector
curvature = getCurvature(xRef2, yRef2);

%% Define vehicle dynamics matrices for MPC Model
% Updated vehicle parameters
tau = 0.5;
Vx = 10;
m = 57; % Updated mass of the vehicle (57 kg)
Iz = 20.52; % Yaw moment of inertia for 57 kg vehicle
Lf = 1.4;
Lr = 1.6;
Cf = 12e3;
Cr = 11e3;

% Longitudinal and lateral models
A1 = [-1/tau 0; 1 0];
B1 = [1/tau; 0];
C1 = [0 1];
D1 = 0;

A2 = [-2*(Cf+Cr)/m/Vx -Vx-2*(Cf*Lf-Cr*Lr)/m/Vx; -2*(Cf*Lf-Cr*Lr)/Iz/Vx -2*(Cf*Lf^2+Cr*Lr^2)/Iz/Vx];
B2 = 2*Cf * [1/m; Lf/Iz];
C2 = [1 0; 0 1];
D2 = [0; 0];

A = [A1 zeros(2,2); zeros(2,2) A2];
B = [B1 zeros(2,1); zeros(2,1) B2];
C = [C1 zeros(1,2); zeros(2,2) C2];
D = [D1 zeros(1,1); zeros(2,1) D2];

%% MPC Pedal Map
rho = 1.21;
Cd = 0.3;
Af = 2;
tire_r = 0.309;

accel_vec = (-4:0.5:4)'; % acceleration range
vel_vec = 0:2:20; % velocity range
torque_map = zeros(length(accel_vec), length(vel_vec));

for i = 1:length(accel_vec)
    for j = 1:length(vel_vec)
        torque_map(i,j) = tire_r * ((m * accel_vec(i)) + (0.5 * rho * Cd * Af * vel_vec(j)^2) + 160);
    end
end

% Convert torque to pedal positions
pedal_map = torque_map;
max_prop_torque = 425 * 9.5;
pedal_map(pedal_map > 0) = pedal_map(pedal_map > 0) / max_prop_torque;
pressure_conv = (0.2 * 7.5e6 * pi * 0.05 * 0.05 * .177 * 2 / 4) * 4 * 1.01;
pedal_map(pedal_map < 0) = pedal_map(pedal_map < 0) / pressure_conv;

%% Curvature Function
function curvature = getCurvature(xRef, yRef)
    DX = gradient(xRef);
    D2X = gradient(DX);
    DY = gradient(yRef);
    D2Y = gradient(DY);
    curvature = (DX .* D2Y - DY .* D2X) ./ (DX.^2 + DY.^2).^(3/2);
end
