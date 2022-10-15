%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%  ASEN 3112 Lab 1               %
%  Author: Caleb Bristol         %
%  Date: 09/16/21                %
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminary
clc
clear

%% Read in Data

T_20_data = importdata("Torsion20inlbf_146.txt");
T_400_data = importdata("Torsion400inlbf_145.txt");
T_20 = T_20_data.data;
T_400 = T_400_data.data;
%Column 1: Time [s]
%Column 2: Angle [degrees]
    % Convert Angle to Radians
    phi_400 = T_400(:,2) * (pi/180);
    phi_20 = T_20(:,2) * (pi/180);
%Column 3: Shear Strain [degrees]
%Column 4: Torque [inlbf]
%Column 5: Axial stress[in]


%% Constants
D_e = 0.75; %[in]
R_e = D_e / 2; %[in]
t = 1/16; %[in]
R_i = R_e - t; %[in]
L_ext = 1; %[in]
L = 10; %[in]
G = 3.75e06; %[psi]

J = (pi / 2) * (R_e^4 - R_i^4);


%% Convert Angle to Max Shear Strain (Prediction from Angle)

%T = (phi * G * J) / L
%tau_max = (T * R) / J = (phi * G * J * R) / (L * J) = (phi * G * R) / L
%gamma_max = (tau_max) / G = (phi * R) / L

gamma_max_400 = (T_400(:,2) - T_400(1,2)) .* (pi/180) .* (R_e / L_ext);



    %% Plot

    figure(1)

    plot(T_400(:,4),T_400(:,3)); hold on
    plot(T_400(:,4),gamma_max_400)
    xlabel("Torque [inlbf]")
    ylabel("Strain [degrees]")
    title("Strain vs. Torque")
    legend("Measured","Predicted")
    hold off
    
    
%% Calculating GJ for Closed Thin Wall

% (T * L) = (G * J) * phi
%
% This forms a linear pattern where TL is a function of phi
%
% The slope of this linear function is GJ

% Redefine TL Vector By Multiplying Torque by Scalar 
TL_400 = T_400(:,4) * L;

% Polyfit to find slope of function
GJ_400 = polyfit(phi_400,TL_400,1);

Gamma_400 = deg2rad(T_400(:,3));

GJ_400_Gamma = polyfit(Gamma_400,T_400(:,4),1) * R_e;
