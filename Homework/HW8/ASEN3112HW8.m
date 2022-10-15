%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 Homework 8
% Author: Caleb Bristol
% Date: 11/02/21
%
%%%%%%%%%%%%%%%%%%%%%%%%

%% 
clc
clear
close all;

%% Establish Constants
%
% All variables are unitless in this problem
E = 3000;
L = 30;
P = 200;
A_1 = 2;
A_2 = 4;
A_3 = 3;
L_1 = L;
L_2 = L * (2 * sqrt(3) / 3);
L_3 = L * (2 / 3);
phi_1 = 0;
phi_2 = 30;
phi_3 = 120;

%% Functions to Help Out
T_block = @(phi) [cosd(phi) sind(phi);-sind(phi) cosd(phi)];
T = @(phi) [T_block(phi) zeros(2);zeros(2) T_block(phi)];
block1 = @(phi) [cosd(phi)^2 sind(phi)*cosd(phi);sind(phi)*cosd(phi) sind(phi)^2];
K_mat = @(E,A,L,phi) (E * A / L) * [block1(phi) -block1(phi);-block1(phi) block1(phi)];

%% Step 1: Globalization
K_1 = K_mat(E,A_1,L,phi_1); %Contains Nodes 1 and 4
K_2 = K_mat(E,A_2,L_1,phi_2); %Nodes 2 and 4
K_3 = K_mat(E,A_3,L_3,phi_3); %Nodes 3 and 4


%% Step 2: Merge
z = zeros(2);
K_1_master = [K_1(1:2,1:2) z z K_1(1:2,3:4);z z z z;z z z z;K_1(3:4,1:2) z z K_1(3:4,3:4)];
K_2_master = [z z z z;z K_2(1:2,1:2) z K_2(1:2,3:4);z z z z;z K_2(3:4,1:2) z K_2(3:4,3:4)];
K_3_master = [z z z z;z z z z;[z z;z z] K_3];

K_master = K_1_master + K_2_master + K_3_master;

%% Step 3: Boundary Conditions
K_reduced = K_master(7:8,7:8);
f_reduced = [0;-P];


%% Step 4: Displacement Solution
u_reduced = linsolve(K_reduced,f_reduced);

%% Step 5: Recovery of Reactions
u = [0;0;0;0;0;0;u_reduced];
f = K_master * u;

u_1 = u([1,2,7,8]);
u_2 = u([3,4,7,8]);
u_3 = u([5,6,7,8]);


%% Step 6: Recovery of Internal Forces
u_bar_1 = T(phi_1) * u_1;
u_bar_2 = T(phi_2) * u_2;
u_bar_3 = T(phi_3) * u_3;

d_1 = u_bar_1(3) - u_bar_1(1);
d_2 = u_bar_2(3) - u_bar_2(1);
d_3 = u_bar_3(3) - u_bar_3(1);

F_1 = (E * A_1 / L_1) * d_1;
F_2 = (E * A_2 / L_2) * d_2;
F_3 = (E * A_3 / L_3) * d_3;


%% Display Results to Command Window
fprintf('K1 \n')
disp(K_1)
fprintf('K2 \n')
disp(K_2)
fprintf('K3 \n')
disp(K_3)
fprintf('Master K Matrix \n')
disp(K_master)
fprintf('Reduced K Matrix \n')
disp(K_reduced)
fprintf('Reduced f vector \n')
disp(f_reduced)
fprintf('Complete u vector \n')
disp(u)
fprintf('Complete f vector \n')
disp(f)
fprintf('Internal Force in Bar 1 \n')
disp(F_1)
fprintf('Internal Force in Bar 2 \n')
disp(F_2)
fprintf('Internal Force in Bar 3 \n')
disp(F_3)
