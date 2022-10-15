%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3112 Lab 4
% Author: Caleb Bristol
% Date: 12/07/21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear 
close all;

%% Read In Data
cubic = load('cubicdata.mat');
rectang = load('rectangdata.mat');
cubic_2020 = load('Squarerod_FA2020.mat');
rectang_2020 = load('ThinBar_FA2020.mat');

cubicdata = cubic.data;
rectangdata = rectang.data;
cubicdata_2020 = cubic_2020.Squarebar_FA2020;
rectangdata_2020 = rectang_2020.Thinbar_FA2020;

%% Convert Voltage into Force
cubicdata(:,4) = cubicdata(:,2) ./ 2.37e-3;
rectangdata(:,4) = rectangdata(:,2) ./ 2.37e-3;
cubicdata_2020(:,4) = cubicdata_2020(:,2) ./ 2.37e-3;
rectangdata_2020(:,4) = rectangdata_2020(:,2) ./ 2.37e-3;

%% Plotting

    %{
    %% Square Rod
    figure()
    plot(cubicdata(:,3),cubicdata(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Square Rod')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    xline(0.07,'--r','LineWidth',2)
    yline(305,'--g','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    grid on
    hold off

    
    %% Thin Bar
    figure()
    plot(rectangdata(:,3),rectangdata(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Thin Bar')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    xline(0.1,'--r','LineWidth',2)
    yline(160,'--g','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    grid on
    hold off
    %}

    
    %% Question 1
    
    % Square Rod
    figure()
    plot(cubicdata_2020(:,3),cubicdata_2020(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Square Rod')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    yline(265.55,'--r','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    grid on
    hold off

    
    % Thin Bar
    figure()
    plot(rectangdata_2020(:,3),rectangdata_2020(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Thin Bar')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    yline(135.24,'--r','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    ylim([0 140])
    grid on
    hold off
    
    
    %% Question 2
    
    % Square Rod
    figure()
    plot(cubicdata_2020(:,3),cubicdata_2020(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Square Rod')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    xline(0.67,'--r','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    grid on
    hold off

    
    % Thin Bar
    figure()
    plot(rectangdata_2020(:,3),rectangdata_2020(:,4),'LineWidth',2); hold on
    title('Applied Load vs. Vertical Deflection for Thin Bar')
    xlabel('Lateral Deflection, \delta [in]')
    ylabel('Load, P [lb]')
    set(gca,'FontSize',14)
    xline(0.32,'--r','LineWidth',2)
    legend('Deflection Curve','Approximate Post Buckling')
    ylim([0 140])
    grid on
    hold off
    
%% Problem 3

    %% Function for P_cr
    P_cr = @(E,I,L) pi^2 .* E .* I ./ L.^2;
    
    %% Define Values
    E = 1e7;
    I_square = 3.1e-4;
    I_rect = 1.628e-4;
    L_real = linspace(10,12,1000);
    L = linspace(3,15,1000);
    L_eff = L/2;
    
    %% Calculate buckling Loads
    P_cr_sq_real = P_cr(E,I_square,L_real);
    P_cr_re_real = P_cr(E,I_rect,L_real);
    P_cr_sq = P_cr(E,I_square,L);
    P_cr_re = P_cr(E,I_rect,L);
    P_cr_sq_eff = P_cr(E,I_square,L_eff);
    P_cr_re_eff = P_cr(E,I_rect,L_eff);
    
    %% Calculate Yielding Loads
    A_square = 0.25^2 - (0.25-0.0625)^2;
    A_rect = 0.125 * 1;
    sigma_max = 35000;
    
    F_yield_sq = sigma_max * A_square;
    F_yield_re = sigma_max * A_rect;
    
    %% Plotting
    
    % Realistic Values (L = 10:12)
    figure()
    plot(L_real,P_cr_sq_real,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Square Rod')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(265.55,'--r','LineWidth',2)
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load')
    grid on
    hold off
    
    figure()
    plot(L_real,P_cr_re_real,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Thin Bar')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(135.24,'--r','LineWidth',2)
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load')
    grid on
    hold off
    
    
    % Full Range Values (L = 3:15)
    
    figure()
    plot(L,P_cr_sq,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Square Rod, Simply Supported')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(265.55,'--r','LineWidth',2)
    yline(F_yield_sq,'--g','LineWidth',2)
    xlim([L(1) L(end)])
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load','Specimen Yielding Load')
    grid on
    hold off
    
    figure()
    plot(L,P_cr_re,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Thin Bar, Simply Supported')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(135.24,'--r','LineWidth',2)
    yline(F_yield_sq,'--g','LineWidth',2)
    xlim([L(1) L(end)])
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load','Specimen Yielding Load')
    grid on
    hold off
    
    
    % Fixed Fixed Load (part b)
    
    figure()
    plot(L,P_cr_sq_eff,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Square Rod, Fixed-Fixed')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(265.55,'--r','LineWidth',2)
    yline(F_yield_sq,'--g','LineWidth',2)
    xlim([L(1) L(end)])
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load','Specimen Yielding Load')
    grid on
    hold off
    
    figure()
    plot(L,P_cr_re_eff,'LineWidth',2); hold on
    title('Buckling Load vs. Specimen Length: Thin Bar, Fixed-Fixed')
    xlabel('Specimen Length [in]')
    ylabel('Buckling Load [lb]')
    yline(135.24,'--r','LineWidth',2)
    yline(F_yield_sq,'--g','LineWidth',2)
    xlim([L(1) L(end)])
    set(gca,'FontSize',14)
    legend('Predicted Load','Specimen Predicted Buckling Load','Specimen Yielding Load')
    grid on
    hold off
    
    

