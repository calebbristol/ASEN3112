%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  For Doing Secret Stuff    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear

%% Part 1: Collected Data

    %% Read in Data
    data = load("Lab 2 - Data for group 2-2");
    Load_case = data(:,1); %[lbs]
    F_0 = data(:,2); %[lbf] %Twin Load Cell
    F_1 = data(:,3); %[lbf] %Twin Load Cell
    F_2 = data(:,4); %[lbf] %Single Load Cell
    F_3D = data(:,5); %[lbf] %Inline Load Cell
    LVDT = data(:,6); %[in] %Linear Variable Differential Transformer

    %% Calculate Linear Models
    lin_0 = polyfit(Load_case,F_0,1);
    lin_1 = polyfit(Load_case,F_1,1);
    lin_2 = polyfit(Load_case,F_2,1);
    lin_3D = polyfit(Load_case,F_3D,1);
    lin_LVDT = polyfit(Load_case,LVDT,1);

    F_0_pred = polyval(lin_0,Load_case);
    F_1_pred = polyval(lin_1,Load_case);
    F_2_pred = polyval(lin_2,Load_case);
    F_3D_pred = polyval(lin_3D,Load_case);
    LVDT_pred = polyval(lin_LVDT,Load_case);

    %% Linear Error Analysis
    resid_0 = F_0 - F_0_pred;
    SS_resid_0 = sum(resid_0.^2);
    SS_total_0 = (length(F_0) - 1) * var(F_0);
    r_sqr_0 = 1 - SS_resid_0 / SS_total_0;
    r_sqr_adj_0 = 1 - (SS_resid_0 / SS_total_0) * (length(F_0) - 1) / (length(F_0) - length(lin_0));

    resid_1 = F_1 - F_1_pred;
    SS_resid_1 = sum(resid_1.^2);
    SS_total_1 = (length(F_1) - 1) * var(F_1);
    r_sqr_1 = 1 - SS_resid_1 / SS_total_1;
    r_sqr_adj_1 = 1 - (SS_resid_1 / SS_total_1) * (length(F_1) - 1) / (length(F_1) - length(lin_1));

    resid_2 = F_2 - F_2_pred;
    SS_resid_2 = sum(resid_2.^2);
    SS_total_2 = (length(F_2) - 1) * var(F_2);
    r_sqr_2 = 1 - SS_resid_2 / SS_total_2;
    r_sqr_adj_2 = 1 - (SS_resid_2 / SS_total_2) * (length(F_2) - 1) / (length(F_2) - length(lin_2));

    resid_3D = F_3D - F_3D_pred;
    SS_resid_3D = sum(resid_3D.^2);
    SS_total_3D = (length(F_3D) - 1) * var(F_3D);
    r_sqr_3D = 1 - SS_resid_3D / SS_total_3D;
    r_sqr_adj_3D = 1 - (SS_resid_3D / SS_total_3D) * (length(F_3D) - 1) / (length(F_3D) - length(lin_3D));

    resid_LVDT = LVDT - LVDT_pred;
    SS_resid_LVDT = sum(resid_LVDT.^2);
    SS_total_LVDT = (length(LVDT) - 1) * var(LVDT);
    r_sqr_LVDT = 1 - SS_resid_LVDT / SS_total_LVDT;
    r_sqr_adj_LVDT = 1 - (SS_resid_LVDT / SS_total_LVDT) * (length(LVDT) - 1) / (length(LVDT) - length(lin_LVDT));

    %% Plotting
    figure()
    plot(Load_case,F_0,'linewidth',2); hold on
    plot(Load_case,F_0_pred,'linewidth',2)
    title("Load Cell 0 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 0 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off

    figure()
    plot(Load_case,F_1,'linewidth',2); hold on
    plot(Load_case,F_1_pred,'linewidth',2)
    title("Load Cell 1 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 1 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off

    figure()
    plot(Load_case,F_2,'linewidth',2); hold on
    plot(Load_case,F_2_pred,'linewidth',2)
    title("Load Cell 2 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 2 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off

    figure()
    plot(Load_case,F_3D,'linewidth',2); hold on
    plot(Load_case,F_3D_pred,'linewidth',2)
    title("Inline Load Cell vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Inline Internal Force [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off

    figure()
    plot(Load_case,LVDT,'linewidth',2); hold on
    plot(Load_case,LVDT_pred,'linewidth',2)
    title("LVDT vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Displacement [in]")
    legend("Data","Linear Approximation",'location','NE')
    grid on
    hold off

    figure()
    subplot(2,2,1)
    plot(Load_case,F_0,'linewidth',2); hold on
    plot(Load_case,F_0_pred,'linewidth',2)
    title("Load Cell 0 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 0 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    subplot(2,2,2)
    plot(Load_case,F_1,'linewidth',2); hold on
    plot(Load_case,F_1_pred,'linewidth',2)
    title("Load Cell 1 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 1 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off
    subplot(2,2,3)
    plot(Load_case,F_2,'linewidth',2); hold on
    plot(Load_case,F_2_pred,'linewidth',2)
    title("Load Cell 2 vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Reaction Force 2 [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off
    subplot(2,2,4)
    plot(Load_case,F_3D,'linewidth',2); hold on
    plot(Load_case,F_3D_pred,'linewidth',2)
    title("Inline Load Cell vs. Central Loading Case")
    xlabel("Load Case [lbs]")
    ylabel("Inline Internal Force [lbf]")
    legend("Data","Linear Approximation",'location','NW')
    grid on
    hold off
    hold off
    
%% Part 2: ANSYS Analytical Model

    %% Read In Data
    ANSYS_dis = load("Nodal_Solution_data");
    Node_dis = ANSYS_dis(:,1);
    u_x = ANSYS_dis(:,2);
    u_y = ANSYS_dis(:,3);
    u_z = ANSYS_dis(:,4);
    u_sum = ANSYS_dis(:,5);
    
    ANSYS_load = load("Reaction_Solution_data");
    Node_load = ANSYS_load(:,1);
    F_x = ANSYS_load(:,2);
    F_y = ANSYS_load(:,3);
    F_z = ANSYS_load(:,4);
    
    ANSYS_int = load("Bar Forces - Luis Alvidrez");
    Element_int = ANSYS_int(:,1);
    F_int = ANSYS_int(:,2);
    
    % Convert internal force to stress
    r_e = 0.375 / 2; %[in]
    t = 1 / 16; %[in]
    r_e = r_e * 0.0254; %[m]
    t = t * 0.0254; %[m]
    r_i = r_e - t;
    A = pi * (r_e^2 - r_i^2);
    Element_stress = F_int ./ A;
    
    %% Plotting
    %
    %
    
    figure()
    plot(Element_int,F_int); hold on
    xlabel('Element Number')
    ylabel('Internal Force [N]')
    title('Internal Bar Forces from ANSYS Base Case')
    xlim([Element_int(1) Element_int(end)])
    grid on
    hold off
    
    figure()
    plot(Element_int,Element_stress); hold on
    xlabel('Element Number')
    ylabel('Internal Stress [Pa]')
    title('Internal Bar Stress from ANSYS Base Case')
    xlim([Element_int(1) Element_int(end)])
    grid on
    hold off
    
    
    %% Factor of Safety Analysis
    sigma_max = 310e6;
    sigma = max(Element_stress);
    
    FS = sigma_max/sigma;
    
    
    %% Error Analysis (Different Cases)
    A_69 = load("InternalForces69Gp");
    
    Element_69 = A_69(:,1);
    F_69 = A_69(:,2);
    
    
    figure()
    plot(Element_int,F_int); hold on
    plot(Element_69,F_69)
    xlabel('Element Number')
    ylabel('Internal Force [N]')
    title('Internal Bar Forces from ANSYS Base Case')
    xlim([Element_int(1) Element_int(end)])
    grid on
    hold off

    A_69dot5 = load("InternalForces695Gp");
    Element_69dot5 = A_69dot5(:,1);
    F_69dot5 = A_69dot5(:,2);
    
        
    figure()
    plot(Element_69dot5,F_69dot5); hold on
    plot(Element_69,F_69)
    xlabel('Element Number')
    ylabel('Internal Force [N]')
    title('Internal Bar Forces from ANSYS Base Case')
    xlim([Element_69dot5(1) Element_69dot5(end)])
    grid on
    hold off
    
    % Beams Case
    A_beams_f = load("reaction_forces_beams");
    A_beams_d = load("Beam_disp");
    
    F_x_beam = A_beams_f(:,2);
    F_y_beam = A_beams_f(:,3);
    F_z_beam = A_beams_f(:,4);
    
    u_x_beam = A_beams_d(:,2);
    u_y_beam = A_beams_d(:,3);
    u_z_beam = A_beams_d(:,4);
    
    error_f_beam = norm([std(F_x_beam - F_x) std(F_y_beam - F_y) std(F_z_beam - F_z)]);
    error_d_beam = norm([std(u_x_beam - u_x) std(u_y_beam - u_y) std(u_z_beam - u_z)]);
