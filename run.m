clc;
clear;
close all;

%% System

Number_of_Spatial_Orbitals = 4;
Number_of_Spin_Up_Electrons = 2; % must be >= 2 to compuate a nonzero spin-up lesser green's function 
Number_of_Spin_Down_Electrons = 2; % must be >= 2 to compute a nonzero spin-down lesser green's function 

system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);

%% Hubbard

U = 1;
t = 1;
connected_ends = true;

system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);

hubbard_model = Hubbard(U,t,connected_ends,system,system_minus_up,system_minus_down);

%% Lesser Green

spin = "up";
spatial_orbital_index_i = 1;
spatial_orbital_index_j = 1;

lesser_green = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model);

t_values = 0:0.1:35;

[lesser_real,lesser_imaginary] = lesser_green.compute(t_values);

%% Greater Green

greater_green = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model);

[greater_real,greater_imaginary] = greater_green.compute(t_values);

%% Plotting

figure;
plot(t_values,greater_imaginary,t_values,greater_real)
title('Greater')

figure;
plot(t_values,lesser_imaginary,t_values,lesser_real)
title('Lesser')

figure;
scatter(greater_green.angular_frequency_differences,greater_green.weights);
hold on;
scatter(lesser_green.angular_frequency_differences,lesser_green.weights);
title('Lesser and Greater')
