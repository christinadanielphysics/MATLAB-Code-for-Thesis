%% Setup

clc;
clear;
close all;

%% Compressive Sensing Parameters

start_time_inclusive = 0;
end_time_inclusive = 50;
number_of_time_values = 500;

%% System

Number_of_Spatial_Orbitals = 4;
Number_of_Spin_Up_Electrons = 2; % must be >= 2 to compuate a spin-up lesser green's function that is nonzero
Number_of_Spin_Down_Electrons = 2; % must be >= 2 to compute a spin-down lesser green's function that is nonzero

system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);

%% Hubbard

U = 1.5;
t = 1.5;
connected_ends = true;

system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);

hubbard_model = Hubbard(U,t,connected_ends,system,system_minus_up,system_minus_down);

%% Lesser Green

spin = "up";
spatial_orbital_index_i = 1;
spatial_orbital_index_j = 1;

lesser_green = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model);

t_values = linspace(start_time_inclusive,end_time_inclusive,number_of_time_values);

[lesser_real,lesser_imaginary] = lesser_green.compute(t_values);

%% Greater Green

greater_green = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model);

[greater_real,greater_imaginary] = greater_green.compute(t_values);

%% Plotting

figure;
plot(t_values,greater_imaginary,'cyan')
hold on;
plot(t_values,greater_real,'red')
hold on;
scatter(t_values,lesser_imaginary,'green')
hold on;
scatter(t_values,lesser_real,'blue')
title('Lesser and Greater')

figure;
scatter(greater_green.angular_frequency_differences,greater_green.weights);
hold on;
scatter(lesser_green.angular_frequency_differences,lesser_green.weights);
title('Lesser and Greater')
