clc;
clear;
close all;

%% System

Number_of_Spatial_Orbitals = 4;
Number_of_Spin_Up_Electrons = 2;
Number_of_Spin_Down_Electrons = 2;

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

t_values = 0:0.1:10;

[data_real,data_imaginary] = lesser_green.compute(t_values);
plot(data_real)


%% Greater Green

system_plus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons+1,Number_of_Spin_Down_Electrons,true,false);
system_plus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons+1,true,false);