clc;
clear;
close all;

%% Systems

Number_of_Spatial_Orbitals = 2;
Number_of_Spin_Up_Electrons = 1;
Number_of_Spin_Down_Electrons = 1;

system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);

system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);

system_plus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons+1,Number_of_Spin_Down_Electrons,true,false);
system_plus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons+1,true,false);

%% Hubbard

U = 1;
t = 1;
connected_ends = true;

hubbard_model = Hubbard(U,t,connected_ends,system,system_minus_up,system_minus_down);

%% Green



