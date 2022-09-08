%% Systems

Number_of_Spatial_Orbitals = 3;
Number_of_Spin_Up_Electrons = 2;
Number_of_Spin_Down_Electrons = 2;

system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons);
system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1);

%% Hubbard

U = 1;
mu = U * 0.5;
t = 1;
connected_ends = true;

hubbard_model = Hubbard(U,mu,t,connected_ends,system,system_minus_up,system_minus_down);

