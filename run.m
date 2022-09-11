%% Setup

clc;
clear;
close all;

%% Compressive Sensing Parameters

n = 4096; % number of time values = length of signal
p = 128;

Fs = 10; % Sampling frequency
T = 1/Fs; % Sampling period
t_values = (0:n-1)*T;
f = Fs*(0:n-1)/n;
w_values = f*pi;

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

[lesser_real,lesser_imaginary] = lesser_green.compute(t_values);

%% Greater Green

greater_green = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model);

[greater_real,greater_imaginary] = greater_green.compute(t_values);



%% Compressive Sensing
% Theta*s == y from Steve Brunton is the same as A*x = b from Emmanuel Candes

% greater
perm = round(rand(p,1) * n);
y = -greater_imaginary(perm)'; % compressed measurement

Psi = dct(eye(n,n)); % build Psi
Theta = Psi(perm,:); % measure rows of Psi

% L1-Minimization using CVX
cvx_begin;
    variable s(n);
    minimize( norm(s,1) );
    subject to
        Theta*s == y;
cvx_end;

% lesser
y_lesser = lesser_imaginary(perm)'; % compressed measurement

% L1-Minimization using CVX
cvx_begin;
    variable s_lesser(n);
    minimize( norm(s_lesser,1) );
    subject to
        Theta*s_lesser == y_lesser;
cvx_end;

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
scatter(greater_green.angular_frequency_differences,greater_green.weights,'red');
hold on;
scatter(lesser_green.angular_frequency_differences,lesser_green.weights,'red');
hold on;
plot(w_values,s,'blue');
hold on;
plot(-w_values,s_lesser,'blue')
title('Lesser and Greater')



