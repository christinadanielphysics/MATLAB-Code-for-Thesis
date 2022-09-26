%% Setup

clc;
clear;
close all;

%% Compressive Sensing Parameters

n = 4096; % number of time values = length of signal
p = 128; % number of random samples
Fs = 10; % Sampling frequency

T = 1/Fs; % Sampling period
t_values = (0:n-1)*T; % or t_values = (0:n-1)*T ...
f = Fs*(0:n-1)/n; % ... with f = Fs*(0:n-1)/n
w_values = f*pi;

delta_f = f(2)-f(1);

%% System

Number_of_Spatial_Orbitals = 4;
Number_of_Spin_Up_Electrons = 2; % must be >= 2 to compuate a spin-up lesser green's function that is nonzero
Number_of_Spin_Down_Electrons = 2; % must be >= 2 to compute a spin-down lesser green's function that is nonzero
Number_of_Electrons = Number_of_Spin_Up_Electrons + Number_of_Spin_Down_Electrons;

system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);

%% Hubbard

U = 1;
mu = 0.5*U;
t_1 = 1;
t_0 = 0;
t_2 = 0;
connected_ends = true;

system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);

hubbard_model = Hubbard(U,t_1,t_0,t_2,connected_ends,system,system_minus_up,system_minus_down);

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

compressive_zero = 1e-8;
compressive_threshold = 0.05;

perm = round(rand(p,1) * n);

% greater
y_greater = -greater_imaginary(perm)'; % compressed measurement
greater_compressive = CompressiveSensing(n,y_greater,perm,compressive_zero,compressive_threshold);
s_greater_full = greater_compressive.compute();
[combined_greater_w_values,combined_greater_weights] = greater_compressive.combine(w_values,s_greater_full);
[new_greater_w_values,new_greater_weights] = greater_compressive.chop(combined_greater_w_values,combined_greater_weights);

% lesser
y_lesser = lesser_imaginary(perm)'; % compressed measurement
lesser_compressive = CompressiveSensing(n,y_lesser,perm,compressive_zero,compressive_threshold);
s_lesser_full = lesser_compressive.compute();
[combined_lesser_w_values,combined_lesser_weights] = lesser_compressive.combine(w_values,s_lesser_full);
[new_lesser_w_values,new_lesser_weights] = lesser_compressive.chop(combined_lesser_w_values,combined_lesser_weights);



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
scatter(lesser_green.angular_frequency_differences,lesser_green.weights,'black');
hold on;
scatter(new_greater_w_values,new_greater_weights/abs(sum(new_greater_weights))/2,'blue','x');
hold on;
scatter(-new_lesser_w_values,new_lesser_weights/abs(sum(new_lesser_weights))/2,'blue','x')
title('Lesser and Greater')



