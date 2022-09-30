%% Setup
clc;
clear;
close all;
%% Compressive Sensing Parameters
n = 4096; % number of time values = length of signal
p = 128; % number of random samples
Fs = 10; % Sampling frequency
combine_zero = 1e-8; % computational value of zero for combine procedure
chop_threshold = 1e-1; % threshold for chopping compressive sensing weights
%% Derived Quantities
T = 1/Fs; % Sampling period
t_values = (0:n-1)*T; % or t_values = (1:n)*T ...
f = Fs*(0:n-1)/n; % ... with f = Fs*(1:n)/n
w_values = f*pi;
delta_f = f(2)-f(1);
perm = round(rand(p,1) * n);
%% System
Number_of_Spatial_Orbitals = 4;
Number_of_Spin_Up_Electrons = 2; % must be >= 2 to compuate a spin-up lesser green's function that is nonzero
Number_of_Spin_Down_Electrons = 2; % must be >= 2 to compute a spin-down lesser green's function that is nonzero
Number_of_Electrons = Number_of_Spin_Up_Electrons + Number_of_Spin_Down_Electrons;
system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);
%% Hubbard
U = 0;
mu = 0.5*U;
t_0 = 0;
t_1 = 1;
t_2 = 0;
connected_ends = true;
system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);
hubbard_model = Hubbard(U,t_1,t_0,t_2,connected_ends,system,system_minus_up,system_minus_down);
%% Lesser and Greater Green

% spin = "up";
% spatial_orbital_index_i = 1;
% spatial_orbital_index_j = 4;
% lesser_green = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold);
% [lesser_real,lesser_imaginary] = lesser_green.compute(t_values);
% greater_green = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold);
% [greater_real,greater_imaginary] = greater_green.compute(t_values);
% 
% %% Compressive Sensing
% 
% greater_compressive_w_differences = greater_green.approximate_angular_frequency_differences;
% greater_compressive_weights = greater_green.approximate_weights;
% lesser_compressive_w_differences = lesser_green.approximate_angular_frequency_differences;
% lesser_compressive_weights = lesser_green.approximate_weights;
% 
% %% Plotting
% 
% figure;
% plot(t_values,greater_imaginary,'cyan')
% hold on;
% plot(t_values,greater_real,'red')
% hold on;
% scatter(t_values,lesser_imaginary,'green')
% hold on;
% scatter(t_values,lesser_real,'blue')
% title('Lesser and Greater')
% figure;
% scatter(greater_green.angular_frequency_differences,greater_green.weights,'black');
% hold on;
% scatter(lesser_green.angular_frequency_differences,lesser_green.weights,'black');
% hold on;
% scatter(greater_compressive_w_differences,greater_compressive_weights,'blue','o','MarkerFaceColor', 'b');
% hold on;
% scatter(-lesser_compressive_w_differences,lesser_compressive_weights,'blue','o','MarkerFaceColor', 'b');
% title('Lesser and Greater')

%% Inverse Retarded Green

my_spin = "up";
k_value = pi;
isexact = false; % true for Lehmann frequencies and weights, false for compressive sensing frequencies and weights
inverse_retarded_green = InverseRetardedGreen(isexact,my_spin,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold,k_value,Number_of_Spatial_Orbitals);
syms z;
inverse_retarded_green.form_denominator(z)
inverse_retarded_green.form_numerator(z)

%% Roots

denominator_roots = [];
for x0 = -3:0.1:3
    fun = @g;
    root = fzero(fun,x0);
    if ~ismember(round(root,4),round(denominator_roots,4))
        denominator_roots(end+1) = root;
    end
end

denominator_roots


numerator_roots = [];
for x0 = -3:0.1:3
    fun = @h;
    root = fzero(fun,x0);
    if ~ismember(round(root,4),round(numerator_roots,4))
        numerator_roots(end+1) = root;
    end
end

numerator_roots


%% Weights of the Self Energy

% Form numerator just based on factors from roots


% Form denominator just based on factors from roots


% Remove a given root r my multiplying N/D by (z-r)


% Set z equal to the value of the root


% Use the subs command to evaluate the weight for that root for N/D


% Do the same for the inverse noninteracting green's function


% Combine the results for the inverse of the noninteracting and the inverse of the interacting green's
% functions, to get the roots and weights of the self energy


% Form the imaginary and real parts of the self energy

%% Weights of the Interacting Retarded Green's Function


% Analyze the roots and weights of the interacting retarded green's
% function


% Form the imaginary and real parts of the interacting retarded green's function


%% Ground State Energy

% derivation

% computation

