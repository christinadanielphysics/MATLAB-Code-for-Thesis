%% Setup

clc;
clear;
close all;

%% Compressive Sensing Parameters (Part 1/2)

n = 4096; % number of time values = length of signal
p = 128; % number of random samples
perm = round(rand(p,1) * n);

%% Roots, Weights, and Symbolic Form of the Self Energy and the Interacting Retarded Green's Function

%U = 0;
%[d_noninteracting_for_file,n_noninteracting_for_file,noninteracting_denominator_roots,noninteracting_weights,noninteracting_symbolic_inverse,noninteracting_numerator_roots,noninteracting_retarded_weights,noninteracting_retarded_symbolic_inverse] = solve(U,n,perm);
U = 5; 
[d_for_file,n_for_file,interacting_denominator_roots,interacting_weights,interacting_symbolic_inverse,interacting_numerator_roots,interacting_retarded_weights,interacting_retarded_symbolic_inverse] = solve(U,n,perm);

%% Plotting

% Self Energy

% figure; 
% fplot(noninteracting_symbolic_inverse-interacting_symbolic_inverse);
% hold on;
% scatter(noninteracting_denominator_roots,noninteracting_weights,'red','o');
% hold on;
% scatter(interacting_denominator_roots,-interacting_weights,'blue','x');
% title("Self Energy")
% xlabel('$\omega$',interpreter='latex') 
% ylabel('$\Sigma^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')

% Interacting Retarded Green's Function

figure;
fplot(interacting_retarded_symbolic_inverse);
hold on;
scatter(interacting_numerator_roots,interacting_retarded_weights);
title("Interacting Retarded Green's Function");
xlabel('$\omega$',interpreter='latex') 
ylabel('$G^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')

% Noninteracting Retarded Green's Function

% figure;
% fplot(noninteracting_retarded_symbolic_inverse);
% hold on;
% scatter(noninteracting_numerator_roots,noninteracting_retarded_weights);
% title("Noninteracting Retarded Green's Function");
% xlabel('$\omega$',interpreter='latex') 
% ylabel('$G0^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')

%% Computations

function [denominator_polynomial, numerator_polynomial, denominator_roots,weights,symbolic_inverse,numerator_roots,retarded_weights,retarded_symbolic_inverse] = solve(U,n,perm)

    % Compressive Sensing Parameters (Part 2/2)

    Fs = 10; % Sampling frequency
    combine_zero = 1e-8; % computational value of zero for combine procedure
    chop_threshold = 1e-1; % threshold for chopping compressive sensing weights
    
    % Derived Quantities
    
    T = 1/Fs; % Sampling period
    t_values = (0:n-1)*T; % or t_values = (1:n)*T ...
    f = Fs*(0:n-1)/n; % ... with f = Fs*(1:n)/n
    w_values = f*pi;
    delta_f = f(2)-f(1);
    
    % System
    
    Number_of_Spatial_Orbitals = 6;
    Number_of_Spin_Up_Electrons = 3; % must be >= 2 to compuate a spin-up lesser green's function that is nonzero
    Number_of_Spin_Down_Electrons = 3; % must be >= 2 to compute a spin-down lesser green's function that is nonzero
    Number_of_Electrons = Number_of_Spin_Up_Electrons + Number_of_Spin_Down_Electrons;
    system = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,true,true);
    
    % Hubbard
    
    mu = 0.5*U; % chemical potential for half-filling
    t_0 = 0;
    t_1 = 1;
    t_2 = 0;
    connected_ends = true;
    system_minus_up = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons-1,Number_of_Spin_Down_Electrons,false,true);
    system_minus_down = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons-1,false,true);
    hubbard_model = Hubbard(U,t_1,t_0,t_2,connected_ends,system,system_minus_up,system_minus_down);    
    
    % Inverse Retarded Green
    
    my_spin = "up";
    k_value = 0;
    isexact = false; % true for Lehmann frequencies and weights, false for compressive sensing frequencies and weights
    inverse_retarded_green = InverseRetardedGreen(isexact,my_spin,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold,k_value,Number_of_Spatial_Orbitals);
    syms z;

    denominator_polynomial = inverse_retarded_green.form_denominator(z);
    numerator_polynomial = inverse_retarded_green.form_numerator(z);


    D = matlabFunction(denominator_polynomial,'Vars',z);
    N = matlabFunction(numerator_polynomial,'Vars',z);

    
    % Roots of the Inverse Retarded Green's Function
    min = -3;
    max = 3;
    step = 0.1;
    denominator_function_handle = D;
    numerator_function_handle = N;
    denominator_roots = get_roots(min,max,step,denominator_function_handle);
    numerator_roots = get_roots(min,max,step,numerator_function_handle);

    % Weights of the Inverse Retarded Green's Function
    weights = inverse_retarded_green.get_weights_inverse(numerator_roots,denominator_roots);

    % Symbolic Form of the Inverse Retarded Green's Function
    symbolic_numerator = inverse_retarded_green.form_factored_numerator(numerator_roots);
    symbolic_denominator = inverse_retarded_green.form_factored_denominator(denominator_roots);
    symbolic_inverse = symbolic_numerator/symbolic_denominator;

    % Roots of the Retarded Green's Function = numerator_roots

    % Weights of the Retarded Green's Function
    retarded_weights = inverse_retarded_green.get_weights_retarded(numerator_roots,denominator_roots);

    % Symbolic Form of the Retarded Green's Function
    retarded_symbolic_inverse = symbolic_denominator/symbolic_numerator;

end
