%% Setup

clc;
clear;
close all;

%% Compressive Sensing Parameters (Part 1/2)

n = 4096; % number of time values = length of signal
p = 128; % number of random samples
perm = round(rand(p,1) * n);

%% Roots, Weights, and Symbolic Form of the Self Energy and the Interacting Retarded Green's Function

U = 0;
[d_noninteracting_for_file,n_noninteracting_for_file,noninteracting_denominator_roots,noninteracting_weights,noninteracting_symbolic_inverse,noninteracting_numerator_roots,noninteracting_retarded_weights,noninteracting_retarded_symbolic_inverse] = solve(U,n,perm);
U = 2; 
[d_for_file,n_for_file,interacting_denominator_roots,interacting_weights,interacting_symbolic_inverse,interacting_numerator_roots,interacting_retarded_weights,interacting_retarded_symbolic_inverse] = solve(U,n,perm);

%% Plotting

% % Self Energy
% 
% figure; 
% fplot(noninteracting_symbolic_inverse-interacting_symbolic_inverse);
% hold on;
% scatter(noninteracting_denominator_roots,noninteracting_weights,'red','o');
% hold on;
% scatter(interacting_denominator_roots,-interacting_weights,'blue','x');
% title("Self Energy")
% xlabel('$\omega$',interpreter='latex') 
% ylabel('$\Sigma^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')
% 
% % Interacting Retarded Green's Function
% 
% figure;
% fplot(interacting_retarded_symbolic_inverse);
% hold on;
% scatter(interacting_numerator_roots,interacting_retarded_weights);
% title("Interacting Retarded Green's Function");
% xlabel('$\omega$',interpreter='latex') 
% ylabel('$G^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')
% 
% % Noninteracting Retarded Green's Function
% 
% figure;
% fplot(noninteracting_retarded_symbolic_inverse);
% hold on;
% scatter(noninteracting_numerator_roots,noninteracting_retarded_weights);
% title("Noninteracting Retarded Green's Function");
% xlabel('$\omega$',interpreter='latex') 
% ylabel('$G0^R_{\uparrow}(k=0,\omega)$',Interpreter='latex')

%% Write Data to File

writematrix(noninteracting_denominator_roots,'noninteracting_denominator_roots.dat','Delimiter',' ')  

writematrix(noninteracting_weights,'noninteracting_weights.dat','Delimiter',' ')  

writematrix(interacting_denominator_roots,'interacting_denominator_roots.dat','Delimiter',' ')  

writematrix(interacting_weights,'interacting_weights.dat','Delimiter',' ')  

writematrix(interacting_numerator_roots,'interacting_numerator_roots.dat','Delimiter',' ')  

writematrix(interacting_retarded_weights,'interacting_retarded_weights.dat','Delimiter',' ')  

writematrix(noninteracting_numerator_roots,'noninteracting_numerator_roots.dat','Delimiter',' ')  

writematrix(noninteracting_retarded_weights,'noninteracting_retarded_weights.dat','Delimiter',' ')  

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
    
    Number_of_Spatial_Orbitals = 4;
    Number_of_Spin_Up_Electrons = 2; % must be >= 2 to compuate a spin-up lesser green's function that is nonzero
    Number_of_Spin_Down_Electrons = 2; % must be >= 2 to compute a spin-down lesser green's function that is nonzero
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
    
%     % Lesser and Greater Green
%     
%     spin = "up";
%     spatial_orbital_index_i = 1;
%     spatial_orbital_index_j = 4;
%     lesser_green = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold);
%     [lesser_real,lesser_imaginary] = lesser_green.compute(t_values);
%     greater_green = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard_model,n,perm,t_values,w_values,combine_zero,chop_threshold);
%     [greater_real,greater_imaginary] = greater_green.compute(t_values);
%     
%     %% Compressive Sensing
%     
%     greater_compressive_w_differences = greater_green.approximate_angular_frequency_differences;
%     greater_compressive_weights = greater_green.approximate_weights;
%     lesser_compressive_w_differences = lesser_green.approximate_angular_frequency_differences;
%     lesser_compressive_weights = lesser_green.approximate_weights;
%     
%     %% Plotting
%     
%     figure;
%     plot(t_values,greater_imaginary,'cyan')
%     hold on;
%     plot(t_values,greater_real,'red')
%     hold on;
%     scatter(t_values,lesser_imaginary,'green')
%     hold on;
%     scatter(t_values,lesser_real,'blue')
%     title('Lesser and Greater')
%     figure;
%     scatter(greater_green.angular_frequency_differences,greater_green.weights,'black');
%     hold on;
%     scatter(lesser_green.angular_frequency_differences,lesser_green.weights,'black');
%     hold on;
%     scatter(greater_compressive_w_differences,greater_compressive_weights,'blue','o','MarkerFaceColor', 'b');
%     hold on;
%     scatter(-lesser_compressive_w_differences,lesser_compressive_weights,'blue','o','MarkerFaceColor', 'b');
%     title('Lesser and Greater')
    
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
