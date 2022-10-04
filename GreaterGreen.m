classdef GreaterGreen

    properties
        spin
        spatial_orbital_index_i
        spatial_orbital_index_j
        site_difference
        hubbard
        eigenvalues_for_hubbard
        eigenvectors_for_hubbard
        hubbard_plus
        eigenvalues_for_hubbard_plus
        eigenvectors_for_hubbard_plus
        system
        system_minus_up
        system_minus_down
        system_for_element
        system_plus_up_for_element
        system_plus_down_for_element
        b_max
        weights % each weight is a product of two matrix elements
        angular_frequency_differences
        approximate_weights % compressive sensing
        approximate_angular_frequency_differences % compressive sensing
    end

    methods
        function obj = GreaterGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold)
            if nargin > 0
                obj.site_difference = abs(spatial_orbital_index_i - spatial_orbital_index_j);

                obj.spin = spin; % operator (i.e. "up" or "down")
                obj.spatial_orbital_index_i = spatial_orbital_index_i;
                obj.spatial_orbital_index_j = spatial_orbital_index_j;
                obj.hubbard = hubbard;
    
                obj.eigenvalues_for_hubbard = hubbard.eigenvalues;
                obj.eigenvectors_for_hubbard = hubbard.eigenvectors;
    
                if obj.spin == "up"
                    obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons+1,hubbard.number_of_spin_down_electrons,true,false);
                    obj.system_minus_up = hubbard.system;
                    obj.system_minus_down = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons+1,hubbard.number_of_spin_down_electrons-1,false,true);
                else % "down"
                    obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons+1,true,false);
                    obj.system_minus_up = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons+1,false,true);
                    obj.system_minus_down = hubbard.system;
                end
    
                hubbard_plus = Hubbard(hubbard.U,hubbard.t,hubbard.t_0,hubbard.t_2,hubbard.connected_ends,obj.system,obj.system_minus_up,obj.system_minus_down);
                obj.hubbard_plus = hubbard_plus;
                obj.eigenvalues_for_hubbard_plus = hubbard_plus.eigenvalues;
                obj.eigenvectors_for_hubbard_plus = hubbard_plus.eigenvectors;
                
                obj.system_for_element = hubbard.system;
                obj.system_plus_up_for_element = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons+1,hubbard.number_of_spin_down_electrons,true,false);
                obj.system_plus_down_for_element = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons+1,true,false);
    
                obj.b_max = length(obj.eigenvalues_for_hubbard_plus);
                weights = zeros(1,obj.b_max);
                angular_frequency_differences = zeros(1,obj.b_max);
                for b = 1:obj.b_max
                    matrix_element_1 = MatrixElement("annihilation",spatial_orbital_index_i,spin,1,hubbard,b,hubbard_plus,obj.system_for_element,obj.system_for_element,obj.system_for_element,obj.system_plus_up_for_element,obj.system_plus_down_for_element,"greater");
                    matrix_element_2 = MatrixElement("creation",spatial_orbital_index_j,spin,b,hubbard_plus,1,hubbard,obj.system_for_element,obj.system_for_element,obj.system_for_element,obj.system_plus_up_for_element,obj.system_plus_down_for_element,"greater");
                    weights(1,b) = matrix_element_1.compute() * matrix_element_2.compute();
                    angular_frequency_differences(1,b) = obj.eigenvalues_for_hubbard_plus(b) - obj.eigenvalues_for_hubbard(1);
                end
                obj.weights = weights;
                obj.angular_frequency_differences = angular_frequency_differences;

                % Write Data
                prefix = "greater"+"_"+string(spatial_orbital_index_i)+"_"+string(spatial_orbital_index_j)+"_";
                writematrix(angular_frequency_differences,prefix+"angular_frequency_differences.dat",'Delimiter',' ')  
                writematrix(weights,prefix+"weights.dat",'Delimiter',' ') 
    
                [~,greater_imaginary] = obj.compute(t_values);
                y_greater = -greater_imaginary(perm)'; % compressed measurement
                greater_compressive = CompressiveSensing(n,y_greater,perm,combine_zero,chop_threshold);
                s_greater_full = greater_compressive.compute();
                [combined_w_values,combined_weights] = greater_compressive.combine(w_values,s_greater_full);
                [compressive_w_differences,compressive_weights] = greater_compressive.chop(combined_w_values,combined_weights);
                compressive_weights = compressive_weights/sum(abs(compressive_weights))/2; % normalize the weights
    
                obj.approximate_angular_frequency_differences = compressive_w_differences;
                obj.approximate_weights = compressive_weights;
            end
        end

        function [data_real,data_imaginary] = compute(obj,t_values)

            data_real = zeros(1,length(t_values));
            data_imaginary = zeros(1,length(t_values));

            for t_index = 1:length(t_values)
                data_complex = 0;
                for b = 1:obj.b_max
                    w_difference = obj.eigenvalues_for_hubbard(1) - obj.eigenvalues_for_hubbard_plus(b);
                    data_complex = data_complex + (-1i) * obj.weights(b) * exp( (1i)*w_difference*t_values(t_index) );
                end
                data_real(1,t_index) = real(data_complex);
                data_imaginary(1,t_index) = imag(data_complex);
            end

            return
        end
    end
end