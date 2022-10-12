classdef LesserGreen

    properties
        spin
        spatial_orbital_index_i
        spatial_orbital_index_j
        site_difference
        hubbard
        eigenvalues_for_hubbard
        eigenvectors_for_hubbard
        hubbard_minus
        eigenvalues_for_hubbard_minus
        eigenvectors_for_hubbard_minus
        system
        system_minus_up
        system_minus_down
        system_for_element
        system_minus_up_for_element
        system_minus_down_for_element
        a_max
        weights % each weight is a product of two matrix elements
        angular_frequency_differences
        approximate_weights % compressive sensing
        approximate_angular_frequency_differences % compressive sensing
    end

    methods
        function obj = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold)
            if nargin > 0
                obj.site_difference = abs(spatial_orbital_index_i - spatial_orbital_index_j);

                obj.spin = spin; % operator (i.e. "up" or "down")
                obj.spatial_orbital_index_i = spatial_orbital_index_i;
                obj.spatial_orbital_index_j = spatial_orbital_index_j;
                obj.hubbard = hubbard;
    
                obj.eigenvalues_for_hubbard = hubbard.eigenvalues;
                obj.eigenvectors_for_hubbard = hubbard.eigenvectors;
    
                if obj.spin == "up"
                    obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons,true,false);
                    obj.system_minus_up = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-2,hubbard.number_of_spin_down_electrons,false,true);
                    obj.system_minus_down = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons-1,false,true);
                else % "down"
                    obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons-1,true,false);
                    obj.system_minus_up = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons-1,false,true);
                    obj.system_minus_down = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons-2,false,true);
                end
    
                hubbard_minus = Hubbard(hubbard.U,hubbard.t,hubbard.t_0,hubbard.t_2,hubbard.connected_ends,obj.system,obj.system_minus_up,obj.system_minus_down);
                obj.hubbard_minus = hubbard_minus;
                obj.eigenvalues_for_hubbard_minus = hubbard_minus.eigenvalues;
                obj.eigenvectors_for_hubbard_minus = hubbard_minus.eigenvectors;
                
                obj.system_for_element = hubbard.system;
                obj.system_minus_up_for_element = hubbard.system_minus_up;
                obj.system_minus_down_for_element = hubbard.system_minus_down;
                
                obj.a_max = length(obj.eigenvalues_for_hubbard_minus);
                weights = zeros(1,obj.a_max); % array to store each product of two matrix elements; one product for each index a
                angular_frequency_differences = zeros(1,obj.a_max); % array to store each difference of two eigenvalues; one difference for each index a
                for a = 1:obj.a_max
                    matrix_element_1 = MatrixElement("creation",spatial_orbital_index_j,spin,1,hubbard,a,hubbard_minus,obj.system_for_element,obj.system_minus_up_for_element,obj.system_minus_down_for_element,obj.system_for_element,obj.system_for_element,"lesser");
                    matrix_element_2 = MatrixElement("annihilation",spatial_orbital_index_i,spin,a,hubbard_minus,1,hubbard,obj.system_for_element,obj.system_minus_up_for_element,obj.system_minus_down_for_element,obj.system_for_element,obj.system_for_element,"lesser");
                    weights(1,a) = matrix_element_1.compute() * matrix_element_2.compute();
                    angular_frequency_differences(1,a) = obj.eigenvalues_for_hubbard(1) - obj.eigenvalues_for_hubbard_minus(a); 
                end
                obj.weights = weights;
                obj.angular_frequency_differences = angular_frequency_differences;

                % Write Data
                folder = "./U="+string(hubbard.U)+"/";
                prefix = folder + "lesser"+"_"+string(spatial_orbital_index_i)+"_"+string(spatial_orbital_index_j)+"_";
                writematrix(angular_frequency_differences,prefix+"angular_frequency_differences.dat",'Delimiter',' ')       
                writematrix(weights,prefix+"weights.dat",'Delimiter',' ') 

%                 obj.weights = importdata(prefix+"weights.dat");
%                 obj.angular_frequency_differences = readmatrix(prefix+"angular_frequency_differences.dat");


                [~,lesser_imaginary] = obj.compute(t_values);  


                y_lesser = lesser_imaginary(perm)'; % compressed measurement
                lesser_compressive = CompressiveSensing(n,y_lesser,perm,combine_zero,chop_threshold);
                s_lesser_full = lesser_compressive.compute();
                [combined_w_values,combined_weights] = lesser_compressive.combine(w_values,s_lesser_full);
                [compressive_w_differences,compressive_weights] = lesser_compressive.chop(combined_w_values,combined_weights);
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
                for a = 1:obj.a_max
                   w_difference = obj.eigenvalues_for_hubbard_minus(a) - obj.eigenvalues_for_hubbard(1);
                   data_complex = data_complex + (1i) * obj.weights(a) * exp( (1i)*w_difference*t_values(t_index) );
                end
                data_real(1,t_index) = real(data_complex);
                data_imaginary(1,t_index) = imag(data_complex);
            end
            
            return
        end

        function result = numerical_derivative(obj,t_1,t_2)
            two_t_values = [t_1,t_2];
            [lesser_real,lesser_imaginary] = obj.compute(two_t_values);
            
            delta_t = t_2 - t_1;
            result_real = ( lesser_real(2) - lesser_real(1) ) / delta_t;
            result_imaginary = ( lesser_imaginary(2) - lesser_imaginary(1) ) / delta_t;
            result = result_real + result_imaginary;
            return
        
        end
    end
end

