classdef LesserGreen

    properties
        spin
        spatial_orbital_index_i
        spatial_orbital_index_j
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
        weights % a product of two matrix elements
    end

    methods
        function obj = LesserGreen(spin,spatial_orbital_index_i,spatial_orbital_index_j,hubbard)
            obj.spin = spin; % operator (i.e. "up" or "down")
            obj.spatial_orbital_index_i = spatial_orbital_index_i;
            obj.spatial_orbital_index_j = spatial_orbital_index_j;
            obj.hubbard = hubbard;

            obj.eigenvalues_for_hubbard = hubbard.eigenvalues;
            obj.eigenvectors_for_hubbard = hubbard.eigenvectors;

            if obj.spin == "up"
                obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons,true,true);
                obj.system_minus_up = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-2,hubbard.number_of_spin_down_electrons,false,true);
                obj.system_minus_down = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons-1,false,true);
            else % "down"
                obj.system = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons-1,true,true);
                obj.system_minus_up = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons-1,hubbard.number_of_spin_down_electrons-1,false,true);
                obj.system_minus_down = System(hubbard.number_of_spatial_orbitals,hubbard.number_of_spin_up_electrons,hubbard.number_of_spin_down_electrons-2,false,true);
            end

            hubbard_minus = Hubbard(hubbard.U,hubbard.t,hubbard.connected_ends,obj.system,obj.system_minus_up,obj.system_minus_down);
            obj.hubbard_minus = hubbard_minus;
            obj.eigenvalues_for_hubbard_minus = hubbard_minus.eigenvalues;
            obj.eigenvectors_for_hubbard_minus = hubbard_minus.eigenvectors;

            
            obj.system_for_element = hubbard.system;
            obj.system_minus_up_for_element = hubbard.system_minus_up;
            obj.system_minus_down_for_element = hubbard.system_minus_down;
            
            obj.a_max = length(obj.eigenvalues_for_hubbard_minus);
            weights = zeros(1,obj.a_max); % array to store each product of two matrix elements; one product for each index a
            for a = 1:obj.a_max
                matrix_element_1 = MatrixElement("creation",spatial_orbital_index_j,spin,1,hubbard,a,hubbard_minus,obj.system_for_element,obj.system_minus_up_for_element,obj.system_minus_down_for_element,obj.system_for_element,obj.system_for_element,"lesser");
                matrix_element_2 = MatrixElement("annihilation",spatial_orbital_index_i,spin,a,hubbard_minus,1,hubbard,obj.system_for_element,obj.system_minus_up_for_element,obj.system_minus_down_for_element,obj.system_for_element,obj.system_for_element,"lesser");
                weights(1,a) = matrix_element_1.compute() * matrix_element_2.compute();
            end
            obj.weights = weights;
        end

        function [data_real,data_imaginary] = compute(obj,t_values)

            data_real = zeros(1,length(t_values));
            data_imaginary = zeros(1,length(t_values));

            for t_index = 1:length(t_values)
                data_complex = 0;
                for a = 1:obj.a_max
                   w_difference = obj.eigenvalues_for_hubbard_minus(a) - obj.eigenvalues_for_hubbard(1);
                   data_complex = data_complex + (1i) * obj.weights(a) * exp( (1i)*(w_difference)*t_values(t_index) );
                end
                data_real(1,t_index) = real(data_complex);
                data_imaginary(1,t_index) = imag(data_complex);
            end
            
            return
        end

        function [angular_frequencies,weights] = compute_angular_frequencies_and_weights(obj)
            angular_frequencies = [];
            weights = [];
            % code goes here
            return
        end
    end
end

