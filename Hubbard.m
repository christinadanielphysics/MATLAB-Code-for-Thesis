classdef Hubbard

    properties
        U
        t
        t_0
        t_2
        connected_ends
        system
        system_minus_up
        system_minus_down
        number_of_spatial_orbitals
        system_dimension
        annihilation_map
        creation_map_for_up
        creation_map_for_down
        spin_values
        number_of_electrons
        basis_states
        eigenvectors
        eigenvalues
        number_of_spin_up_electrons
        number_of_spin_down_electrons
    end

    methods
        function obj = Hubbard(U,t_1,t_0,t_2,connected_ends,system,system_minus_up,system_minus_down)
            obj.U = U;
            obj.t = t_1;
            obj.t_0 = t_0;
            obj.t_2 = t_2;
            obj.connected_ends = connected_ends;
            obj.system = system;
            obj.system_minus_up = system_minus_up;
            obj.system_minus_down = system_minus_down;
            obj.number_of_spatial_orbitals = system.Number_of_Spatial_Orbitals;
            obj.system_dimension = system.Dimension;
            obj.annihilation_map = system.Annihilation_Map;
            obj.creation_map_for_up = system_minus_up.Creation_Map;
            obj.creation_map_for_down = system_minus_down.Creation_Map;
            obj.spin_values = system.Spin_Values;
            obj.number_of_electrons = system.Number_of_Electrons;
            obj.basis_states = system.Basis_States;
            [eigenvectors,eigenvalues] = obj.diagonalize_hamiltonian_matrix();
            obj.eigenvectors = eigenvectors;
            obj.eigenvalues = diag(eigenvalues);
            obj.number_of_spin_up_electrons = system.Number_of_Spin_Up_Electrons;
            obj.number_of_spin_down_electrons = system.Number_of_Spin_Down_Electrons;
        end

        function hopping_matrix = get_hopping_matrix(obj)

            hopping_matrix = zeros(obj.number_of_spatial_orbitals, obj.number_of_spatial_orbitals);
            for i = 1:obj.number_of_spatial_orbitals
                for j = 1:obj.number_of_spatial_orbitals
                    if abs(i-j)==1
                        hopping_matrix(i,j) = obj.t;
                    end
                    if abs(i-j)==2
                        hopping_matrix(i,j) = obj.t_2;
                    end
                    if obj.connected_ends == true
                        if abs(i-j) == obj.number_of_spatial_orbitals - 1
                            hopping_matrix(i,j) = obj.t;
                        end
                        if abs(i-j) == obj.number_of_spatial_orbitals - 2
                            hopping_matrix(i,j) = obj.t_2;
                        end
                    end
                    if i == j
                        hopping_matrix(i,j) = obj.t_0;
                    end
                end
            end
            return
        end

        function final_states = apply_kinetic_operator(obj,basis_state)

            [up_string_1,down_string_1] = basis_state.get_strings();

            hopping_matrix = obj.get_hopping_matrix();

            final_states = OrderedOccupationState.empty;
            counter = 1;
            for i = 1:obj.number_of_spatial_orbitals
                for j = 1:obj.number_of_spatial_orbitals
                    for spin_index = 1:length(obj.spin_values)

                        key_1 = char("c" + string(j) + obj.spin_values(spin_index) + " | " + up_string_1 + ";" + down_string_1 + " >");
                        result_1 = obj.annihilation_map(key_1);

                        [up_string_2,down_string_2] = result_1.get_strings();
                        key_2 = char("c???" + string(i) + obj.spin_values(spin_index)  + " | " + up_string_2 + ";" + down_string_2 + " >");

                        if obj.spin_values(spin_index) == "up"
                            result_2 = obj.creation_map_for_up(key_2);
                            result_2.Coefficient = result_2.Coefficient * result_1.Coefficient * hopping_matrix(i,j);
                            final_states(counter) = result_2;
                            counter = counter + 1;
                        else
                            result_2 = obj.creation_map_for_down(key_2);
                            result_2.Coefficient = result_2.Coefficient * result_1.Coefficient * hopping_matrix(i,j);
                            final_states(counter) = result_2;
                            counter = counter + 1;
                        end
                    end
                end
            end

            return
        end

        function final_states = apply_interaction_operator(obj,basis_state)
            final_states = OrderedOccupationState.empty;
            counter = 1;
            for i = 1:obj.number_of_spatial_orbitals

                [up_string,down_string] = basis_state.get_strings();
                key_1 = char("c" + string(i) + "down" + " | " + up_string + ";" + down_string + " >");
                result_1 = obj.annihilation_map(key_1);

                
                [up_string,down_string] = result_1.get_strings();
                key_2 = char("c???" + string(i) + "down"  + " | " + up_string + ";" + down_string + " >");
                result_2 = obj.creation_map_for_down(key_2);

                

                [up_string,down_string] = result_2.get_strings();
                key_3 = char("c" + string(i) + "up" + " | " + up_string + ";" + down_string + " >");
                result_3 = obj.annihilation_map(key_3);

                

                [up_string,down_string] = result_3.get_strings();
                key_4 = char("c???" + string(i) + "up"  + " | " + up_string + ";" + down_string + " >");
                result_4 = obj.creation_map_for_up(key_4);

                result_4.Coefficient = result_4.Coefficient * result_3.Coefficient * result_2.Coefficient * result_1.Coefficient * obj.U;

                final_states(counter) = result_4;
                counter = counter + 1;
                
            end

            return

        end

        function kinetic_matrix = get_kinetic_matrix(obj)
            
            kinetic_matrix = zeros(obj.system_dimension,obj.system_dimension);
            for right = 1:obj.system_dimension
                final_states = obj.apply_kinetic_operator(obj.basis_states(right));
                for left = 1:obj.system_dimension
                    kinetic_matrix(left,right) = obj.basis_states(left).scalar_product_one_with_many(final_states);
                end
            end
            return
        end

        function interaction_matrix = get_interaction_matrix(obj)

            interaction_matrix = zeros(obj.system_dimension,obj.system_dimension);
            for right = 1:obj.system_dimension
                final_states = obj.apply_interaction_operator(obj.basis_states(right));
                for left = 1:obj.system_dimension
                    interaction_matrix(left,right) = obj.basis_states(left).scalar_product_one_with_many(final_states);
                    if left == right
                        interaction_matrix(left,right)  = interaction_matrix(left,right) - 0.5 * obj.U * obj.number_of_electrons;
                    end
                end
            end
            return
        end

        function hamiltonian_matrix = get_hamiltonian_matrix(obj)
            hamiltonian_matrix = obj.get_kinetic_matrix() + obj.get_interaction_matrix();
            return
        end

        function [eigenvectors,eigenvalues] = diagonalize_hamiltonian_matrix(obj)
            [eigenvectors,eigenvalues] = eig( obj.get_hamiltonian_matrix() );
            return
        end
    end
end