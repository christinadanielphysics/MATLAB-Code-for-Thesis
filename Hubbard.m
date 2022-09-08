classdef Hubbard

    properties
        U
        mu
        t
        connected_ends
        hopping_matrix
        system
        system_minus
        system_plus
    end

    methods
        function obj = Hubbard(U,mu,t,connected_ends,system,system_minus,system_plus)
            obj.U = U;
            obj.mu = mu;
            obj.t = t;
            obj.connected_ends = connected_ends;
            obj.hopping_matrix = obj.get_hopping_matrix();
            obj.system = system;
            obj.system_minus = system_minus;
            obj.system_plus = system_plus;
        end

        function hopping_matrix = get_hopping_matrix(obj)

            hopping_matrix = zeros(obj.system.Number_of_Spatial_Orbitals, obj.system.Number_of_Spatial_Orbitals);
            for i = 1:obj.system.Number_of_Spatial_Orbitals
                for j = 1:obj.system.Number_of_Spatial_Orbitals
                    if abs(i-j)==1
                        hopping_matrix(i,j) = obj.t;
                    end
                    if obj.connected_ends == true
                        if abs(i-j) == obj.system.Number_of_Spatial_Orbitals - 1
                            hopping_matrix(i,j) = obj.t;
                        end
                    end
                end
            end
            return
        end

        function final_states = apply_kinetic_operator(obj,basis_state)
            final_states = OrderedOccupationState.empty;
            counter = 1;
            for i = 1:obj.system.Number_of_Spatial_Orbitals
                for j = 1:obj.system.Number_of_Spatial_Orbitals
                    for spin_index = 1:length(obj.system.Spin_Values)

                        [up_string,down_string] = basis_state.get_strings();
                        key_1 = char("c" + string(j) + string(obj.system.Spin_Values(spin_index)) + " | " + up_string + ";" + down_string + " >");
                        result_1 = obj.system.Annihilation_Map(key_1);

                        [up_string,down_string] = result_1.get_strings();
                        key_2 = char("c†" + string(i) + string(obj.system.Spin_Values(spin_index))  + " | " + up_string + ";" + down_string + " >");
                        result_2 = obj.system_minus.Creation_Map(key_2);

                        result_2.Coefficient = result_2.Coefficient * result_1.Coefficient * obj.hopping_matrix(i,j);

                        final_states(counter) = result_2;
                        counter = counter + 1;
                    end
                end
            end

            return
        end

        function final_states = apply_interaction_operator(obj,basis_state)
            final_states = OrderedOccupationState.empty;
            counter = 1;
            for i = 1:obj.system.Number_of_Spatial_Orbitals

                [up_string,down_string] = basis_state.get_strings();
                key_1 = char("c" + string(i) + "down" + " | " + up_string + ";" + down_string + " >");
                result_1 = obj.system.Annihilation_Map(key_1);

                % 
                
                [up_string,down_string] = result_1.get_strings();
                key_2 = char("c†" + string(i) + "down"  + " | " + up_string + ";" + down_string + " >");
                result_2 = obj.system_minus.Creation_Map(key_2);

                

                [up_string,down_string] = result_2.get_strings();
                key_3 = char("c" + string(i) + "up" + " | " + up_string + ";" + down_string + " >");
                result_3 = obj.system.Annihilation_Map(key_3);

                

                [up_string,down_string] = result_3.get_strings();
                key_4 = char("c†" + string(i) + "up"  + " | " + up_string + ";" + down_string + " >");
                result_4 = obj.system_minus.Creation_Map(key_4);

                result_4.Coefficient = result_4.Coefficient * result_3.Coefficient * result_2.Coefficient * result_1.Coefficient * obj.U;

                final_states(counter) = result_4;
                counter = counter + 1;
                
            end

            return

        end

        function kinetic_matrix = get_kinetic_matrix(obj)
            
            kinetic_matrix = zeros(obj.system.Dimension,obj.system.Dimension);
            for right = 1:obj.system.Dimension
                for left = 1:obj.system.Dimension
                    final_states = obj.apply_kinetic_operator(obj.system.Basis_States(right));
                    kinetic_matrix(left,right) = obj.system.Basis_States(left).scalar_product_one_with_many(final_states);
                end
            end
            return
        end

        function interaction_matrix = get_interaction_matrix(obj)

            interaction_matrix = zeros(obj.system.Dimension,obj.system.Dimension);
            for right = 1:obj.system.Dimension
                for left = 1:obj.system.Dimension
                    final_states = obj.apply_interaction_operator(obj.system.Basis_States(right));
                    interaction_matrix(left,right) = obj.system.Basis_States(left).scalar_product_one_with_many(final_states);
                    if left == right
                        interaction_matrix(left,right) = interaction_matrix(left,right) - obj.mu * obj.system.Number_of_Electrons;
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