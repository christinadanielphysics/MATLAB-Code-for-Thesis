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
        function obj = Hubbard(U,mu,t,connected_ends,hopping_matrix,system,system_minus,system_plus)
            obj.U = U;
            obj.mu = mu;
            obj.t = t;
            obj.connected_ends = connected_ends;
            obj.hopping_matrix = hopping_matrix;
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
            final_states = OderedOccupationState.empty;
            counter = 1;
            for i = 1:obj.system.Number_of_Spatial_Orbitals
                for j = 1:obj.system.Number_of_Spatial_Orbitals
                    for spin_index = 1:length(obj.system.Spin_Values)

                        [up_string,down_string] = basis_state.get_strings();
                        key_1 = char("c" + string(j) + string(obj.system.Spin_Values(spin_index)) + " | " + string(up_string) + ";" + string(down_string) + " >");
                        result_1 = obj.system.Annihilation_Map(key_1);

                        [up_string,down_string] = result_1.get_strings();
                        key_2 = char("c†" + string(i) + string(obj.system.Spin_Values(spin_index))  + " | " + string(up_string) + ";" + string(down_string) + " >");
                        result_2 = obj.system_minus.Creation_Map(key_2);

                        result_2.Coefficient = result_2.Coefficient * result_1.Coefficient * obj.hopping_matrix(i,j);

                        final_states(counter) = result_2;
                        counter = counter + 1;
                    end
                end
            end

            return
        end

    end
end