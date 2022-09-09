classdef System

    properties
        Number_of_Spatial_Orbitals
        Number_of_Spin_Up_Electrons 
        Number_of_Spin_Down_Electrons
        Number_of_Electrons
        Annihilation_Map
        Creation_Map
        Spin_Values
        Basis_States
        Dimension
    end

    methods
        function obj = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons,make_annihilation_map,make_creation_map)
            obj.Number_of_Spatial_Orbitals = Number_of_Spatial_Orbitals;
            obj.Number_of_Spin_Up_Electrons = Number_of_Spin_Up_Electrons;
            obj.Number_of_Spin_Down_Electrons = Number_of_Spin_Down_Electrons;
            if make_annihilation_map == true
                obj.Creation_Map = obj.get_creation_map();
            end
            if make_creation_map == true
                obj.Annihilation_Map = obj.get_annihilation_map();
            end
            obj.Spin_Values = ["up","down"];
            obj.Basis_States = obj.get_basis_states();
            obj.Dimension = length(obj.Basis_States);
            obj.Number_of_Electrons = Number_of_Spin_Up_Electrons + Number_of_Spin_Down_Electrons;
        end

        function basis_states = get_basis_states(obj)
            up_array_matrix = nchoosek(1:obj.Number_of_Spatial_Orbitals,obj.Number_of_Spin_Up_Electrons);
            down_array_matrix = nchoosek(1:obj.Number_of_Spatial_Orbitals,obj.Number_of_Spin_Down_Electrons);

            [numRows_for_up,numCols_for_up] = size(up_array_matrix);
            [numRows_for_down,numCols_for_down] = size(down_array_matrix);

            basis_states = OrderedOccupationState.empty;

            counter = 1;

            for up_row = 1:numRows_for_up

                up_array = zeros(numCols_for_up,1);
                for up_col = 1:numCols_for_up
                    up_array(up_col,1) = up_array_matrix(up_row,up_col);
                end

                for down_row = 1:numRows_for_down

                    down_array = zeros(numCols_for_down,1);
                    for down_col = 1:numCols_for_down
                        down_array(down_col,1) = down_array_matrix(down_row,down_col);
                    end

                    state = OrderedOccupationState(1,up_array,down_array);
                    basis_states(1,counter) = state;
                    counter = counter + 1;
                end
            end

            return

        end

        function annihilation_operators = get_annihilation_operators(obj)

            annihilation_operators = AnnihilationOperator.empty;
            spin_values = ["up","down"];
            counter = 1;
            for number = 1:obj.Number_of_Spatial_Orbitals
                for spin_index = 1:length(spin_values)
                    annihilation_operators(1,counter) = AnnihilationOperator(number,spin_values(spin_index));
                    counter = counter + 1;
                end
            end

            return
        end

        function creation_operators = get_creation_operators(obj)

            creation_operators = CreationOperator.empty;
            spin_values = ["up","down"];
            counter = 1;
            for number = 1:obj.Number_of_Spatial_Orbitals 
                for spin_index = 1:length(spin_values)
                    creation_operators(1,counter) = CreationOperator(number,spin_values(spin_index));
                    counter = counter + 1;
                end
            end

            return
   
        end

        function annihilation_map = get_annihilation_map(obj)
            
            annihilation_operators = obj.get_annihilation_operators();

            keySet = {};
            basis_states = obj.get_basis_states();
            valueSet = {};

            map_counter = 1;
            for operator_index = 1:length(annihilation_operators)
                annihilation_operator = annihilation_operators(operator_index);
                for counter = 1:length(basis_states)
                    initial_state = basis_states(1,counter);
                    [up_string,down_string] = initial_state.get_strings();
                    key_for_initial_state = char("c" + string(annihilation_operator.Spatial_Orbital_Index) + annihilation_operator.Spin + " | " + up_string + ";" + down_string + " >");
                    keySet{map_counter} = key_for_initial_state;
                    final_state = annihilation_operator.apply(initial_state);
                    valueSet{map_counter} = final_state;
                    map_counter = map_counter + 1;
                end
            end
            
            annihilation_map = containers.Map(keySet,valueSet);
            return
        end

        function creation_map = get_creation_map(obj)

            creation_operators = obj.get_creation_operators();
            
            keySet = {};
            basis_states = obj.get_basis_states();
            valueSet = {};

            map_counter = 1;
            for operator_index = 1:length(creation_operators)
                creation_operator = creation_operators(operator_index);
                for counter = 1:length(basis_states)
                    initial_state = basis_states(1,counter);
                    [up_string,down_string] = initial_state.get_strings();
                    key_for_initial_state = char("c†" + string(creation_operator.Spatial_Orbital_Index) + creation_operator.Spin  + " | " + up_string + ";" + down_string + " >");
                    keySet{map_counter} = key_for_initial_state;
                    final_state = creation_operator.apply(initial_state);
                    valueSet{map_counter} = final_state;
                    map_counter = map_counter + 1;
                end
            end
            
            creation_map = containers.Map(keySet,valueSet);
            return
        end
    end
end

