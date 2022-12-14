classdef MatrixElement

    properties
        type % creation or annihilation
        spatial_orbital_index
        spin
        left_subscript
        left_hubbard
        basis_left
        right_subscript
        right_hubbard
        basis_right
        y_max
        x_max
        left_eigenvectors
        right_eigenvectors
        left_number_of_electrons
        right_number_of_electrons
        system
        system_minus_up
        system_minus_down
        system_plus_up
        system_plus_down
        lesser_or_greater
        map
    end

    methods
        function obj = MatrixElement(type,spatial_orbital_index,spin,left_subscript,left_hubbard,right_subscript,right_hubbard,system,system_minus_up,system_minus_down,system_plus_up,system_plus_down,lesser_or_greater)
            obj.type = type; % operator (i.e. "creation" or "annihilation")
            obj.spatial_orbital_index = spatial_orbital_index; % operator
            obj.spin = spin; % operator
            obj.left_subscript = left_subscript; % subscript for left eigenstate (e.g. subscript for ground state is 1)
            obj.left_hubbard = left_hubbard; % hubbard object for left eigenstate
            obj.right_subscript = right_subscript; % subscript for right eigenstate (e.g. subscript for ground state is 1)
            obj.right_hubbard = right_hubbard; % hubbard object for right eigenstate
            obj.basis_left = left_hubbard.basis_states; % occupation basis for left eigenstate
            obj.basis_right = right_hubbard.basis_states; % occupation basis for right eigenstate
            obj.y_max = length(obj.basis_left);
            obj.x_max = length(obj.basis_right);
            obj.left_eigenvectors = left_hubbard.eigenvectors;
            obj.right_eigenvectors = right_hubbard.eigenvectors;
            obj.left_number_of_electrons = left_hubbard.number_of_electrons;
            obj.right_number_of_electrons = right_hubbard.number_of_electrons;
            
            obj.system = system; % N electrons
            obj.system_minus_up = system_minus_up; % one up electron missing
            obj.system_minus_down = system_minus_down; % one down electron missing
            obj.system_plus_up = system_plus_up; % one extra up electron
            obj.system_plus_down = system_plus_down; % one extra down electron

            obj.lesser_or_greater = lesser_or_greater; % "lesser" or "greater" 

            % obtain the correct map for the matrix element
            if obj.type == "creation"
                if obj.lesser_or_greater == "lesser"
                    if obj.spin == "up"
                        obj.map = system_minus_up.Creation_Map;
                    else
                        obj.map = system_minus_down.Creation_Map;
                    end
                else % "greater"
                    obj.map = system.Creation_Map;
                end

            else % "annihilation"
                if obj.lesser_or_greater == "lesser"
                    obj.map = system.Annihilation_Map;
                else % "greater"
                    if obj.spin == "up"
                        obj.map = system_plus_up.Annihilation_Map;
                    else
                        obj.map = system_plus_down.Annihilation_Map;
                    end
                end
            end


        end

        function matrix_element = compute(obj)
            
            left_eigenvector = obj.left_eigenvectors(:,obj.left_subscript);
            right_eigenvector = obj.right_eigenvectors(:,obj.right_subscript);

            operator_string = "";
            if obj.type == "creation"
                operator_string = "c???";
            else
                operator_string = "c";
            end

            matrix_element = 0;
            for x = 1:obj.x_max
                bracket_x = right_eigenvector(x);
                for y = 1:obj.y_max
                    bracket_y = conj(left_eigenvector(y));
                        
                    % get right basis state
                    right_basis_state = obj.basis_right(x);

                    % get key
                    [up_string,down_string] = right_basis_state.get_strings();
                    key = char(operator_string + string(obj.spatial_orbital_index) + obj.spin + " | " + up_string + ";" + down_string + " >");
                    
                    resulting_state = obj.map(key);
                    
                    % get left basis state
                    left_basis_state = obj.basis_left(y);

                    inner_matrix_element = left_basis_state.scalar_product(resulting_state);

                    matrix_element = matrix_element + bracket_y*bracket_x*inner_matrix_element;
                    
                end
            end

            return
        end
    end
end