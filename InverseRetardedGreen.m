classdef InverseRetardedGreen
    properties
        isexact % true for Lehmann, false for compressive sensing
        perm % for consistency, use the same indices for all applications of compressive sensing 
        k_value % value that can be negative, zero, positive
        V % number of sites = number of spatial orbitals
        spin % spin for fermionic operators
        hubbard % for N electrons
        n
        t_values
        w_values
        combine_zero
        chop_threshold
        lesser_green_matrix
        greater_green_matrix
    end

    methods
        function obj = InverseRetardedGreen(isexact,spin,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold,k_value,V)
            
            lesser_green_matrix(V,V) = LesserGreen(spin,1,1,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold);
            greater_green_matrix(V,V) = GreaterGreen(spin,1,1,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold);
            
            for i = 1:V
                for j = 1:V
                    lesser_object = LesserGreen(spin,i,j,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold);
                    lesser_green_matrix(i,j) = lesser_object;
                    greater_object = GreaterGreen(spin,i,j,hubbard,n,perm,t_values,w_values,combine_zero,chop_threshold);
                    greater_green_matrix(i,j) = greater_object;
                end
            end

            obj.lesser_green_matrix = lesser_green_matrix;
            obj.greater_green_matrix = greater_green_matrix;

            obj.isexact = isexact; % true to use Lehmann weights and frequencies, false to use compressive sensing weights and frequencies
            obj.perm = perm; % indices for compressive sensing

            obj.k_value = k_value;
            obj.V = V;
            obj.spin = spin;
            obj.hubbard = hubbard;
            obj.n = n;
            obj.t_values = t_values;
            obj.w_values = w_values;
            obj.combine_zero = combine_zero;
            obj.chop_threshold = chop_threshold;

        end

        function product_over_alpha_lehmann = form_product_over_alpha_lehmann(obj,z)
            matrix_element_lesser = obj.lesser_green_matrix(1,1); % arbitrary site indices
            product_over_alpha_lehmann = 1;
            for alpha = 1:length(matrix_element_lesser.angular_frequency_differences)
                product_over_alpha_lehmann = product_over_alpha_lehmann * (z - matrix_element_lesser.angular_frequency_differences(alpha));
            end
            return
        end

        function product_over_beta_lehmann = form_product_over_beta_lehmann(obj,z)
            matrix_element_greater = obj.greater_green_matrix(1,1); % arbitrary site indices
            product_over_beta_lehmann = 1;
            for beta = 1:length(matrix_element_greater.angular_frequency_differences)
                product_over_beta_lehmann = product_over_beta_lehmann * (z - matrix_element_greater.angular_frequency_differences(beta));
            end
            return
        end

        function product_over_alpha_compressive = form_product_over_alpha_compressive(obj,z)
            w_differences = [];
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element_lesser = obj.lesser_green_matrix(i,j);
                    for index = 1:length(matrix_element_lesser.approximate_angular_frequency_differences)
                        if ~ismember(round(matrix_element_lesser.approximate_angular_frequency_differences(index),1),round(w_differences,1))
                            w_differences(end+1) = matrix_element_lesser.approximate_angular_frequency_differences(index);
                        end
                    end
                end
            end
            product_over_alpha_compressive = 1;
            for alpha = 1:length(w_differences)
                product_over_alpha_compressive = product_over_alpha_compressive * (z + w_differences(alpha));
            end
            return
        end

        function product_over_beta_compressive = form_product_over_beta_compressive(obj,z)
            w_differences = [];
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element_greater = obj.greater_green_matrix(i,j);
                    for index = 1:length(matrix_element_greater.approximate_angular_frequency_differences)
                        if ~ismember(round(matrix_element_greater.approximate_angular_frequency_differences(index),1),round(w_differences,1))
                            w_differences(end+1) = matrix_element_greater.approximate_angular_frequency_differences(index);
                        end
                    end
                end
            end
            product_over_beta_compressive = 1;
            for beta = 1:length(w_differences)
                product_over_beta_compressive = product_over_beta_compressive * (z - w_differences(beta));
            end
            return
        end

        function N = form_numerator(obj,z)
            if obj.isexact == true % use frequencies from Lehmann
                product_over_beta = obj.form_product_over_beta_lehmann(z);
                product_over_alpha = obj.form_product_over_alpha_lehmann(z);
            else % use frequencies from compressive sensing
                product_over_beta = obj.form_product_over_beta_compressive(z);
                product_over_alpha = obj.form_product_over_alpha_compressive(z);
            end
            N = product_over_alpha * product_over_beta;
            return
        end

        function X_a_sum = form_X_a_sum_compressive(obj,matrix_element,i,j,z)
            X_a_sum = 0;
            for a = 1:length(matrix_element.approximate_angular_frequency_differences)
                lesser_weights = matrix_element.approximate_weights;
                X_a_sum = X_a_sum + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * lesser_weights(a);
                
                product_over_alpha_not_a = 1;
                for alpha = 1:length(matrix_element.approximate_angular_frequency_differences)
                    if alpha ~= a
                        product_over_alpha_not_a = product_over_alpha_not_a * (z + matrix_element.approximate_angular_frequency_differences(alpha));
                    end
                end
                X_a_sum = X_a_sum * product_over_alpha_not_a;
            end
            product_over_beta = obj.form_product_over_beta_compressive(z);
            X_a_sum = X_a_sum * product_over_beta;
            return
        end

        function X_a_sum = form_X_a_sum_lehmann(obj,matrix_element,i,j,z)
            X_a_sum = 0;
            for a = 1:length(matrix_element.angular_frequency_differences)
                lesser_weights = matrix_element.weights;
                X_a_sum = X_a_sum + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * lesser_weights(a);
                product_over_alpha_not_a = 1;
                for alpha = 1:length(matrix_element.angular_frequency_differences)
                    if alpha ~= a
                        product_over_alpha_not_a = product_over_alpha_not_a * (z - matrix_element.angular_frequency_differences(alpha));
                    end
                end
                X_a_sum = X_a_sum * product_over_alpha_not_a;
            end
            product_over_beta = obj.form_product_over_beta_lehmann(z);
            X_a_sum = X_a_sum * product_over_beta;
            return
        end

        function X_b_sum = form_X_b_sum_compressive(obj,matrix_element,i,j,z)
            X_b_sum = 0;
            for b = 1:length(matrix_element.approximate_angular_frequency_differences)
                greater_weights = matrix_element.approximate_weights;
                X_b_sum = X_b_sum + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * greater_weights(b);
                product_over_beta_not_b = 1;
                for beta = 1:length(matrix_element.approximate_angular_frequency_differences)
                    if beta ~= b
                        product_over_beta_not_b = product_over_beta_not_b * (z - matrix_element.approximate_angular_frequency_differences(beta));
                    end
                end
                X_b_sum = X_b_sum * product_over_beta_not_b;
            end
            product_over_alpha = obj.form_product_over_alpha_compressive(z);
            X_b_sum = X_b_sum * product_over_alpha;
            return
        end

        function X_b_sum = form_X_b_sum_lehmann(obj,matrix_element,i,j,z)
            X_b_sum = 0;
            for b = 1:length(matrix_element.angular_frequency_differences)
                greater_weights = matrix_element.weights;
                X_b_sum = X_b_sum + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * greater_weights(b);
                product_over_beta_not_b = 1;
                for beta = 1:length(matrix_element.angular_frequency_differences)
                    if beta ~= b
                        product_over_beta_not_b = product_over_beta_not_b * (z - matrix_element.angular_frequency_differences(beta));
                    end
                end
                X_b_sum = X_b_sum * product_over_beta_not_b;
            end
            product_over_alpha = obj.form_product_over_alpha_lehmann(z);
            X_b_sum = X_b_sum * product_over_alpha;
            return
        end

        function X_a_sum = select_X_a_sum(obj,matrix_element,i,j,z)
            if obj.isexact == true 
                X_a_sum = obj.form_X_a_sum_lehmann(matrix_element,i,j,z);
            else 
                X_a_sum = obj.form_X_a_sum_compressive(matrix_element,i,j,z);
            end
            return
        end

        function X_b_sum = select_X_b_sum(obj,matrix_element,i,j,z)
            if obj.isexact == true 
                X_b_sum = obj.form_X_b_sum_lehmann(matrix_element,i,j,z);
            else 
                X_b_sum = obj.form_X_b_sum_compressive(matrix_element,i,j,z);
            end
            return
        end

        function D = form_denominator(obj,z)

            X_a_sum = 0;
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element = obj.lesser_green_matrix(i,j);
                    X_a_sum = obj.select_X_a_sum(matrix_element,i,j,z);
                end
            end

            X_b_sum = 0;
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element = obj.greater_green_matrix(i,j);
                    X_b_sum = obj.select_X_b_sum(matrix_element,i,j,z);
                end
            end

            D = X_a_sum + X_b_sum;
            return

        end

        function inverse_retarded_green_as_a_ratio = compute(obj,z)
            N = obj.form_numerator(z);
            D = obj.form_denominator(z);
            inverse_retarded_green_as_a_ratio = N/D;
            return
        end

        function retarded_green_as_a_ratio = retarded(obj,z)
            N = obj.form_numerator(z);
            D = obj.form_denominator(z);
            retarded_green_as_a_ratio = D/N;
            return
        end
    end
end