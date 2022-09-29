classdef InverseRetardedGreen
    properties
        isexact % true for Lehmann, false for compressive sensing
        lesser_angular_frequency_differences % no site dependence
        greater_angular_frequency_differences % no site dependence
        perm % for consistency, use the same indices for all applications of compressive sensing 
        alpha_max % for product over alpha
        beta_max % for product over beta
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

            lesser_green = lesser_green_matrix(1,1);
            greater_green = greater_green_matrix(1,1);

            obj.alpha_max = length(lesser_green.eigenvalues_for_hubbard_minus);
            obj.beta_max = length(greater_green.eigenvalues_for_hubbard_plus);

            obj.isexact = isexact; % true to use Lehmann weights and frequencies, false to use compressive sensing weights and frequencies
            obj.perm = perm; % indices for compressive sensing

            if obj.isexact == true % Lehmann
                obj.lesser_angular_frequency_differences = lesser_green.angular_frequency_differences;
                obj.greater_angular_frequency_differences = greater_green.angular_frequency_differences;
            else % compressive sensing
                obj.lesser_angular_frequency_differences = lesser_green.approximate_angular_frequency_differences;
                obj.greater_angular_frequency_differences = greater_green.approximate_angular_frequency_differences;
            end

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

        function product_over_alpha = form_product_over_alpha(obj,z)
            product_over_alpha = 1;
            for alpha = 1:obj.alpha_max
                product_over_alpha = product_over_alpha * (z - obj.lesser_angular_frequency_differences(alpha));
            end
            return
        end

        function product_over_beta = form_product_over_beta(obj,z)
            product_over_beta = 1;
            for beta = 1:obj.beta_max
                product_over_beta = product_over_beta * (z - obj.greater_angular_frequency_differences(beta));
            end
            return
        end

        function product_over_alpha_not_a = form_product_over_alpha_not_a(obj,z,a)
            product_over_alpha_not_a = 1;
            for alpha = 1:obj.alpha_max
                if alpha ~= a
                    product_over_alpha_not_a = product_over_alpha_not_a * (z - obj.lesser_angular_frequency_differences(alpha));
                end
            end
            return
        end

        function product_over_beta_not_b = form_product_over_beta_not_b(obj,z,b)
            product_over_beta_not_b = 1;
            for beta = 1:obj.beta_max
                if beta ~= b
                    product_over_beta_not_b = product_over_beta_not_b * (z - obj.greater_angular_frequency_differences(beta));
                end
            end
            return
        end

        function X_a = form_X_a(obj,a)
            X_a = 0;
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element = obj.lesser_green_matrix(i,j);
                    if obj.isexact == true % use weight from Lehmann
                        lesser_weights = matrix_element.weights;
                        X_a = X_a + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * lesser_weights(a);
                    else % use weight from compressive sensing
                        lesser_weights = matrix_element.approximate_weights;
                        X_a = X_a + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * lesser_weights(a);
                    end
                end
            end
            return
        end

        function X_b = form_X_b(obj,b)
            X_b = 0;
            for i = 1:obj.V
                for j = 1:obj.V
                    matrix_element = obj.greater_green_matrix(i,j);
                    if obj.isexact == true % use weight from Lehmann
                        greater_weights = matrix_element.weights;
                        X_b = X_b + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * greater_weights(b);
                    else % use weight from compressive sensing
                        greater_weights = matrix_element.approximate_weights;
                        X_b = X_b + (1/obj.V) * cos( abs(obj.k_value) * abs(i - j) ) * greater_weights(b);
                    end
                end
            end
            return
        end

        function N = form_numerator(obj,z)
            product_over_alpha = obj.form_product_over_alpha(z);
            product_over_beta = obj.form_product_over_beta(z);
            N =  product_over_alpha * product_over_beta;
            return
        end

        function D = form_denominator(obj,z)
            product_over_alpha = obj.form_product_over_alpha(z);
            product_over_beta = obj.form_product_over_beta(z);

            sum_over_b = 0;
            for b = 1:obj.beta_max
                X_b = obj.form_X_b(b);
                product_over_beta_not_b = obj.form_product_over_beta_not_b(z,b);
                sum_over_b = sum_over_b + X_b * product_over_beta_not_b;
            end
            sum_over_b = sum_over_b * product_over_alpha;

            sum_over_a = 0;
            for a = 1:obj.alpha_max
                X_a = obj.form_X_a(a);
                product_over_alpha_not_a = obj.form_product_over_alpha_not_a(z,a);
                sum_over_a = sum_over_a + X_a * product_over_alpha_not_a;
            end
            sum_over_a = sum_over_a * product_over_beta;

            D = sum_over_a + sum_over_b;
            return
        end

        function inverse_retarded_green_as_a_ratio = compute(obj,z)
            N = obj.form_numerator(z);
            D = obj.form_denominator(z);
            inverse_retarded_green_as_a_ratio = N/D;
            return
        end
    end
end