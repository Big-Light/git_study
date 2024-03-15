function Matrix_g_two = function_g_two(Matrix_g , Matrix_X , test_points , candidate_Site)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%初始化
Matrix_g_two= zeros(test_points,candidate_Site);

for t = 1 : test_points
    for j = 1 : candidate_Site
        for n = 1 :candidate_Site
            if n == j
                continue
            else
                Matrix_g_two(t , j) = Matrix_g_two(t , j) + Matrix_g(t , n) * Matrix_X(n , n);
            end
        end
    end
end

end