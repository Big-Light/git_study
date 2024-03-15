function Matrix_g_ = function_g_(Matrix_g , Vector_x , test_points , candidate_Site)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%初始化
Matrix_g_= zeros(test_points,candidate_Site);

for t = 1 : test_points
    for j = 1 : candidate_Site
        for n = 1 :candidate_Site
            if n == j
                continue
            else
                Matrix_g_(t , j) = Matrix_g_(t , j) + Matrix_g(t , n) * Vector_x(n);
            end
        end
    end
end

end