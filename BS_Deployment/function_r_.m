function [Vector_r_] = function_r_(Matrix_g,Matrix_g_,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site)


%初始化
Vector_r_ = zeros(test_points , 1);  
temp = zeros(candidate_Site , 1);
ls = zeros(candidate_Site , 1) + 1;
for t = 1 : test_points
    for j = 1 : candidate_Site
        temp(j,1)= W*A(t,j)*(1/log(2)*log(P*Matrix_g_(t,j)+P*Matrix_g(t,j)+N0*W))+ ...
            1/log(2)*log(auxiliar_Vector_y(t,j))...
        +1-auxiliar_Vector_y(t,j)*(P*Matrix_g_(t,j)+N0*W);
    end
    Vector_r_(t,1) = temp' * ls;
end
end