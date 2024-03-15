function [result] = function_Matrix_F(v1,v2,angle1,angle2)


result = v1 * v2 *(sin(angle1 - angle2))^2;
end