function [result] = function_angle(vector1, vector2)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%计算两个向量的点积
    dot_product = dot(vector1, vector2);
    
    %计算向量的模
    norm_a = norm(vector1);
    norm_b = norm(vector2);
    
    %计算向量的余弦值
    cosine_theta = dot_product / (norm_a * norm_b);

    %计算向量的弧度
    theta_radians = acos(cosine_theta);

    %弧度转化为角度
    result = rad2deg(theta_radians);
end