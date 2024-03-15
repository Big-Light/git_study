function random_vector = generate_random_vector(n, m)
    % 创建全零向量
    random_vector = zeros(n, 1);
    
    % 在前 m 个元素中设置为 1
    random_vector(1:m) = 1;
end

