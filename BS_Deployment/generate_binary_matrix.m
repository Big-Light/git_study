function random_matrix = generate_random_matrix(n, s, m)
    % 创建全零矩阵
    random_matrix = zeros(n, s);
    
    % 在前 m 行中，随机选择一个位置为 1，其余位置为 0
    for i = 1:m
        random_col = randi([1, s]);  % 随机选择列索引
        random_matrix(i, random_col) = 1;
    end
    
    % 剩余的 n-m 行全部为 0
    
end

