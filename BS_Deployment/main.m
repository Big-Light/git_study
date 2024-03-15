candidate_Site = 5;     %基站候选点S
test_points = 3;        %测试点t
deploy_num = 3;         %实际部署基站数目，部署预算

%关联决策矩阵A
A = generate_binary_matrix(test_points, candidate_Site, deploy_num);
% disp("A:")
% disp(A);


%生成测试点坐标
t_x = rand(test_points, 1)*10;
t_y = rand(test_points, 1)*10;

%生成候选站点坐标
Cs_x = rand(candidate_Site, 1)*10;
Cs_y = rand(candidate_Site, 1)*10;

pt = [t_x,t_y];      %测试点坐标
Cs = [Cs_x,Cs_y];    %候选站点坐标

%一些常量初始化
alpha_value = 3;        %g
fk = 3.5;
v_lightspeed = 3.00e8;
sigma_sk = 7;
P = 20;                  %r'
N0 = 1;
W = 100;
lambda = 1;             %b
delta = 0.8;            %误差阈值


%一些向量矩阵初始化
Vector_v = zeros(test_points,candidate_Site);               %计算bt的向量vt
Matrix_g = zeros(test_points,candidate_Site);               %存储gtj的矩阵
Vector_x = generate_random_vector(candidate_Site,deploy_num); %部署决策向量
% disp('Vector_x:')
% disp(Vector_x);
Matrix_X = Vector_x * Vector_x';                            %求解Problem4的矩阵X
% [U, S, V] = svd(Matrix_X);
% 
% % 获取向量y
% y = sqrt(S(1,1)) * U(:,1);
% disp(y)
% disp('Matrix_X:')
% disp(Matrix_X)
auxiliar_Vector_y = zeros(test_points,candidate_Site)+ 1;    %初始化辅助变量y的矩阵
Vector_lt = zeros(test_points,1) + 1;
Vector_ls = zeros(candidate_Site,1) + 1;


%矩阵g
for t=1 : test_points       
    for j=1 : candidate_Site
        Matrix_g(t,j) = function_g(pt(t,1) , pt(t,2) , Cs(j,1) , Cs(j,2) , fk , v_lightspeed , alpha_value , sigma_sk);
    end
end
% disp("Matrix_g:")
% disp(Matrix_g)

%矩阵g'
Matrix_g_ = function_g_(Matrix_g , Vector_x , test_points , candidate_Site);
% disp("Matrix_g_:")
% disp(Matrix_g_)

%优化目标r'
Vector_r_t = function_r_(Matrix_g,Matrix_g_,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site);
% disp("Vector_r_t:")
% disp(Vector_r_t)

%向量v
for t=1 : test_points
    for j=1 : candidate_Site
        d = sqrt((pt(t,1) - Cs(j,1))^2 + (pt(t,2) - Cs(j,2))^2);
        %disp(d);
        fsd = function_sigma(d,alpha_value);
        %disp(fsd);
        Vector_v(t,j) = function_v(d,lambda,fsd,alpha_value);
    end
end
% disp("Vector_v:") 
% disp(Vector_v)          

%创建t个F矩阵
Matrix_F = zeros(candidate_Site, candidate_Site, t)+1.0;

for t=1 : test_points
    for i=1 : candidate_Site
        for j=1 : candidate_Site
            if i<j
                theta_j = function_angle(pt(t,:) , Cs(j,:));
                theta_i = function_angle(pt(t,:) , Cs(i,:));
                Matrix_F(i,j,t) = function_Matrix_F(Vector_v(t , i),Vector_v(t,j),theta_j,theta_i);
            end
        end
    end
end

% for t = 1:test_points
%     disp(['F matrix for t = ', num2str(t)]);
%     disp(Matrix_F(: , : , t));
% end


%tolerance = 1.0e-6;   %收敛阈值

%——————————————————————————————————————————————————————————————————————————
%在考虑定位性能要求下总速率与迭代次数的关系
max_iter = 15;     %最大迭代次数
Vector_iter = zeros(max_iter , 1);  %判断收敛向量
cvx_clear
% 迭代求解
% 
for iter = 1:max_iter 
    fprintf('第%d次迭代' , iter);
    % Step1————初始化y，A，求解更新 X
    cvx_begin 
        vector_w = zeros(test_points,1)+1.0;
        variable Matrix_X(candidate_Site,candidate_Site)
        expression totalRate
        expression gg(test_points,candidate_Site)
        expression Vector_r_t(test_points,1)
        expression temp(candidate_Site,1)
        %gg= function_g_(Matrix_g , Matrix_X ,test_points ,candidate_Site);
        for t = 1 : test_points
            for j = 1 : candidate_Site
                for n = 1 :candidate_Site
                    if n == j
                        continue
                    else
                        gg(t , j) = gg(t , j) + Matrix_g(t , n) * Matrix_X(n,n);
                    end
                end
            end
        end
%       Vector_r_t = function_r_(Matrix_g,gg,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site);
        ls = zeros(candidate_Site , 1) + 1;
        for t = 1 : test_points
            for j = 1 : candidate_Site
                temp(j,1)= W*A(t,j)*(1/log(2)*log(P*gg(t,j)+P*Matrix_g(t,j)+N0*W))+ ...
                    1/log(2)*log(auxiliar_Vector_y(t,j))...
                +1-auxiliar_Vector_y(t,j)*(P*gg(t,j)+N0*W);
            end
            Vector_r_t(t,1) = temp' * ls;
        end
        totalRate = Vector_r_t' * vector_w;
        maximize(totalRate )
        subject to
            for t=1:test_points
                trace(Matrix_X .* diag(Vector_v(t,:))) <= delta*trace(Matrix_F(: , : , t)*Matrix_X);
            end
            A <= Vector_lt*diag(Matrix_X)';
            trace(Matrix_X) <= deploy_num;
            Matrix_X >= 0;
            Matrix_X <= 1;
    cvx_end
     
    disp('Solution Matrix_X:');
    disp(Matrix_X);
    
%     Matrix_X = setMatrix(Matrix_X);
%     disp('after Solution Matrix_X:');
%     disp(Matrix_X);

    % Step2————根据上述求得的X推出y
    Matrix_g_ = function_g_two(Matrix_g , Matrix_X ,test_points ,candidate_Site);
    for t=1 : test_points
        for j=1 : candidate_Site
            auxiliar_Vector_y(t,j) = 1/(P*Matrix_g_(t , j) + N0*W); 
        end
    end
%     disp('auxiliar_Vector_y:')
%     disp(auxiliar_Vector_y)

    % Step3————求A
    cvx_begin quiet
        variable A(test_points,candidate_Site) 
        expression totalRate
        expression Vector_r_t(test_points,1)
        expression temp(candidate_Site , 1)
        %Vector_r_t = function_r_(Matrix_g,Matrix_g_,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site);
        ls = zeros(candidate_Site , 1) + 1;
        for t = 1 : test_points
            for j = 1 : candidate_Site
                temp(j,1)= W*A(t,j)*(1/log(2)*log(P*Matrix_g_(t,j)+P*Matrix_g(t,j)+N0*W))+ ...
                    1/log(2)*log(auxiliar_Vector_y(t,j))...
                +1-auxiliar_Vector_y(t,j)*(P*Matrix_g(t,j)+N0*W);
            end
            Vector_r_t(t,1) = temp' * ls;
        end
        totalRate = Vector_r_t' * vector_w;
        maximize(totalRate)
        subject to
            A <= Vector_lt*diag(Matrix_X)';
            A*Vector_ls <= Vector_lt;
            A >= 0;
            A <= 1;
    cvx_end

    disp('Solution A:');
    disp(A);

    disp("Vector_r_t:")
    disp(Vector_r_t)

%     检查收敛条件
    Vector_iter(iter,1) = Vector_r_t' * vector_w;
    
    fprintf('第%d次的总速率为：' , iter);
    disp(Vector_iter(iter,1))
%     if iter >= 2
%         if Vector_iter(iter,1) - Vector_iter(iter-1,1) < tolerance
%             disp(['Converged at iteration ', num2str(iter)]);
%             break;
%         end
%     end
% 
%     if iter == max_iter
%         disp('未收敛')
%     end
end

hold on;
x_coordinate = 1:1:max_iter ;    %绘图横坐标
plot(x_coordinate,Vector_iter,'-*b'); %线性，颜色，标记
title('下总速率与迭代次数的关系');
xlabel('Iteration Index');
ylabel('totalRate');

%——————————————————————————————————————————————————————————————————————————
%仅最大化速率不考虑定位性能
cvx_clear
% 迭代求解
max_iter2 = 15;
for iter = 1:max_iter2 
    fprintf('第%d次迭代' , iter);
    % Step1————初始化y，A，求解更新 X
    cvx_begin quiet
        vector_w = zeros(test_points,1)+1.0;
        variable Matrix_X(candidate_Site,candidate_Site)
        expression totalRate
        expression gg(test_points,candidate_Site)
        expression Vector_r_t(test_points,1)
        expression temp(candidate_Site,1)
        %gg= function_g_(Matrix_g , Matrix_X ,test_points ,candidate_Site);
        for t = 1 : test_points
            for j = 1 : candidate_Site
                for n = 1 :candidate_Site
                    if n == j
                        continue
                    else
                        gg(t , j) = gg(t , j) + Matrix_g(t , n) * Matrix_X(n,n);
                    end
                end
            end
        end
%       Vector_r_t = function_r_(Matrix_g,gg,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site);
        ls = zeros(candidate_Site , 1) + 1;
        for t = 1 : test_points
            for j = 1 : candidate_Site
                temp(j,1)= W*A(t,j)*(1/log(2)*log(P*gg(t,j)+P*Matrix_g(t,j)+N0*W))+ ...
                    1/log(2)*log(auxiliar_Vector_y(t,j))...
                +1-auxiliar_Vector_y(t,j)*(P*gg(t,j)+N0*W);
            end
            Vector_r_t(t,1) = temp' * ls;
        end
        totalRate = Vector_r_t' * vector_w;
        maximize(totalRate )
        subject to
            A <= Vector_lt*diag(Matrix_X)';
            trace(Matrix_X) <= deploy_num;
            Matrix_X >= 0;
            Matrix_X <= 1;
    cvx_end
     
    disp('Solution Matrix_X:');
    disp(Matrix_X);
    
%     Matrix_X = setMatrix(Matrix_X);
%     disp('after Solution Matrix_X:');
%     disp(Matrix_X);

    % Step2————根据上述求得的X推出y
    Matrix_g_ = function_g_two(Matrix_g , Matrix_X ,test_points ,candidate_Site);
    for t=1 : test_points
        for j=1 : candidate_Site
            auxiliar_Vector_y(t,j) = 1/(P*Matrix_g_(t , j) + N0*W); 
        end
    end
%     disp('auxiliar_Vector_y:')
%     disp(auxiliar_Vector_y)

    % Step3————求A
    cvx_begin quiet
        variable A(test_points,candidate_Site) 
        expression totalRate
        expression Vector_r_t(test_points,1)
        expression temp(candidate_Site , 1)
        %Vector_r_t = function_r_(Matrix_g,Matrix_g_,auxiliar_Vector_y,A,P,N0,W,test_points,candidate_Site);
        ls = zeros(candidate_Site , 1) + 1;
        for t = 1 : test_points
            for j = 1 : candidate_Site
                temp(j,1)= W*A(t,j)*(1/log(2)*log(P*Matrix_g_(t,j)+P*Matrix_g(t,j)+N0*W))+ ...
                    1/log(2)*log(auxiliar_Vector_y(t,j))...
                +1-auxiliar_Vector_y(t,j)*(P*Matrix_g(t,j)+N0*W);
            end
            Vector_r_t(t,1) = temp' * ls;
        end
        totalRate = Vector_r_t' * vector_w;
        maximize(totalRate)
        subject to
            A <= Vector_lt*diag(Matrix_X)';
            A*Vector_ls <= Vector_lt;
            A >= 0;
            A <= 1;
    cvx_end

    disp('Solution A:');
    disp(A);

    disp("Vector_r_t:")
    disp(Vector_r_t)

%     检查收敛条件
    Vector_iter2(iter,1) = Vector_r_t' * vector_w;
    fprintf('第%d次的总速率为：' , iter);
    disp(Vector_iter(iter,1))

end
%figure();

x_coordinate = 1:1:max_iter2 ;    %绘图横坐标
plot(x_coordinate,Vector_iter2,'-^g'); %线性，颜色，标记
title('不考虑定位性能要求下总速率与迭代次数的关系');
xlabel('Iteration Index');
ylabel('totalRate');
legend('考虑定位性能','不考虑定位性能');
hold off;
%——————————————————————————————————————————————————————————————————————————
