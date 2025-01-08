function IntegratedCSMA
    % 输入参数
    Input = [20; 2000; 10000; 1];
    N = Input(1); % 节点数量是20
    Packet_size = Input(2) * 8; % 数据包大小（位）
    T = Input(3) * 10^-3; % 模拟时间（秒）
    Backoff_St = Input(4); % 随机退避策略
    iterations = 25000; % 模拟轮数

    % 数据率和初始化
    data_rate = 6 * 10^6; % 6 Mbps
    packet_time = Packet_size / data_rate; % 数据包时间
    slot_size = 9 * 10^-6; % 时隙大小

    % 初始化退避窗口和随机计时器
    CW_min = 15;
    r = [6; 6 * ones(N-1, 1)];

    % 干扰图和收益函数
    interference_graph = create_bipartite_interference_graph(N);
    U = @(s) log(s + 1e-6); % 收益函数
    beta = 1.0;
    alpha = 0.5;

    % 模拟状态变量
    total_time = 0; good_time = 0;
    collision_flag = false; 

    % 缓存计算的服务率
    exp_r = exp(r); 

    while total_time < T
        % 找到最小计时器节点
        [M, I] = min(r);
        collision_nodes = find(r == M);
        collision_flag = numel(collision_nodes) > 1;

        if ~collision_flag
            % 无碰撞情况
            good_time = good_time + packet_time;
            exp_r = exp(r); % 更新缓存
            s = calculate_service_rate(exp_r, interference_graph);
            r(I) =update_strategy(r(I), s(I), U, beta,alpha);
        else
            % 碰撞处理
            for node = collision_nodes'
                s = calculate_service_rate(exp_r, interference_graph);
                r(node) = update_strategy(r(node), s(node), U, beta,alpha);
            end
        end

        % 更新计时器
        r = max(0, r - slot_size);
        total_time = total_time + packet_time;
            % 打印最优的 r 值
        fprintf('Optimal r values for each node:\n');
        for i = 1:N
            fprintf('Node %d: %f\n', i, r(i));
        end
    end

    % 结果输出
    disp_results(r, interference_graph, good_time, total_time, N);
end

function interference_graph = create_bipartite_interference_graph(N)
    % N 是节点数量，对于二部图，N 应该是 20
    if N ~= 20
        error('节点数量必须是 20');
    end
    
    % 初始化干扰图
    interference_graph = zeros(N, N); % 创建一个全 0 的矩阵
    
    % 将节点分为两组，每组 10 个节点
    group1 = 1:10;
    group2 = 11:20;
    
    % 连接两组中的节点
    for i = group1
        for j = group2
            interference_graph(i, j) = 1;
            interference_graph(j, i) = 1;
        end
    end
    
    % 确保有100个链接
    num_links = sum(interference_graph(:));
    while num_links < 100
        % 随机选择两个节点
        node1 = randi(N);
        node2 = randi(N);
        % 如果这两个节点之间没有链接，则添加一个链接
        if node1 ~= node2 && interference_graph(node1, node2) == 0
            interference_graph(node1, node2) = 1;
            interference_graph(node2, node1) = 1;
            num_links = num_links + 1;
        end
    end
end


function s = calculate_service_rate(exp_r, interference_graph)
    num_nodes = size(interference_graph, 1);
    s = zeros(num_nodes, 1);
    
    for i = 1:num_nodes
        conflicting_nodes = find(interference_graph(i, :) == 1);
        feasible = exp_r(i) / (exp_r(i) + sum(exp_r(conflicting_nodes)));
        s(i) = feasible;
    end
end

function r_new = update_strategy(r_i, s_i, U, beta, alpha)
    % 更新策略，使用Jacobi dynamics
    % r_i: 当前策略
    % s_i: 当前服务率
    % U: 效用函数
    % beta: 调整参数
    % alpha: 平滑参数

    % 计算最佳响应策略
    BR = best_response(r_i, s_i, U, beta);

    % 更新策略
    r_new = r_i + alpha * (BR - r_i);
end

function BR = best_response(r_i, s_i, U, beta)
    % 最佳响应策略
    % r_i: 当前策略
    % s_i: 当前服务率
    % U: 效用函数
    % beta: 调整参数

    % 计算梯度
    gradient = beta * U(s_i);

    % 更新最佳响应策略
    BR = r_i + gradient;
end
function disp_results(r, interference_graph, good_time, total_time, N)
    utility = good_time / total_time;
    final_service_rates = calculate_service_rate(exp(r), interference_graph);
    GAT = exp(sum(log(final_service_rates)) / N);
    
    fprintf('Utilization: %f\nGAT: %f\n', utility, GAT);
end



% 调用函数1e+07轮
for i = 1:1e+07
    IntegratedCSMA;
end