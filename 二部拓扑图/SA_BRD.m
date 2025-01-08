function IntegratedCSMA
    % 输入参数
    Input = [20; 2000; 10000; 1];
    N = Input(1); % 节点数量是20
    Packet_size = Input(2) * 8; % 数据包大小（位）
    T = Input(3) * 10^-3; % 模拟时间（秒）
    iterations = 25000; % 模拟轮数

    % 数据率和初始化
    data_rate = 6 * 10^6; % 6 Mbps
    packet_time = Packet_size / data_rate; % 数据包时间
    slot_size = 9 * 10^-6; % 时隙大小

    % 初始化退避窗口和随机计时器
    CW_min = 15;
    r = rand(N, 1) * CW_min; % 初始化计时器随机值

    % 干扰图和收益函数
    interference_graph = create_bipartite_interference_graph(N);
    U = @(s) log(s + 1e-6); % 收益函数
    r_min = 0; % r 的最小值
    r_max = 10; % r 的最大值
    beta = 1.0; % 调整参数
    alpha = 0.5;

    % 模拟状态变量
    total_time = 0; 
    good_time = 0;

    while total_time < T
        % 找到最小计时器节点
        [M, I] = min(r);
        collision_nodes = find(r == M);
        collision_flag = numel(collision_nodes) > 1;

        if ~collision_flag
            % 无碰撞情况
            good_time = good_time + packet_time;
        end

        % 计算服务率并更新计时器
        exp_r = exp(r); 
        s = calculate_service_rate(exp_r, interference_graph);

        for node = collision_nodes'
            % 对每个冲突节点或成功传输节点执行更新
            r(node) = best_response(r(node), s(node), alpha, r_min, r_max, U, beta);
        end

        % 更新计时器并防止为负值
        r = max(0, r - M);
        total_time = total_time + slot_size;
    end

    % 输出结果
    disp_results(r, interference_graph, good_time, total_time, N);
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

function r_new = best_response(r_i, s_i, alpha, r_min, r_max, U, beta)
    % 计算最优响应
    utility = U(s_i);
    service_rate_derivative = s_i * (1 - s_i);
    gradient = service_rate_derivative * (utility - r_i / beta);
    r_new = r_i + alpha * gradient;
    r_new = max(r_min, min(r_new, r_max));
end

function disp_results(r, interference_graph, good_time, total_time, N)
    % 输出结果
    utility = good_time / total_time;
    exp_r = exp(r);
    final_service_rates = calculate_service_rate(exp_r, interference_graph);
    GAT = exp(sum(log(final_service_rates)) / N);
    fprintf('Utilization: %f\nGAT: %f\n', utility, GAT);
end

% 调用函数3e+07轮
for i = 1:3e+07
    IntegratedCSMA;
end