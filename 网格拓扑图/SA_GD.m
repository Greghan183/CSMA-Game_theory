function interference_graph = create_grid_interference_graph(N)
    if N ~= 25
        error('节点数量必须是 25');
    end
    
    % 初始化干扰图
    interference_graph = zeros(N, N);
    
    % 生成 5x5 网格
    for i = 1:N
        % 计算节点的行和列
        row = ceil(i / 5);
        col = mod(i - 1, 5) + 1;
        
        % 左边节点
        if col > 1
            left_node = i - 1;
            interference_graph(i, left_node) = 1;
            interference_graph(left_node, i) = 1;
        end
        
        % 右边节点
        if col < 5
            right_node = i + 1;
            interference_graph(i, right_node) = 1;
            interference_graph(right_node, i) = 1;
        end
        
        % 上边节点
        if row > 1
            top_node = i - 5;
            interference_graph(i, top_node) = 1;
            interference_graph(top_node, i) = 1;
        end
        
        % 下边节点
        if row < 5
            bottom_node = i + 5;
            interference_graph(i, bottom_node) = 1;
            interference_graph(bottom_node, i) = 1;
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

function IntegratedCSMA
    % 输入参数
    Input = [25; 2000; 10000; 1];
    N = Input(1); % 节点数量
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
    r = [5; 5 * ones(N-1, 1)];

    % 干扰图和收益函数
    interference_graph = create_grid_interference_graph(N);
    U = @(s) log(s + 1e-6); % 收益函数
    r_min = 0; % r 的最小值
    r_max = 10; % r 的最大值
    beta = 1.0; % 调整参数
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
            r(I) = best_response(r(I), s(I), alpha, r_min, r_max, U, beta);
        else
            % 碰撞处理
            for node = collision_nodes'
                s = calculate_service_rate(exp_r, interference_graph);
                r(node) = best_response(r(node), s(node),alpha, r_min, r_max, U, beta);
            end
        end

        % 更新计时器
        r = max(0, r - slot_size);
        total_time = total_time + packet_time;
            % 打印最优的 r 值
        %fprintf('Optimal r values for each node:\n');
        %for i = 1:N
            %fprintf('Node %d: %f\n', i, r(i));
        %end
    end

    % 结果输出
    disp_results(r, interference_graph, good_time, total_time, N);
end

function disp_results(r, interference_graph, good_time, total_time, N)
    utility = good_time / total_time;
    final_service_rates = calculate_service_rate(exp(r), interference_graph);
    GAT = exp(sum(log(final_service_rates)) / N);
    fprintf('Utilization: %f\nGAT: %f\n', utility, GAT);
end

function r_new = best_response(r_i, s_i, alpha, r_min, r_max, U, beta)
    % r_i: 当前的 r_i 值
    % s_i: 服务率
    % alpha: 步长
    % r_min: r 的最小值
    % r_max: r 的最大值
    % U: 效用函数
    % beta: 平滑参数
    
    % 计算效用函数 U
    utility = U(s_i);
    
    % 计算服务率 s_i 的偏导数
    service_rate_derivative = s_i * (1 - s_i);

    % 计算梯度
    gradient = service_rate_derivative * (utility - r_i / beta);

    % 更新 r_i
    r_new = r_i + alpha * gradient;

    % 确保 r_new 在合理的范围内
    r_new = max(r_min, min(r_new, r_max));
end

% 调用函数
IntegratedCSMA;