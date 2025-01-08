function IntegratedCSMA
    % 输入参数
    Input = [20; 2000; 10000; 1]; 
    N = Input(1); % 节点数量
    L = 28;
    Packet_size = Input(2) * 8; % 数据包大小（位）
    T = Input(3) * 10^-3; % 模拟时间（秒）
    Backoff_St = Input(4); % 随机退避策略
    iterations = 25000; % 模拟轮数

    % 数据率和数据包时间
    data_rate = 6 * 10^6; % 6 Mbps
    packet_time = Packet_size / data_rate; % 数据包时间（秒）
    slot_size = 9 * 10^-6; % 时隙大小（秒）

    % 初始化退避窗口和随机数
    CW_min = 15; % 最小退避窗口
    r = [5.5; 2 * ones(N-1, 1)]; % 调整退避计时器的初始值

    % 定义星型拓扑的干扰图
    interference_graph =  create_random_interference_graph(N, L);

    % 定义收益函数和初始策略
    U = @(s) log(s + 1e-6); % Proportional fairness，加小量避免log(0)
    beta = 1.0; % 调整参数

    % 定义效用函数的逆函数
    U_inv = @(U) 2.^U - 1e-6;
    % 定义步长函数
    a_i = @(t) 1 / t; % 示例步长函数

    % 模拟CSMA网络
    total_time = 0;
    count = 0;
    good_time = 0;
    CW = CW_min;
    collision_flag = 0;
    j = 1;
    simulation_count = 0;
    collision_index = zeros(N, 1); % 初始化碰撞索引数组

    while total_time < T
        [M, I] = min(r); % find the node with the minimum counter and index of node
        simulation_count = exp(sqrt(simulation_count + 1));%设置步长为exp(sqrt(t))

        for i = 1:N   % check if there are more than one nodes with same minimum counter
            if (M == r(i))
                count = count + 1;
                collision_index(j) = i;
                j = j + 1;
            end
        end

        if count > 1
            collision_flag = 1; % collision occured
        end

        if collision_flag ~= 1
            good_time = good_time + packet_time;
            % Calculate service rate and update strategy using best response dynamics
            s = calculate_service_rate(r, interference_graph);
            
            r(I) = best_response(r(I), s(I), U_inv, beta, a_i(simulation_count));
        else
            CW = CW + 2;
            for i = 1:N
                if (M == r(i))
                    s = calculate_service_rate(r, interference_graph);
                    r(i) = best_response(r(i), s(i), U_inv, beta, a_i(simulation_count));
                end
            end
        end

        total_time = total_time + packet_time;
        for i = 1:N
            for j = 1:length(collision_index)
                if (i ~= collision_index(j))
                    r(i) = max(0, r(i) - slot_size); % 确保 r(i) 不为负
                end
            end
        end
        count = 0;
        collision_flag = 0;
        j = 1; % 重置碰撞索引计数器
    end

    utility = good_time / total_time;

    % 计算GAT
    s = calculate_service_rate(r, interference_graph); % 计算最终的服务率
    GAT = calculate_GAT(s); % 计算几何平均值

    % 计算最终的服务率
    final_service_rates = calculate_service_rate(r, interference_graph);

    % 打印各节点的服务率
    fprintf('Service rates for each node:\n');
    for i = 1:N
        fprintf('Node %d: %f\n', i, final_service_rates(i));
    end

    % 打印最优的 r 值
    fprintf('Optimal r values for each node:\n');
    for i = 1:N
        fprintf('Node %d: %f\n', i, r(i));
    end

    % 输出结果
    fprintf('Number of Nodes: %d ; Packet Size: %d ; Simulation Time(s): %d ; Backoff Strategy: %d \n Utilization: %f\n', N, Input(2), T, Backoff_St, utility);
    fprintf('Geometric Average of user Throughput (GAT): %f\n', GAT);
end

function interference_graph = create_random_interference_graph(N, L)
    % N 是节点数量，L 是边的数量
    if N ~= 20 || L ~= 28
        error('节点数量必须是 20，边的数量必须是 28');
    end
    
    % 初始化干扰图
    interference_graph = zeros(N, N); % 创建一个全 0 的矩阵
    
    % 随机生成28条边
    while sum(interference_graph(:)) < 2*L
        % 随机选择两个不同的节点
        u = randi(N);
        v = randi(N);
        while u == v || interference_graph(u, v) == 1
            u = randi(N);
            v = randi(N);
        end
        
        % 添加边
        interference_graph(u, v) = 1;
        interference_graph(v, u) = 1;
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

function GAT = calculate_GAT(s)
    % 确保服务率数组中的所有值都是正数
    s = max(s, 1e-10); % 添加一个小的正数以避免零值

    % 计算几何平均值
    GAT = exp(sum(log(s)) / length(s));
end


function r_new = best_response(r_i, s_i_current, U_inv, beta, a_i)
    % r_i: 当前的 r_i 值
    % s_i_current: 瞬时服务率
    % U_inv: 效用函数的逆函数
    % beta: 调整参数
    % a_i: 步长函数

    % 计算效用函数的逆
    U_inv_value = U_inv(s_i_current/beta);

    % 更新 r_i
    r_new = r_i + a_i * (U_inv_value - r_i);

     r_new = max(0, min(r_new, 10)); % 确保 r_new 不为负
end

% 调用函数1e+07轮
for i = 1:1e+07
    IntegratedCSMA;
end