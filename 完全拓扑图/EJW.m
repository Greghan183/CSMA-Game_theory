function IntegratedCSMA
    % 输入参数
    Input = [5; 2000; 10000; 1];
    N = Input(1); % 节点数量
    Packet_size = Input(2) * 8; % 数据包大小（位）
    T = Input(3) * 10^-3; % 模拟时间（秒）
    Backoff_St = Input(4); % 随机退避策略
    iterations = 1e+07; % 模拟轮数

    % 数据率和数据包时间
    data_rate = 6 * 10^6; % 6 Mbps
    packet_time = Packet_size / data_rate; % 数据包时间（秒）
    slot_size = 9 * 10^-6; % 时隙大小（秒）

    % 初始化退避窗口和随机数
    CW_min = 15; % 最小退避窗口
    r = [6; 2 * ones(N-1, 1)]; % 调整退避计时器的初始值

    % 定义完全型拓扑的干扰图
    interference_graph = create_complete_interference_graph(N);

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
        simulation_count = simulation_count + 1;%设置固定步长是1

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

function interference_graph = create_complete_interference_graph(N)
    % N 是节点数量，对于完整图，N 应该是 5
    if N ~= 5
        error('节点数量必须是 5');
    end
    
    % 初始化干扰图
    interference_graph = ones(N, N); % 创建一个全 1 的矩阵
    
    % 将对角线元素设置为 0，因为节点不与自身相连
    for i = 1:N
        interference_graph(i, i) = 0;
    end
end

function s = calculate_service_rate(r, interference_graph)
    num_nodes = size(interference_graph, 1);
    exp_r = exp(r); % 确保exp_r是一个行向量
    feasible_schedules = get_feasible_schedules(interference_graph);
    
    % 初始化分子和分母
    numerator = zeros(num_nodes, 1);
    denominator = 0;
    
    % 计算分子和分母
    for i = 1:size(feasible_schedules, 1)
        schedule = feasible_schedules(i, :);
        
        % 只计算非零元素的乘积
        non_zero_indices = schedule > 0;
        prod_term = prod(exp_r(non_zero_indices) .* schedule(non_zero_indices)');
        
        numerator = numerator + prod_term * schedule'; % 分别更新每个节点的贡献
        denominator = denominator + prod_term;
    end
    
    % 计算服务率
    if denominator > 0
        s = numerator / denominator; % 添加一个小的正数以避免0值
    else
        s = ones(num_nodes, 1); % 防止分母为0
    end
end

function GAT = calculate_GAT(s)
    % 确保服务率数组中的所有值都是正数
    s = max(s, 1e-10); % 添加一个小的正数以避免零值

    % 计算几何平均值
    GAT = exp(sum(log(s)) / length(s));
end

function feasible_schedules = get_feasible_schedules(interference_graph)
    num_nodes = size(interference_graph, 1);
    feasible_schedules = []; % 初始化一个空矩阵来存储可行的调度

    for i = 0:(2^num_nodes - 1)
        binary_i = dec2bin(i, num_nodes) - '0'; % 将整数转换为二进制向量
        feasible = true; % 假设当前调度是可行的

        % 检查是否有干扰
        for j = 1:num_nodes
            for k = 1:num_nodes
                if binary_i(j) == 1 && binary_i(k) == 1 && interference_graph(j, k) == 1
                    feasible = false; % 如果发现干扰，则当前调度不可行
                    break;
                end
            end
            if ~feasible
                break;
            end
        end

        if feasible
            feasible_schedules = [feasible_schedules; binary_i]; % 存储可行的调度
        end
    end
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


% 调用函数3e+07轮
for i = 1:3e+07
    IntegratedCSMA;
end