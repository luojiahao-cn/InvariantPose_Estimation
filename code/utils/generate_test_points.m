function test_points = generate_test_points(params, num_points, method)
% GENERATE_TEST_POINTS 生成多个测试点（p_true 和 theta_true 的组合）
% 输入：
%   params     - 实验参数结构体
%   num_points - 测试点数量
%   method     - 生成方法：
%                'random': 在工作空间内随机生成
%                'grid': 在工作空间内网格采样（如果num_points是立方数）
%                'uniform_sphere': 在球面上均匀采样位置，随机旋转
% 输出：
%   test_points - 结构体数组，每个元素包含：
%                 - p_true: 3×1 真实位置
%                 - theta_true: 3×1 真实旋转向量

if nargin < 3
    method = 'random';
end

workspace_center = params.workspace.center;
workspace_radius = params.workspace.radius - 0.01;

test_points = struct('p_true', {}, 'theta_true', {});

switch method
    case 'fixed'
        % 将位姿选用params中的ground_truth
        for i = 1:num_points
            test_points(i).p_true = params.ground_truth.p_true;
            test_points(i).theta_true = params.ground_truth.theta_true;
        end
    case 'random'
        % 在工作空间内随机生成位置和旋转
        i = 1;
        while i <= num_points
            % 随机位置（在工作空间球体z>=0部分内）
            r = workspace_radius * (rand()^(1/3)); % 均匀分布在球体内
            theta = 2*pi*rand();
            phi = acos(rand()); % phi限制在[0, pi/2]保证z>=0
            p_true_candidate = workspace_center + r * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)];
            if p_true_candidate(3) >= workspace_center(3)
                p_true = p_true_candidate;
                % 随机旋转向量（限制在合理范围内）
                theta_true = pi * (2*rand(3,1) - 1); % [-pi, pi]
                test_points(i).p_true = p_true;
                test_points(i).theta_true = theta_true;
                i = i + 1;
            end
        end
        
    case 'grid'
        % 网格采样（适用于立方数）
        n = round(num_points^(1/3));
        if n^3 ~= num_points
            warning('num_points不是立方数，将使用最接近的立方数: %d', n^3);
        end
        
        x = linspace(-workspace_radius, workspace_radius, n) + workspace_center(1);
        y = linspace(-workspace_radius, workspace_radius, n) + workspace_center(2);
        z = linspace(-workspace_radius, workspace_radius, n) + workspace_center(3);
        
        idx = 1;
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if idx > num_points
                        break;
                    end
                    p_true = [x(i); y(j); z(k)];
                    theta_true = pi * (2*rand(3,1) - 1);
                    test_points(idx).p_true = p_true;
                    test_points(idx).theta_true = theta_true;
                    idx = idx + 1;
                end
                if idx > num_points
                    break;
                end
            end
            if idx > num_points
                break;
            end
        end
        
    case 'uniform_sphere'
        % 在球面上均匀采样位置，随机旋转
        for i = 1:num_points
            % 球面上均匀采样
            theta = 2*pi*rand();
            phi = acos(2*rand() - 1);
            p_true = workspace_center + workspace_radius * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)];
            
            % 随机旋转向量
            theta_true = pi * (2*rand(3,1) - 1);
            
            test_points(i).p_true = p_true;
            test_points(i).theta_true = theta_true;
        end
        
    otherwise
        error('未知的生成方法: %s', method);
end

end

