function [B, gradB] = dipole_b_and_gradb(position_vec, moment_unit, moment_mag, mu_0)
% DIPOLE_B_AND_GRADB 计算磁偶极子的磁场及其梯度
%
% 输入参数：
%   position_vec  - 观测点位置矢量 [x, y, z] (m)
%   moment_unit   - 磁偶极矩单位方向矢量 [mx, my, mz] (无量纲)
%   moment_mag    - 磁偶极矩标量大小 (A·m²)
%   mu_0          - (可选) 真空磁导率 (T·m/A)，默认4π×10⁻⁷
%
% 输出参数：
%   B      - 磁场强度矢量 [Bx, By, Bz] (T)
%   gradB  - 磁场梯度张量 (3×3) (T/m)
%
% 公式依据：
%   B = (μ₀/4π) * [3(m·r)r/r⁵ - m/r³]
%   ∇B = (3μ₀/4πr⁵)[(m·r)I + r⊗m + m⊗r - (5(m·r)/r²)(r⊗r)]

% 输入验证与默认值处理
arguments
    position_vec (3,1) double
    moment_unit (3,1) double
    moment_mag (1,1) double
    mu_0 (1,1) double = 4*pi*1e-7  % 默认真空磁导率
end

% 检查磁矩方向是否为矢量
if abs(norm(moment_unit) - 1) > 1e-6
    warning('单位磁矩方向矢量未归一化，将自动归一化');
    moment_unit = moment_unit / norm(moment_unit);
end

% 计算磁矩矢量
moment_vec = moment_unit * moment_mag;

% 计算位置向量的模和平方模
r = position_vec;
r_norm = norm(r);
r_norm_sq = r_norm^2;

% 检查观测点距离
if r_norm < 1e-12  % 1皮米级阈值
    error('观测点距离偶极子太近（r < 1 pm），模型不适用');
end

% 计算中间变量
m_dot_r = dot(moment_vec, r);
common_factor = mu_0/(4*pi);  % 前置比例系数

% ===== 1. 计算磁场 =====
inv_r3 = 1/(r_norm^3);
inv_r5 = 1/(r_norm^5);

B = common_factor * ( (3 * m_dot_r) * r * inv_r5 - moment_vec * inv_r3 );

% ===== 2. 计算磁场梯度 =====
if nargout > 1
    % 预计算复用项
    factor = (3 * common_factor) / (r_norm^5);
    term_5 = 5 * m_dot_r / r_norm_sq;

    % 直接计算各分量（避免中间矩阵分配）
    gradB = zeros(3, 3);

    for i = 1:3
        for j = 1:3
            % 张量分量公式
            gradB(i, j) = factor * ( ...
                m_dot_r * (i == j) ...  % (m·r)I δ_ij
                + r(i)*moment_vec(j) ... % (r⊗m)
                + moment_vec(i)*r(j) ... % (m⊗r)
                - term_5 * r(i)*r(j) ... % -(5(m·r)/r²)(r⊗r)
            );
        end
    end

    % 替代方法：使用MATLAB矩阵运算
    % identity = eye(3);
    % outer_rr = r' * r;
    % outer_rm = r' * moment_vec;
    % outer_mr = moment_vec' * r;
    % gradB = factor * (m_dot_r*identity + outer_rm + outer_mr - term_5*outer_rr);
end
end

