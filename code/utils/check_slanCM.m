% 快速查看 slanCM_Data 变量
% 运行此脚本前，请确保已加载: load('D:\Softwares\MATLAB\toolbox\slan_CM\slanCM_Data.mat')

fprintf('=== 检查 slanCM_Data 变量 ===\n\n');

% 检查变量是否存在
if ~exist('slanCM_Data', 'var')
    fprintf('错误: 工作空间中不存在 slanCM_Data 变量\n');
    fprintf('请先运行: load(''D:\\Softwares\\MATLAB\\toolbox\\slan_CM\\slanCM_Data.mat'')\n');
    return;
end

% 显示基本信息
fprintf('变量名: slanCM_Data\n');
fprintf('类型: %s\n', class(slanCM_Data));
fprintf('大小: %s\n\n', mat2str(size(slanCM_Data)));

% 如果是结构体
if isstruct(slanCM_Data)
    field_names = fieldnames(slanCM_Data);
    fprintf('结构体字段 (%d个):\n', length(field_names));
    fprintf('----------------------------------------\n');
    
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = slanCM_Data.(field_name);
        
        fprintf('\n[%d] %s\n', i, field_name);
        fprintf('    类型: %s\n', class(field_data));
        
        if isnumeric(field_data)
            fprintf('    大小: %s\n', mat2str(size(field_data)));
            if size(field_data, 2) == 3 && size(field_data, 1) > 1
                fprintf('    ✓ 这是颜色映射 (Nx3 RGB)\n');
                fprintf('    颜色点数: %d\n', size(field_data, 1));
                fprintf('    起始RGB: [%.4f, %.4f, %.4f]\n', ...
                    field_data(1, 1), field_data(1, 2), field_data(1, 3));
                fprintf('    结束RGB: [%.4f, %.4f, %.4f]\n', ...
                    field_data(end, 1), field_data(end, 2), field_data(end, 3));
            end
        end
    end
    
    fprintf('\n=== 所有颜色映射列表 ===\n');
    for i = 1:length(field_names)
        fprintf('  %d. %s\n', i, field_names{i});
    end
    
elseif isnumeric(slanCM_Data)
    fprintf('数值数组:\n');
    if size(slanCM_Data, 2) == 3
        fprintf('可能是颜色映射 (Nx3 RGB数组)\n');
        fprintf('颜色点数: %d\n', size(slanCM_Data, 1));
    end
end

fprintf('\n=== 检查完成 ===\n');

