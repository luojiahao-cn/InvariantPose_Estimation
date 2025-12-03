function inspect_slanCM_data()
% INSPECT_SLANCM_DATA 检查当前工作空间中的 slanCM_Data 变量
% 
% 使用说明：
%   1. 先运行: load('D:\Softwares\MATLAB\toolbox\slan_CM\slanCM_Data.mat')
%   2. 然后运行此函数查看变量结构

fprintf('=== 检查 slanCM_Data 变量 ===\n\n');

%% 检查变量是否存在
if ~exist('slanCM_Data', 'var')
    fprintf('错误: 工作空间中不存在 slanCM_Data 变量\n');
    fprintf('请先运行: load(''D:\\Softwares\\MATLAB\\toolbox\\slan_CM\\slanCM_Data.mat'')\n');
    return;
end

%% 显示基本信息
fprintf('变量名: slanCM_Data\n');
fprintf('类型: %s\n', class(slanCM_Data));
fprintf('大小: %s\n', mat2str(size(slanCM_Data)));

%% 根据类型显示详细信息
if isstruct(slanCM_Data)
    fprintf('\n结构体字段 (%d个):\n', numel(fieldnames(slanCM_Data)));
    fprintf('----------------------------------------\n');
    
    field_names = fieldnames(slanCM_Data);
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = slanCM_Data.(field_name);
        
        fprintf('\n[%d] %s\n', i, field_name);
        fprintf('    类型: %s\n', class(field_data));
        
        if isnumeric(field_data)
            fprintf('    大小: %s\n', mat2str(size(field_data)));
            if size(field_data, 2) == 3 && size(field_data, 1) > 1
                fprintf('    格式: Nx3 RGB颜色映射\n');
                fprintf('    点数: %d\n', size(field_data, 1));
                fprintf('    起始RGB: [%.4f, %.4f, %.4f]\n', ...
                    field_data(1, 1), field_data(1, 2), field_data(1, 3));
                fprintf('    结束RGB: [%.4f, %.4f, %.4f]\n', ...
                    field_data(end, 1), field_data(end, 2), field_data(end, 3));
            end
        elseif isstruct(field_data)
            fprintf('    子字段数: %d\n', numel(fieldnames(field_data)));
        elseif iscell(field_data)
            fprintf('    单元数组大小: %s\n', mat2str(size(field_data)));
        end
    end
    
elseif isa(slanCM_Data, 'containers.Map')
    fprintf('\nMap容器:\n');
    fprintf('----------------------------------------\n');
    keys_list = keys(slanCM_Data);
    fprintf('键数量: %d\n', length(keys_list));
    fprintf('\n键列表:\n');
    for i = 1:length(keys_list)
        key = keys_list{i};
        value = slanCM_Data(key);
        fprintf('  [%d] %s\n', i, key);
        if isnumeric(value) && size(value, 2) == 3
            fprintf('       类型: Nx3 RGB颜色映射\n');
            fprintf('       大小: %dx3\n', size(value, 1));
        end
    end
    
elseif isnumeric(slanCM_Data)
    fprintf('\n数值数组:\n');
    fprintf('----------------------------------------\n');
    if size(slanCM_Data, 2) == 3
        fprintf('可能是颜色映射 (Nx3 RGB数组)\n');
        fprintf('颜色点数: %d\n', size(slanCM_Data, 1));
        fprintf('前5行:\n');
        disp(slanCM_Data(1:min(5, size(slanCM_Data, 1)), :));
    end
    
else
    fprintf('\n其他类型，详细信息:\n');
    fprintf('----------------------------------------\n');
    disp(slanCM_Data);
end

fprintf('\n=== 检查完成 ===\n');

%% 提供添加颜色的建议
fprintf('\n提示: 如果您想添加新颜色，可以使用以下方式：\n');
if isstruct(slanCM_Data)
    fprintf('  slanCM_Data.新颜色名 = 您的颜色映射数组;\n');
    fprintf('  然后保存: save(''D:\\Softwares\\MATLAB\\toolbox\\slan_CM\\slanCM_Data.mat'', ''slanCM_Data'')\n');
end

end

