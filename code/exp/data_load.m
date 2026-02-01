clc; clear; close all;

%% Environment setup
this_file = mfilename('fullpath');
exp_dir = fileparts(this_file);
repo_root = fileparts(exp_dir);
addpath(fullfile(repo_root, 'utils'));
addpath(fullfile(repo_root, 'Functions'));
addpath(fullfile(repo_root, 'tools'));

%% Locate scan_results JSON file
% scan_filename = 'scan_results_20260103_145914.json'; % Specify filename or leave empty to auto-select
scan_filename = './data/scan_results_rotation.json';
data_pattern = 'scan_results_*.json';
if ~(exist('scan_filename', 'var') && ischar(scan_filename) && ~isempty(scan_filename))
	scan_filename = select_latest_scan(exp_dir, data_pattern);
end
scan_filepath = fullfile(exp_dir, scan_filename);

if ~isfile(scan_filepath)
	error('Scan file not found: %s', scan_filepath);
end

fprintf('Loading scan results from %s\n', scan_filepath);
scan_records = parse_scan_results(scan_filepath);

num_records = numel(scan_records);
num_sensors = unique([scan_records.num_sensors]);
fprintf('Loaded %d records; sensors per record: %s\n', num_records, mat2str(num_sensors));

overview = table(...
	[scan_records.timestamp]', ...
	string({scan_records.target_id})', ...
	[scan_records.num_samples]', ...
	[scan_records.num_sensors]', ...
	'VariableNames', {'timestamp', 'target_id', 'num_samples', 'num_sensors'});

head_count = min(5, height(overview));
disp('Preview of decoded records:');
disp(overview(1:head_count, :));

save(fullfile(exp_dir, 'scan_records_rotation.mat'), 'scan_records');
%% Helper functions
function latest_name = select_latest_scan(data_dir, pattern)
	files = dir(fullfile(data_dir, pattern));
	if isempty(files)
		error('No files matching pattern %s in %s', pattern, data_dir);
	end
	[~, idx] = max([files.datenum]);
	latest_name = files(idx).name;
end

function records = parse_scan_results(file_path)
	raw_text = fileread(file_path);
	json_data = jsondecode(raw_text);
	num_entries = numel(json_data);
	template = struct('timestamp', [], 'target_id', '', 'num_samples', [], ...
		'poses', [], 'hall_matrix', [], 'num_sensors', []);
	records = repmat(template, num_entries, 1);
	for idx = 1:num_entries
		entry = json_data(idx);
		records(idx).timestamp = entry.timestamp;
		records(idx).target_id = string(entry.target_id);
		records(idx).num_samples = entry.num_samples;
		records(idx).poses = normalize_pose_struct(entry.poses);
		records(idx).hall_matrix = hall_struct_to_matrix(entry.hall_data);
		records(idx).num_sensors = size(records(idx).hall_matrix, 2);
	end
end

function pose_struct = normalize_pose_struct(raw_poses)
	pose_names = fieldnames(raw_poses);
	pose_struct = struct();
	for i = 1:numel(pose_names)
		name = pose_names{i};
		pose_struct.(name).position = raw_poses.(name).position(:);
		pose_struct.(name).orientation = raw_poses.(name).orientation(:);
	end
end

function hall_matrix = hall_struct_to_matrix(hall_data)
	x = [hall_data.x];
	y = [hall_data.y];
	z = [hall_data.z];
	hall_matrix = [x; y; z];
end
