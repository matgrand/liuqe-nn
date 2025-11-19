% this script converts all json files in the current directory to mat files
% usage: json2mat

curr_dir = fileparts(mfilename('fullpath'));
json_files = dir(fullfile(curr_dir, '*.json'));

for k = 1:length(json_files)
    json_file = fullfile(curr_dir, json_files(k).name);
    mat_file = strrep(json_file, '.json', '.mat');
    
    fprintf('Converting %s to %s...\n', json_file, mat_file);
    
    % Read JSON file
    fid = fopen(json_file, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    
    % Decode JSON
    data = jsondecode(str);
    
    % Convert all numeric fields to single precision
    fields = fieldnames(data);
    for i = 1:numel(fields)
        if isnumeric(data.(fields{i}))
            data.(fields{i}) = single(data.(fields{i}));
            % fix dimensions: if ends with _bias or _beta make it (1, N)
            if endsWith(fields{i}, '_bias') || endsWith(fields{i}, '_beta')
                data.(fields{i}) = reshape(data.(fields{i}), 1, []);
            end
        end
    end
    
    % Save to MAT file
    save(mat_file, '-struct', 'data');
    
    fprintf('Saved %s\n', mat_file);
end