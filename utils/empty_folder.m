function empty_folder(folderPath)
    % EMPTYFOLDER  Delete all files and subfolders inside folderPath.
    %   EMPTYFOLDER(folderPath) removes all contents of the folder, but keeps
    %   the folder itself.

    if nargin < 1 || ~isfolder(folderPath)
        error('emptyFolder:InvalidPath', 'Folder does not exist: %s', folderPath);
    end

    % List everything in the folder
    contents = dir(folderPath);

    % Remove '.' and '..'
    contents(ismember({contents.name}, {'.', '..'})) = [];

    for k = 1:numel(contents)
        itemPath = fullfile(folderPath, contents(k).name);

        if contents(k).isdir
            % Delete subfolder and its contents
            rmdir(itemPath, 's');
        else
            % Delete file
            delete(itemPath);
        end
    end
end