function [hasGPU, g] = if_has_gpu()
    % HASGPUFORGPUARRAY  Check if there is a usable NVIDIA GPU for gpuArray.
    %
    %   hasGPU = HASGPUFORGPUARRAY() returns true if a CUDA-capable GPU is
    %   available and gpuArray can run on it, false otherwise.
    %
    %   [hasGPU, g] = HASGPUFORGPUARRAY() also returns the gpuDevice object g
    %   if available; otherwise g is [].

    hasGPU = false;
    g = [];

    % 1) Check how many GPUs MATLAB sees
    try
        n = gpuDeviceCount("available");
    catch
        % Older MATLAB versions don’t support the "available" option
        n = gpuDeviceCount;
    end

    if n <= 0
        return;  % no GPU visible
    end

    % 2) Try to select the default GPU and run a tiny gpuArray test
    try
        g = gpuDevice();          % select default GPU
        xg = gpuArray(1);         % simple gpuArray allocation
        yg = xg + 1;              % simple operation
        gather(yg);               % pull back to CPU
        hasGPU = true;            % all good
    catch
        % Something went wrong; leave hasGPU = false and g = []
        g = [];
    end
end
