function buffer_data = set_buffer(tag, buffer_data, iter_data, iter_data_next, log_data)
    if (tag == "reset")
        buffer_data.Z_diff_norm_buffer = [];
        buffer_data.Z_diff_ang_buffer = [];
        buffer_data.Z_buffer = [];
        buffer_data.Z_diff_buffer = [];
        buffer_data.pobj_buffer = [];
        buffer_data.dobj_buffer = [];
        buffer_data.current_size = 0;
    elseif (tag == "append")
        Z_diff = iter_data_next.Z - iter_data.Z;
        if buffer_data.current_size < buffer_data.max_size
            buffer_data.Z_diff_norm_buffer = [buffer_data.Z_diff_norm_buffer, log_data.Z_diff_norm];
            buffer_data.Z_diff_ang_buffer = [buffer_data.Z_diff_ang_buffer, log_data.Z_diff_ang];
            buffer_data.Z_buffer = [buffer_data.Z_buffer, iter_data.Z];
            buffer_data.Z_diff_buffer = [buffer_data.Z_diff_buffer, Z_diff];
            buffer_data.pobj_buffer = [buffer_data.pobj_buffer, log_data.pobj];
            buffer_data.dobj_buffer = [buffer_data.dobj_buffer, log_data.dobj];
            buffer_data.current_size = buffer_data.current_size + 1;
        else
            buffer_data.Z_diff_norm_buffer = [buffer_data.Z_diff_norm_buffer(2: end), log_data.Z_diff_norm];
            buffer_data.Z_diff_ang_buffer = [buffer_data.Z_diff_ang_buffer(2: end), log_data.Z_diff_ang];
            buffer_data.Z_buffer = [buffer_data.Z_buffer(:, 2: end), iter_data.Z];
            buffer_data.Z_diff_buffer = [buffer_data.Z_diff_buffer(:, 2: end), Z_diff];
            buffer_data.pobj_buffer = [buffer_data.pobj_buffer(2: end), log_data.pobj];
            buffer_data.dobj_buffer = [buffer_data.dobj_buffer(2: end), log_data.dobj];
        end
    else 
        error("Tag undefined! \n");
    end
end