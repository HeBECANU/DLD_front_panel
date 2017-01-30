function = wait_till_written(watching_dir,new_file,wait_for_mod)

    while is_early(monitor_dir,new_files{k},wait_for_mod)
                    fprintf('\b\b\b')
                    pause(0.05*wait_for_mod)
                    fprintf('.')
                    pause(0.05*wait_for_mod)
                    fprintf('.')
                    pause(0.05*wait_for_mod)
                    fprintf('.')
                    pause(0.05*wait_for_mod)
    end

end




function[out]= is_early(watching_dir,new_file,wait_for_mod)
%this function simply checks that the modification date of a file is
%wait_for_mod seconds in the past
    file_pointer=dir(fullfile(watching_dir,new_file));
    out=addtodate(datenum(file_pointer.date),wait_for_mod,'second')>now;
end