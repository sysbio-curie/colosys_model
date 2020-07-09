function full_filename_with_path=fcn_save_fig(file_name_prefix,save_folder,fig_file_extension,overwrite_flag,resolution)

full_filename_with_path=strcat(save_folder,file_name_prefix,fig_file_extension);

if isempty(overwrite_flag)
    if exist(full_filename_with_path,'file')==0
        if isempty(resolution)
            export_fig(full_filename_with_path,'-transparent','-nocrop')
        else
            export_fig(full_filename_with_path,'-transparent','-nocrop',resolution)
        end
    else
        disp('Figure already exists, give different name or set overwrite flag to non-empty.')
    end

else
        if isempty(resolution)
            export_fig(full_filename_with_path,'-transparent','-nocrop')
        else
            export_fig(full_filename_with_path,'-transparent','-nocrop',resolution)
        end
        disp('existing figure file overwritten')
end
