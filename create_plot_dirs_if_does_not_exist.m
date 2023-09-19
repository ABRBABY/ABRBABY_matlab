function [] = create_plot_dirs_if_does_not_exist(plot_dir)
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

png_folder = fullfile(plot_dir,'png_folder') ;                          % path to save png files of plots
svg_folder =  strrep(png_folder,'png','svg') ;
fig_folder = strrep(png_folder,'png','fig') ;

if ~exist(svg_folder,'dir')
    mkdir(svg_folder)

end

if ~exist(png_folder,'dir')
    mkdir(png_folder)

end

if ~exist(fig_folder,'dir')
    mkdir(fig_folder)

end

