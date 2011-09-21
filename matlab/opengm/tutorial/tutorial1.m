clear all
close all
clc

grid_size = 10;
alpha = 0.2;
filename = 'gm.h5';

% graphical model
number_of_variables = grid_size * grid_size;
numbers_of_states = 2 * ones(number_of_variables, 1);
gm = opengm_create(numbers_of_states);

% single site factors
for j = 1:number_of_variables
    variable_indices = j;
    table = rand(2, 1);
    opengm_add_factor(gm, variable_indices, table);
end

% Potts factors
for x0 = 1:grid_size
for y0 = 1:grid_size
    if x0 < grid_size
        variable_indices = [sub2ind([grid_size grid_size], x0, y0)
                            sub2ind([grid_size grid_size], x0+1, y0)];
        table = [0 alpha;
                 alpha 0];
        opengm_add_factor(gm, variable_indices, table);
    end
    if y0 < grid_size
        variable_indices = [sub2ind([grid_size grid_size], x0, y0)
                            sub2ind([grid_size grid_size], x0, y0+1)];
        table = [0 alpha;
                 alpha 0];
        opengm_add_factor(gm, variable_indices, table);
    end   
end
end

% serialize and save graphical model
ser = opengm_serialize(gm);
marray_save(filename, 'graphical-model', ser);

