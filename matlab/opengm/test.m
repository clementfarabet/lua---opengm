clear all
close all
clc


%% build graphical model
numbers_of_states = [3 2 2 4];
gm = opengm_create(numbers_of_states);

v = cell(3, 1);
t = cell(3, 1);

v{1} = 1;
t{1} = [1; 2; 3];
opengm_add_factor(gm, v{1}, t{1});

v{2} = [2 3];
t{2} = [1 2; 3 4];
opengm_add_factor(gm, v{2}, t{2});

v{3} = [2 3 4];
t{3} = rand(2, 2, 4);
opengm_add_factor(gm, v{3}, t{3});


%% test build
myassert(opengm_number_of_variables(gm) == 4);
for j = 1:opengm_number_of_variables(gm)
    myassert(gm.numbers_of_states(j) == numbers_of_states(j));
end

myassert(opengm_number_of_factors(gm) == 3);
for j = 1:opengm_number_of_factors(gm)
    [vi, tab] = opengm_factor(gm, j);
    myassert(all( vi == v{j}(:) ));
    myassert(all( tab == t{j}   ));
end


%% edit graphical model
vn = [1 4];
tn = [1 2 3 4; 5 6 7 8; 9 10 11 12];
opengm_update_factor(gm, 2, vn, tn);


%% test edit
myassert(opengm_number_of_factors(gm) == 3);
[vi, tab] = opengm_factor(gm, 2);
myassert(all( vi == vn(:) ));
myassert(all( tab == tn   ));


%% serialize and save
ser = opengm_serialize(gm);
marray_save('test-gm.h5', 'graphical-model', ser);


%% de-serialize
gm2 = opengm_deserialize(ser);


%% test de-serialization
ser2 = opengm_serialize(gm2);
myassert(all(ser == ser2));

