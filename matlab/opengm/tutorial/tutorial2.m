clear all
close all
clc

filename = 'result.h5';

state = marray_load(filename, 'state');
runtimes = marray_load(filename, 'runtimes') / 1000; % seconds
values = marray_load(filename, 'values');
distances = marray_load(filename, 'distances');

fprintf('Energy E=%f at iteration %d, after t=%fs.\n', values(end), length(values), runtimes(end));

figure(1);
plot(runtimes, values);
title('Energy');
xlabel('Runtime [s]');
ylabel('Energy');

figure(2);
plot(runtimes, distances);
title('Convergence');
xlabel('Runtime [s]');
ylabel('Message Norm');

figure(3);
imagesc(reshape(state, 10, 10));
colormap gray;
axis equal tight;

