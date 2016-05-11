%bayes opt impl

F = @(params) saewrapper(params(1), params(2), params(3), params(4), params(5));
addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\DeepLearnToolbox-master_modification'));
startup;

% params = [lr, decay for rmsprop, scaling_lr, epochs]
opt = defaultopt; 
opt.dims = 5;
opt.mins = [5, 0.5, 0.5, 2, 0];
opt.maxes = [0, 1, 1, 6, 0.999];
opt.max_iters = 300;
opt.grid_size = 6*6*6*5*5;

lr = [0, 1, 2, 3, 4, 5];
dr = [0.5, 0.9, 0.99, 0.999, 0.9999, 1];
slr = [0.5, 0.9, 0.99, 0.999, 0.9999, 1];
epo = [2, 3, 4, 5, 6];
mom = [0, 0.5, 0.9, 0.99, 0.999];


gred = zeros(6*6*6*5*5, opt.dims);
c = 1;
for i = 1 : 6
    for j = 1: 6
        for k = 1: 6
            for l = 1 : 5
                for m = 1 : 5
                    gred(c, :) = [lr(i), dr(j), slr(k), epo(l), mom(m)];
                    c = c + 1;
                end
            end
        end
    end
end

gred = gred(randperm(size(gred,1)),:);
opt.grid = gred;

[ min_sample , min_value , botrace ] = bayesopt(F, opt);
disp('over');