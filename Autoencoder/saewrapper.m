function L = saewrapper(lr, dr, slr, epo, mom)
disp(sprintf(['\nlr ', num2str(power(10,-1*lr)), '\ndecay ', num2str(dr), '\nscaling ', num2str(slr), '\nepochs ', num2str(epo)]));

addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\DeepLearnToolbox-master_modification'));

data = load('21_WTSI_BRCA_whole_genome_substitutions.mat');
train_x = data.originalGenomes';
%normalize
%for i = 1 : size(train_x, 1)
%    train_x(i, :) = train_x(i, :) / sum(train_x(i, :));
%end
   
%  Setup and train a stacked denoising autoencoder (SDAE)
rand('state',0);
sae = saesetup([96 5]);
sae.ae{1}.activation_function       = 'relu';
sae.ae{1}.output       = 'linear';
sae.ae{1}.momentum       = mom;
sae.ae{1}.learningRate              =  power(10,-1*lr);
sae.ae{1}.rmsProp              = 0;
sae.ae{1}.decayRate              = dr;
sae.ae{1}.inputZeroMaskedFraction   = 0;
sae.ae{1}.scaling_learningRate      = slr;
%sae.ae{1}.dropoutFraction = 0.5;
%sae.ae{1}.momentum = 0.1;
opts.numepochs =   power(10, epo);
opts.batchsize = 7;
opts.print = 0;
%opts.plot = 1;

sae = saetrain(sae, train_x, opts);
[sae2 final1] = nnloss2(sae.ae{1}, train_x, train_x);
R2 = sae2.a{3};
L = norm(R2-train_x, 'fro');
%disp(['Result', num2str(L)]);

