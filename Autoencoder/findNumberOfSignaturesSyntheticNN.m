%path to deeplear toolbox
addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\DeepLearnToolbox-master_modification\'));
%path to pour code
addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\Coding\'));


inputFileNames = {...
%    'synthetic_5sigs_96mut_500samples_2000tot_2000samp'...
%    'synthetic_5sigs_96mut_500samples_2000tot_200samp'...
    'synthetic_5sigs_96mut_500samples_2000tot_20samp'};
%inputFileNames = {...
%    'EasierHsynthetic_5sigs_96mut_500samples_2000tot_2000samp'};
for index = 1:length(inputFileNames)
    data = load(strcat(inputFileNames{index}, '.mat'));

    [genomes mutationTypesToRemoveSet] = removeWeak(data.originalGenomes, 0.01);
    train_x = genomes';
    totalMutationTypes = size(train_x, 2);
    totalGenomes = size(train_x, 1);
    nrOfTimes = 100;
    opts.numepochs =   20000;
    opts.batchsize = size(train_x, 1);
    %if want to print continuous output of ae
    opts.print = 1;
    for nrOfSignatures = 3:7

        disp(nrOfSignatures);
        architecture = [totalMutationTypes nrOfSignatures]; 

        sae = saesetup(architecture);
        for z = 1:nrOfTimes
            disp(z);
            tmp = max(bootstrapCancerGenomes( genomes ), eps);
            train_x_boot = tmp';

            rand('state',z);
            if z > 1
                opts.numepochs =   6000;
            end

            sae.ae{1}.activation_function       = 'relu';
            sae.ae{1}.output       = 'linear';
            sae.ae{1}.momentum       = 0.9;
            sae.ae{1}.learningRate              = 0.0001;
            sae.ae{1}.rmsProp              = 1;
            sae.ae{1}.decayRate              = 0.999;
            sae.ae{1}.inputZeroMaskedFraction   = 0.00  ;
            sae.ae{1}.scaling_learningRate      = 0.9999;
            %sae.ae{1}.dropoutFraction = 0.01;

            if size(architecture, 2) == 3
                sae.ae{2}.activation_function       = 'relu';
                sae.ae{2}.output       = 'linear';
                sae.ae{2}.momentum       = 0.9;
                sae.ae{2}.learningRate              = 0.00001;
                sae.ae{2}.rmsProp              = 1;
                sae.ae{2}.decayRate              = 0.9;
                sae.ae{2}.inputZeroMaskedFraction   = 0;
                sae.ae{2}.scaling_learningRate      = 1;
            end

            %opts.plot = 1;

            sae = saetrain(sae, train_x_boot, opts);
            W_first= sae.ae{1}.W{1};
            H = max(0,W_first * train_x_boot');
            W = sae.ae{1}.W{2};

            W = addWeakMutations(mutationTypesToRemoveSet, W);
            R =  max(0, W * H);
            disp(norm(R - data.originalGenomes, 'fro'));
        end
        
        if size(architecture, 2) == 3
            %for stacked ae. didn't analyse too much.
            %stackedreconstruction
            %W1= sae.ae{1}.W{1};
            H1 = max(0,train_x * sae.ae{1}.W{1}');
            H2 =  max(0,H1 * sae.ae{2}.W{1}');
            H3 =  max(0,H2 * sae.ae{2}.W{2}');
            R =  max(0,H3 * sae.ae{1}.W{2}');
            W10 = sae.ae{1}.W{2};
            W5 = sae.ae{2}.W{2};
            H10 = H3;
            H5 = H2;
            disp(norm(R-train_x, 'fro'));
            
            save('C:\Users\Tamme\Desktop\Yli\Magistritöö\Autoencoder\ae_stacked.mat', 'R', 'W10', 'W5', 'H10', 'H5');
            keyboard;
        else
            %W = addWeakMutations(mutationTypesToRemoveSet, W);
            R =  max(0, W * H);
            disp('main');
            disp(norm(R - data.originalGenomes, 'fro'));
            allDist = squareform( pdist(W', 'cosine') );
            sigCorr = corr(W);
            N = nrOfSignatures + 1;
            foundSimilar = 0;
            if max(max(corr(W) - eye(size(W, 2)))) > 0.7
                disp('Found similar signatures using correlation');
            end
            if max(max(1 - squareform(pdist(W', 'cosine'))- eye(size(W, 2)))) > 0.7 
                %and spars below 0.3 or sth
                disp('Found similar signatures using cosine similarity');
            end

            W_first= sae.ae{1}.W{1};
            sum(sum(W == 0))/(96*size(W, 2))
            destination = strcat('C:\Users\Tamme\Desktop\Yli\Magistritöö\Coding\Autoencoder\models\', ...
                inputFileNames{index}, '_2second_toy_AE_sigs_', ...
            num2str(nrOfSignatures));
            save(destination, 'R', 'W', 'W_first', 'H', 'allDist', 'sigCorr');

            destination = strcat('C:\Users\Tamme\Desktop\Yli\Magistritöö\Autoencoder\ae_testihng119_lr', ...
            num2str(sae.ae{1}.learningRate), '_epoch_', num2str(opts.numepochs), '_batch_', ...
            num2str(opts.batchsize), '_decay_', num2str(sae.ae{1}.decayRate), '_rmsprop_', ...
            num2str(sae.ae{1}.rmsProp), '_slr_', ...
            num2str(sae.ae{1}.scaling_learningRate), '.mat');
            %save(destination, 'R', 'W', 'H', 'W_first');
        end
    end
end
