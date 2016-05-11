addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\DeepLearnToolbox-master_modification\'));
addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\matlab_original_code_DONTCHANGE\'));
%addpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\matlab_original_code_original\');
addpath(genpath('C:\Users\Tamme\Desktop\Yli\Magistritöö\Coding\'));


%inputFileNames = {...
 %   'synthetic_5sigs_96mut_500samples_2000tot_2000samp'...
 %   'synthetic_5sigs_96mut_500samples_2000tot_200samp'...
 %   'synthetic_5sigs_96mut_500samples_2000tot_20samp'};
inputFileNames = {...
    'synthetic_5sigs_96mut_500samples_2000tot_2000samp'};
for index = 1:length(inputFileNames)
    data = load(strcat(inputFileNames{index}, '.mat'));

    [genomes mutationTypesToRemoveSet] = removeWeak(data.originalGenomes, 0.01);
    train_x = genomes';
    totalMutationTypes = size(train_x, 2);
    totalGenomes = size(train_x, 1);
    nrOfTimes = 40;
    
    opts.batchsize = size(train_x, 1);
    %opts.batchsize = 200;
    opts.print = 1;
    for nrOfSignatures = 3:4
        opts.numepochs =   15000;
        disp(nrOfSignatures);
        architecture = [totalMutationTypes nrOfSignatures]; 
        
        genomeErrors = zeros(totalMutationTypes, totalGenomes, nrOfTimes);
        genomesReconstructed = zeros(totalMutationTypes, totalGenomes, nrOfTimes);
    
        Wall = zeros(totalMutationTypes, nrOfTimes * nrOfSignatures);
        for z = 1:nrOfTimes
            disp(z)
            sae = saesetup(architecture);
            tmp = max(bootstrapCancerGenomes( genomes ), eps);
            train_x_boot = tmp';

            rand('state',z);
            
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

            startAll = (z-1)*nrOfSignatures + 1;
            endAll = startAll + nrOfSignatures - 1;
            for j = 1 : nrOfSignatures
                total = sum( W(:, j) );
                W(:, j) = W(:, j) / total;
                H(j, :) = H(j, :) * total;
            end
            Wall(:, startAll:endAll) = W;
            
            genomeErrors(:, :, z) = train_x_boot' -  W * H;
            genomesReconstructed(:, :, z) = W * H;
        end
        %keyboard;
        [ Wall, genomeErrors, genomesReconstructed ] = ...
        myFilterOut( Wall, genomeErrors, nrOfSignatures, ...
                           genomesReconstructed, 0.07 );
        %keyboard;  
        [centroids, centroidStd, idx, processStab processStabAvg clusterCompactness] = ...
                                myEvaluateStability( Wall, nrOfSignatures, nrOfTimes);
        %keyboard;
        

        if size(architecture, 2) == 3
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
            W = addWeakMutations(mutationTypesToRemoveSet, W);
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
                inputFileNames{index}, '_WALL_lotsforall_toy_AE_sigs_', ...
            num2str(nrOfSignatures));
            save(destination, 'R', 'W', 'W_first', 'H', 'allDist', 'sigCorr', 'Wall', 'processStabAvg', 'idx', 'centroids');


            destination = strcat('C:\Users\Tamme\Desktop\Yli\Magistritöö\Autoencoder\ae_testihng119_lr', ...
            num2str(sae.ae{1}.learningRate), '_epoch_', num2str(opts.numepochs), '_batch_', ...
            num2str(opts.batchsize), '_decay_', num2str(sae.ae{1}.decayRate), '_rmsprop_', ...
            num2str(sae.ae{1}.rmsProp), '_slr_', ...
            num2str(sae.ae{1}.scaling_learningRate), '.mat');
            %save(destination, 'R', 'W', 'H', 'W_first');
        end
    end
end
