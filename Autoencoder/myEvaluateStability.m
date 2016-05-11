%% Evalute process stability
function [centroids centroidStd idx processStab processStabAvg clusterCompactness] = ... 
               myEvaluateStability( Wall, totalProcesses, nrOfIterations)
% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.

                                          
  % Defining function constants
  BIG_NUMBER = 100;
  CONVERG_ITER = 10;
  CONVERG_CUTOFF = 0.005; % cosine distance
  TOTAL_INIT_CONDITIONS = 5;
  processesDist = 'cosine';
  % Clustering mutational processes using custom clustering procedure
  minClusterDist = BIG_NUMBER;
  totalIter = size(Wall, 2) / totalProcesses;
  idx = zeros(size(Wall, 2), 1);
  %keyboard;
  clusterCompactness = zeros(totalProcesses, totalIter);
  iStartigDataSet = 1 : totalProcesses : size(Wall,2);
  iStartingDataSet = iStartigDataSet(randperm(totalIter));
  totalReplicates = 100;
  
  for iInitData = 1 : min(TOTAL_INIT_CONDITIONS, totalIter) % size(Wall, 2)
      iStartingData = iStartingDataSet(iInitData);
      iEnd = iStartingData + totalProcesses - 1;
      centroids = Wall(:, iStartingData:iEnd);
      centroidsTest = rand(size(centroids));
      countIRep = 0;
      
      for iRep = 1 : totalReplicates
          allDist = squareform( pdist(cat(2, centroids, Wall)', processesDist) );
          centroidDist = allDist( 1:size(centroids, 2), (size(centroids, 2)+1): size(allDist, 2) )';

          jRange = randperm(totalProcesses);
          for jIndex = 1 : totalProcesses
            j = jRange(jIndex);
            for i = 1 : totalProcesses : size(Wall, 2)
              iRange = i: (i + totalProcesses - 1);
                  [noUsed, Ind] = min( centroidDist(iRange, j) );
                  centroidDist(iRange(Ind), : ) = BIG_NUMBER;
                  idx( iRange(Ind) ) = j;      
             end
          end

          maxDistToNewCentroids = 0;
          for i = 1 : size(centroids, 2);
            centroids(:, i) = mean( Wall(:, idx == i), 2 );
            maxDistToNewCentroids = max(maxDistToNewCentroids, pdist( cat(2, centroids(:, i), centroidsTest(:, i))', 'cosine'));
          end

          if (  maxDistToNewCentroids < CONVERG_CUTOFF )
              countIRep =  countIRep + 1;
          else
              countIRep = 0;
              centroidsTest = centroids;
          end
          
          if (  countIRep == CONVERG_ITER )
             break;
          end

      end
      
      for i = 1 : size(centroids, 2);
          clusterDist = squareform(pdist( cat(2, centroids(:,i), Wall(:, idx ==i) )', processesDist));
          clusterCompactness(i, :) = clusterDist(1, 2:size(clusterDist, 2));
      end
     

      if ( minClusterDist > mean(clusterCompactness(:)) )
          centroidsFinal = centroids;
          idxFinal = idx;
          clusterCompactnessFinal = clusterCompactness;
      end
  end
  
  centroids = centroidsFinal';
  idx = idxFinal;
  clusterCompactness = clusterCompactnessFinal;
    
  centDist = mean( clusterCompactness, 2 );
  
  [noUsed, centDistInd] = sort(centDist, 'ascend'); % centDistSorted not used
  clusterCompactness = clusterCompactness(centDistInd,:);
  centroids = centroids(centDistInd, :);
  idxNew = idx;
  
  for i = 1 : totalProcesses
    idxNew(idx == centDistInd(i)) = i;
  end
  idx = idxNew;   
  %keyboard;
  if ( totalProcesses > 1)
      processStab = silhouette(Wall', idx, processesDist);
      processStabAvg = zeros(1, totalProcesses);
      %keyboard;
      for i = 1 : totalProcesses
        processStabAvg(i) = mean ( processStab(idx==i) ); 
      end
  else
      allDist = squareform( pdist(cat(2, centroids', Wall)', processesDist) );
      processStab = 1 - allDist( 1:size(centroids', 2), (size(centroids', 2)+1): size(allDist, 2) )';
      processStabAvg = mean(processStab);
  end
  
  centroidStd = zeros( size(centroids) );
  for i = 1 : totalProcesses
      centroidStd(i,:) = std( Wall(:, idx == i), [], 2 );
  end
  
  centroids = centroids';
  centroidStd = centroidStd';





end

