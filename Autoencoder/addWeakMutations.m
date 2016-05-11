function [ W ] = addWeak( mutationTypesToAddSet,  W_small)

% Credit also to:
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
                                                   
   totalMutTypes = size(W_small, 1) + length(mutationTypesToAddSet);
   W = zeros( totalMutTypes, size(W_small,2) );
   
   origArrayIndex = 1;
   for i = 1 : totalMutTypes
     if ( isempty( find(mutationTypesToAddSet == i) ) == 1 )        
         W(i, :) = W_small(origArrayIndex, :);
         origArrayIndex = origArrayIndex + 1;
     end
     
   end    
end