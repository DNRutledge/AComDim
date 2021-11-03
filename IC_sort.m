function [Signal, Scores, Scores_Variance] = IC_sort(Signal,Scores);
% Sort ICs based on decreasing Scores Variances
% USAGE
% [Signal, Scores, Scores_Variance] = IC_sort(Signal,Scores);

% INPUT :
% Signal : From ICA
% Scores : From ICA
%
% OUTPUT :
% Signal : After sorting 
% Scores : After sorting 
% Scores_Variance : Sorted Variance values for Scores vectors

Scores_Variance=var(Scores);
[Scores_Variance,Index]=sort(Scores_Variance,2,'descend');

Scores=Scores(:,Index);
Signal=Signal(:,Index);

