function [Signal,Scores]=IC_signs(Signal,Scores,Options);
% Adjust signs of loadings and scores
% USAGE :
% [Signal,Scores]=IC_signs(Signal,Scores,Options);
%
% INPUT variables :
%   Signal = IC signals ('Loadings')
%   Scores = IC proportions ('Scores')
%   Options.Criterion = 1 or 2 % Max intensity (1) / Sum surface (2)
%   Options.Basis = 1 or 2 % Scores (1) Signal (2)
%
% OUTPUT variables :
%   Signal = IC signals ('Loadings')
%   Scores = IC proportions ('Scores')


if exist('Options','var')
    if isfield(Options,'Criterion')
        Criterion=Options.Criterion;
    else
        Criterion=1;
    end
    
    if isfield(Options,'Basis')
        Basis=Options.Basis;
    else
        Basis=1;
    end
end

ICs=size(Signal,2);

for nIC=1:ICs
    switch Criterion
        case 1 % Max Intensity
            switch Basis
                case 1 % Scores
                    [maxVal,maxInd]=max(abs(Scores(:,nIC)));
                    if Scores(maxInd,nIC) < 0
                        Scores(:,nIC)=-1*Scores(:,nIC);
                        Signal(:,nIC)=-1*Signal(:,nIC);
                    end
                case 2 % Signals
                    [maxVal,maxInd]=max(abs(Signal(:,nIC)));
                    if Signal(maxInd,nIC) < 0
                        Signal(:,nIC)=-1*Signal(:,nIC);
                        Scores(:,nIC)=-1*Scores(:,nIC);
                    end
            end
        case 2 % Surface
            switch Basis
                case 1 % Scores
                    SumVals=sum(Scores(:,nIC));
                    if SumVals < 0
                        Scores(:,nIC)=-1*Scores(:,nIC);
                        Signal(:,nIC)=-1*Signal(:,nIC);
                    end
                case 2 % Signals
                    SumVals=sum(Signal(:,nIC));
                    if SumVals < 0
                        Scores(:,nIC)=-1*Scores(:,nIC);
                        Signal(:,nIC)=-1*Signal(:,nIC);
                    end
            end
    end
    
end
