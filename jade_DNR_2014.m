function [Signal,Scores] =  jade_DNR_2014(X, ICs, Options);

% INPUT variables :
%
%   X(nxT) = mixed signal dataset
%   ICs=  Sources (pure signal components) =underlying signals =pure spectra
%   where :
%   n=  Sensors (mixed signal inputs) = rows =spectra
%   T=  Signal data points = columns =variables
%
% Options.Method = 'Tall_Wide', 'Tall', 'Wide', 'Kernel', 'ComDim', 'Normal'
% Options.Partitions = number of partition if Options.Method = 'Tall' or 'Wide'
% Options.rows_par=number of row partition if Options.Method = 'Tall_Wide'
% Options.cols_par=number of col partition if Options.Method = 'Tall_Wide'
% Options.Data= 'Data', 'Loadings'
% Options.SortICs =0/1; % Sort ICs by Scores Variance or not (Default=0 : Not)
%
%   Options.Centred = Matrix centred by column
%
% OUTPUT variables :
%
%   Signal(nxT) or Signal(ICsxT)  = pure signal dataset
%   Scores = "scores" of samples calculated from ICs
%          = old_X*Signal*inv(Signal'*Signal);
%
%
% CALLS :
%
% jadeR_2005
% ColumnsPartition
% PCA_Tall_PCT_DNR
% PCA_TW_DNR
% ColumnsPartition
% RowsPartition
% IC_sort

[nR,nC]	= size(X);
MaxRank=min(nR,nC);

old_X=X;

if exist('Options','var')
    if isfield(Options,'SortICs')
        SortICs= Options.SortICs;
    else
        SortICs=0;
    end
    
    if ~isfield(Options,'Data')
        Options.Data='Data';
    end

    switch Options.Data
        case 'Loadings'

            Options.Method=('Normal');
            B =  jadeR_2005(X,ICs);

            % Calculate pure signals from unmixing matrix and X
            Signal=(B*old_X)';
%             Signal=(B*X)';

        otherwise
            if isfield(Options,'Method');
                switch Options.Method
                    case 'Tall'
                        if isfield(Options,'Partitions');
                            partitions=Options.Partitions;
                        else
                            Options.Partitions=1;
                            partitions=Options.Partitions;
                        end

%##########################
% %                         Xi = ColumnsPartition(X', partitions);
% % %                         T_all=zeros(size(Xi,1),partitions);
% %                         T_all=[];
% %                         for p = 1:partitions
% %                             [U_in,S_in,V_in]= svd(Xi{p}', 'econ');
% %                             T_in=X*V_in;
% %                             T_all=[T_all T_in];
% % %                             T_all(:,p)= T_in;
% %                         end
% %                         %% Concatenated Scores of Segments
% %                         [u,s,v]= svd(T_all, 'econ');
%##########################
                        
%##########################
                        %% Do Tall segmented PCT
                        [t] = PCA_Tall_PCT_DNR(X, MaxRank, partitions);

                        [u,s,v]=svd(t, 'econ');
%##########################

                        % Calculate PCT Loadings
%                         Techelle=inv(u'*u)*u';
                        Techelle=(u'*u)\u';
                        Vpct=Techelle*old_X;

                        B =  jadeR_2005(Vpct,ICs);

                        % Calculate pure signals from unmixing matrix and X
                        Signal=(B*Vpct)';

                    case 'Wide'
                        if isfield(Options,'Partitions');
                            partitions=Options.Partitions;
                        else
                            Options.Partitions=1;
                            partitions=Options.Partitions;
                        end

                        Xi = ColumnsPartition(X, partitions);
%                         T_all=zeros(size(Xi,1),partitions);
                        T_all=[];
                        for p = 1:partitions
                            [U_in,S_in,V_in]= svd(Xi{p}, 'econ');
                            minICS=min(ICs,size(U_in,2));
                            T_all=[T_all U_in];
%                             T_all(:,p)= U_in;
                        end
                        %% Do PCT on concatenated Scores of Segments
                        [u,s,v]= svd(T_all, 'econ');

                        %% Calculate PCT Loadings
%                         Techelle=inv(u'*u)*u';
                        Techelle=(u'*u)\u';
                        Vpct=Techelle*old_X;

                        B =  jadeR_2005(Vpct,ICs);

                        % Calculate pure signals from unmixing matrix and X
                        Signal=(B*Vpct)';

                    case 'Kernel'
                        if nR>nC
                            Kernel=(X'*X)/nC;

                            [u,s,v]= svd(Kernel,'econ')	;
                            D=s*s;
                            % D= Eigenvalues

                            Ut=u; % Ut here == v from SVD on X
                            % Ut= Eigenvectors

                            [puiss,k]= sort(diag(D),'descend')	;
                            rangeW= 1:ICs			;   % indices to the ICs  most significant directions

                            Ut=Ut(:,k(rangeW));

                            W=(X*Ut)';% whitener
                            X= W*old_X;

                            B =  jadeR_2005(X,ICs);

                            % Calculate pure signals from unmixing matrix and X
                            Signal=(B*X)';
                        else
                            Kernel=(X*X')/nR;

                            [u,s,v]= svd(Kernel,'econ')	;
                            D=s*s;
                            % D= Eigenvalues

                            Ut=u; % Ut here == v from SVD on X
                            % Ut= Eigenvectors

                            [puiss,k]= sort(diag(D),'descend')	;
                            rangeW= 1:ICs;   % indices to the ICs  most significant directions

                            Ut=Ut(:,k(rangeW));

                            %% Calculate PCT Loadings
%                             Techelle=inv(Ut'*Ut)*Ut';
                            Techelle=(Ut'*Ut)\Ut';
                            Vpct=Techelle*old_X;

                            B =  jadeR_2005(Vpct,ICs);

                            % Calculate pure signals from unmixing matrix and X
                            Signal=(B*Vpct)';
                        end

                    case 'Tall_Wide'
                        % Options.rows_par=number of row partition if Options.Method = 'Tall_Wide'
                        % Options.cols_par=number of col partition if Options.Method = 'Tall_Wide'

                        if isfield(Options,'rows_par')
                            rows_par=Options.rows_par;
                        else
                            Options.rows_par=1;
                            partitions=Options.rows_par;
                        end

                        if isfield(Options,'cols_par')
                            cols_par=Options.cols_par;
                        else
                            Options.cols_par=1;
                            partitions=Options.cols_par;
                        end

                        %
                        PCs=ICs;
                        fullRank=ICs;
                        [T, Vq] = PCA_TW_DNR(X, PCs, rows_par, cols_par, fullRank);
                        
                        % Verify NaNs & Inf in Vq
                        B =  jadeR_2005(Vq',ICs);

                        % Calculate pure signals from unmixing matrix and X
                        Signal=(B*Vq')';
                        %

                    case {'Normal' 'PCT'}
                        B =  jadeR_2005(X,ICs);

                        % Calculate pure signals from unmixing matrix and X
                        Signal=(B*old_X)';

                end
            else
                Options.Method='Normal';
                B =  jadeR_2005(X,ICs);

                % Calculate pure signals from unmixing matrix and X
                Signal=(B*old_X)';

            end

    end
end

% Calculate IC_Scores (proportions) of Jade
% Scores=old_X*Signal*inv(Signal'*Signal);
Scores=old_X*Signal/(Signal'*Signal);

% Sort ICs by Variance of Scores
if SortICs~=0
    [Signal, Scores, Scores_Variance] = IC_sort(Signal,Scores);
end


%% For information
% Scores=old_X*Signal;
% Scores=X*Signal*inv(Signal'*Signal);

