function [adecomp_res] = adecomp_DNR(X, Y, Options)
%
% The function adecomp performs the multivariate ANalysis Of VAriance
% (ANOVA) decomposition of the X matrix. X contains experiments in the rows
% and variables in the columns. It is usually the result of designed and
% organized experiments aiming to tackle specific research questions.
%
% To do so, a set of independent variables (factors) which are expected to
% contribute significantly to the experimental variability is chosen and
% levels of these factors are defined in order to perform experiments under
% controlled conditions. The results of the experiments, usually
% multivariate, are investigated to assert whether such variables have an
% effect on the system under study.
%
% The decomposition can be done in two ways :
% (1) classical multivariate ANOVA decomposition
% into matrix effects for factors/interactions by
% successively extracting the grand mean matrix, main factor effects and
% interactions based on the means of factors, levels, and interactions;
% (2) multivariate ANOVA decomposition based on
% the General Linear Models (GLM) presented by Thiel et al. (2017).
%
% Reference :
% ===========
% Thiel, M., Féraud, B., & Govaerts, B. (2017). ASCA+ and APCA+:
% Extensions of ASCA and APCA in the analysis of unbalanced multifactorial
% designs: Analyzing unbalanced multifactorial designs with ASCA+ and APCA+.
% Journal of Chemometrics, 31(6), e2895. https://doi.org/10.1002/cem.2895
%
% Input arguments :
% =================
% X : data matrix with samples in the rows and variables in the columns
%     matrix (n x p) to be decomposed into ANOVA matrices
%
% Y : matrix of factors in the columns and levels for each sample
%     identified by integers (n x q)
%
% Options : field structure containing optional parameters
%   Options.decomp : ANOVA decomposition method 'glm' or 'classical'
%   Options.interactions : the maximum number of interactions to consider
%
% Output arguments :
% ==================
% adecomp_res : structure containing results of the APCA
%   adecomp_res.avgmat : cell structure of factor/interaction effects
%       (see adecomp_res.matid for information on factor/interactions)
%
%   adecomp_res.avgmat{1,1}=Global means
%   adecomp_res.avgmat{1,2}=Means for F1
%
%   adecomp_res.anovamat : ANOVA-deflated matrices contained in cells
%       (the first cell contains the original X)
%   adecomp_res.anovamat{1,7}=adecomp_res.avgmat{1,8}+adecomp_res.anovamat{1,8}
%
%   adecomp_res.anovamat{1,1}=X
%   adecomp_res.anovamat{1,2}=X deflated by Means for F1
%   adecomp_res.anovamat{1,2}=X-adecomp_res.avgmat{1,2}
%
%   adecomp_res.Xs{1,1}=Means for F1+Residuals
%   adecomp_res.Xs{1,1}=adecomp_res.avgmat{1,2}+Residuals
%
%   adecomp_res.matid : identity of ANOVA matrices
%       (grand mean, effects, interactions and residual pure error)
%
%   adecomp_res.ssq : sum of squares for each effect based on adecomp_res.avgmat
%   adecomp_res.ssqvarexp : variance explained by the sum of squares
%
%   adecomp_res.Yint : stores the design matrix for interactions
%   adecomp_res.Grps : Values for for each Factor & Interaction
%   adecomp_res.Options : Options used to perform APCA
%
% Usage :
% =======
% Options.decomp = 'glm'; % or 'classical'
% Options.interactions = 2;
% [adecomp_res] = adecomp_DNR(X, Y, Options)
%
% Related function :
% ==================
% apca.m
% asca.m
% acomdim.m
% anovatp.m
% sumcoding.m (creates design matrix according to GLM methodology)
%
% Author :
% Miguel de Figueiredo
% @ : defiguem@gmail.com
%
% Modifications:
% ==============
% 31/01/21
% include 'R' as identifier for residuals in corresponding position for
% adecomp.matid
% LMS
% =========================================================================

%% Fail-safe section

[n,p] = size(X); % size of X
[~,q] = size(Y); % size of Y

% Checks if a Options structure exists
if exist('Options','var') == 0
    Options = struct;
end

% Checks if the maximum number of interactions to consider was defined
if isfield(Options,'interactions') == 0
    Options.interactions = 2;
elseif isfield(Options,'interactions') == 1 && Options.interactions > q
    Options.interactions = q; % maximum number of interactions cannot be > q
end


%% Multivariate ANOVA decomposition

%%%%% Put X in first cell
adecomp_res.anovamat{1,1}=X;

% Switches cases depending on the multivariate ANOVA decomposition method
switch Options.decomp
    
    case 'classical'
        
        % Grand mean matrix : initialization of decomposition into ANOVA matrices
        adecomp_res.avgmat{1,1} = repmat(mean(X,1), n,1);
        adecomp_res.matid{1,1} = 'Mean'; % grand mean matrix as index 0
        
        adecomp_res.anovamat{1,2} = X - adecomp_res.avgmat{1,1};
        
        % anovamat for variance of centred matrix
        adecomp_res.ssq{1,1} = sum(adecomp_res.anovamat{1,2}(:).^2); % total sum of squares
        adecomp_res.ssqvarexp{1,1} = (adecomp_res.ssq{1,1} / adecomp_res.ssq{1,1}) * 100; % explained variance
        
        % Calculation of the main factor effects : iterates over columns of Y
        for i = 1 : q
            adecomp_res.Yint{1,i} = Y(:,i); % stores factor/level classes
            Y_inter{i}=Y(:,i);
            
            niveaux=unique(Y(:,i)); % levels for factor Y(:,i)
            
            % Creates matrix of averages for each level j of factor i
            adecomp_res.anovamat{1,i+2} = zeros(n,p);
            adecomp_res.avgmat{1,i+1} = zeros(n,p);
            
            for j = 1 : length(niveaux)
                idx = find(Y(:,i) == niveaux(j)); % indices of samples at factor i and level j
                levmean = mean(adecomp_res.anovamat{1,i+1}(idx,:),1); % mean of factor i at level j
                
                adecomp_res.avgmat{1,i+1}(idx,:) = repmat(levmean, length(idx),1);
            end
            
            % Replaces inplace the matrix of averages by decomposition from previous one
            adecomp_res.anovamat{1,i+2} = adecomp_res.anovamat{1,i+1} - adecomp_res.avgmat{1,i+1};
            adecomp_res.matid{1,i+1} = num2str(i); % factor index associated to Y column
            
            % avgmat for variance extracted from centred matrix
            adecomp_res.ssq{1,i+1} = sum(adecomp_res.avgmat{1,i+1}(:).^2); % sum of squares
            
            %             temp1 = sum(adecomp_res.anovamat{1,i+1}(:).^2); % sum of squares
            %             temp2 = sum(adecomp_res.anovamat{1,i+2}(:).^2); % sum of squares
            %             adecomp_res.ssq{1,i+1} = temp1-temp2; % sum of squares
            
            adecomp_res.ssqvarexp{1,i+1} = (adecomp_res.ssq{1,i+1} / adecomp_res.ssq{1,1}) * 100; % explained variance
            
        end
        
        % Calculates interactions only if Options.interactions >= 2
        if Options.interactions >= 2
            
            % Defines interactions matrix according to Options.interactions
            % Brilliant !
            combs = {};
            for i = 2 : Options.interactions
                combs = [combs; num2cell(nchoosek(1:q, i),2)];
            end
            
            tab = length(adecomp_res.avgmat)+1; % Need to add one position to last table
            
            for i = 1 : length(combs)
                Y_inter{i+q}=Y(:,combs{i});
                Ycomb = cellstr(num2str(Y_inter{i+q})); % interaction of factors
                
                adecomp_res.Yint{1,i+q} = Ycomb; % stores factor/level classes
                
                uni = unique(Ycomb); % unique combinations of levels in factors
                
                % Creates matrix of averages for each level j of factor i
                adecomp_res.anovamat{1,tab+1} = zeros(n,p);
                adecomp_res.avgmat{1,tab} = zeros(n,p);
                
                for j = 1 : length(uni)
                    idx = find(strcmp(Ycomb, uni{j}) == 1); % indices of samples at factor i and level j
                    
                    %%%%% MDF
                    % levmean = mean(adecomp_res.anovamat{1,tab-1}(idx,:),1); % mean of factor i at level j
                    levmean = mean(adecomp_res.anovamat{1,tab}(idx,:),1); % mean of factor i at level j
                    
                    adecomp_res.avgmat{1,tab}(idx,:) = repmat(levmean, length(idx),1);
                end
                
                % Replaces in place the matrix of averages by decomposition from previous one
                adecomp_res.anovamat{1,tab+1} = adecomp_res.anovamat{1,tab} - adecomp_res.avgmat{1,tab};
                adecomp_res.matid{1,tab} = num2str(combs{i}); % factor index associated to Y column
                
                % avgmat for variance extracted from centred matrix
                adecomp_res.ssq{1,tab} = sum(adecomp_res.avgmat{1,tab}(:).^2); % sum of squares
                
                %                 temp1 = sum(adecomp_res.anovamat{1,tab}(:).^2); % sum of squares
                %                 temp2 = sum(adecomp_res.anovamat{1,tab+1}(:).^2); % sum of squares
                %                 adecomp_res.ssq{1,tab} = temp1-temp2; % sum of squares
                
                adecomp_res.ssqvarexp{1,tab} = (adecomp_res.ssq{1,tab} / adecomp_res.ssq{1,1}) * 100; % explained variance
                
                tab = tab + 1; % updates table position
            end
            
        end
        
        
    case 'glm'
        
        % Sum coding of the input design matrix according to Thiel et al. (2017)
        [sumcod, sumcodempty] = sumcoding(Y, Options);
        adecomp_res.sumcod = sumcod;
        
        % Creates empty B (it is juste useful for later calculations...)
        bempty = {};
        for i = 1 : length(sumcod)
            bempty{i} = zeros(p,size(sumcod{i},2));
        end
        
        % Grand mean matrix : initialization of decomposition into GLM matrices
        B = bempty;
        B{1} = (pinv(sumcod{1}' * sumcod{1}) * sumcod{1}' * X)';
        Xf = sumcodempty;
        Xf{1} = sumcod{1};
        
        adecomp_res.avgmat{1,1} = cat(2,Xf{:}) * cat(2,B{:})';
        
        adecomp_res.anovamat{1,2} = X - adecomp_res.avgmat{1,1};
        adecomp_res.matid{1,1} = 'Mean'; % grand mean matrix as index 0
        
        % anovamat for variance of centred matrix
        adecomp_res.ssq{1,1} = sum(adecomp_res.anovamat{1,2}(:).^2); % total sum of squares
        adecomp_res.ssqvarexp{1,1} = (adecomp_res.ssq{1,1} / adecomp_res.ssq{1,1}) * 100; % explained variance
        
        % Calculation of the main factor effects : iterates over columns of Y
        for i = 1 : q
            adecomp_res.Yint{1,i} = Y(:,i); % stores factor/level classes
            Y_inter{i}=Y(:,i);
            
            B = bempty;
            B{i+1} = (pinv(sumcod{i+1}' * sumcod{i+1}) * sumcod{i+1}' * adecomp_res.anovamat{1,i+1})';
            
            Xf = sumcodempty;
            Xf{i+1} = sumcod{i+1};
            
            adecomp_res.avgmat{1,i+1} = cat(2,Xf{:}) * cat(2,B{:})';
            
            % Replaces inplace the matrix of averages by decomposition from previous one
            adecomp_res.anovamat{1,i+2} = adecomp_res.anovamat{1,i+1} - adecomp_res.avgmat{1,i+1};
            adecomp_res.matid{1,i+1} = num2str(i); % factor index associated to Y column
            
            % avgmat for variance extracted from centred matrix
            adecomp_res.ssq{1,i+1} = sum(adecomp_res.avgmat{1,i+1}(:).^2); % sum of squares
            
            %             temp1 = sum(adecomp_res.anovamat{1,i+1}(:).^2); % sum of squares
            %             temp2 = sum(adecomp_res.anovamat{1,i+2}(:).^2); % sum of squares
            %             adecomp_res.ssq{1,i+1} = temp1-temp2; % sum of squares
            
            adecomp_res.ssqvarexp{1,i+1} = (adecomp_res.ssq{1,i+1} / adecomp_res.ssq{1,1}) * 100; % explained variance
            
        end
        
        % Calculates interactions only if Options.interactions >= 2
        if Options.interactions >= 2
            
            % Defines interactions matrix according to Options.interactions
            combs = {};
            for i = 2 : Options.interactions
                combs = [combs; num2cell(nchoosek(1:q, i),2)];
            end
            
            tab = length(adecomp_res.avgmat)+1; % Need to add one position to last table
            
            for i = 1 : length(combs)
                
                Y_inter{i+q}=Y(:,combs{i});
                Ycomb = cellstr(num2str(Y_inter{i+q})); % interaction of factors
                
                adecomp_res.Yint{1,i+q} = Ycomb; % stores factor/level classes
                
                B = bempty;
                
                %%%%% MDF
                % B{tab} = (pinv(sumcod{tab}' * sumcod{tab}) * sumcod{tab}' * adecomp_res.anovamat{1,tab-1})';
                B{tab} = (pinv(sumcod{tab}' * sumcod{tab}) * sumcod{tab}' * adecomp_res.anovamat{1,tab})';
                
                Xf = sumcodempty;
                Xf{tab} = sumcod{tab};
                
                adecomp_res.avgmat{1,tab} = cat(2,Xf{:}) * cat(2,B{:})';
                
                % Replaces in place the matrix of averages by decomposition from previous one
                adecomp_res.anovamat{1,tab+1} = adecomp_res.anovamat{1,tab} - adecomp_res.avgmat{1,tab};
                adecomp_res.matid{1,tab} = num2str(combs{i}); % factor index associated to Y column
                
                % avgmat for variance extracted from centred matrix
                adecomp_res.ssq{1,tab} = sum(adecomp_res.avgmat{1,tab}(:).^2); % sum of squares
                
                %                 temp1 = sum(adecomp_res.anovamat{1,tab}(:).^2); % sum of squares
                %                 temp2 = sum(adecomp_res.anovamat{1,tab+1}(:).^2); % sum of squares
                %                 adecomp_res.ssq{1,tab} = temp1-temp2; % sum of squares
                
                adecomp_res.ssqvarexp{1,tab} = (adecomp_res.ssq{1,tab} / adecomp_res.ssq{1,1}) * 100; % explained variance
                
                tab = tab + 1; % updates table position
                
            end
        end
end
ntab=size(adecomp_res.avgmat,2);


for i = 1 : ntab-1
    Xs{1,i} = adecomp_res.avgmat{1,i+1} + adecomp_res.anovamat{1,end};
    
    Grps_tmp=0;
    for j=1:size(Y_inter{i},2)
        Temp=Grps2y(y2Grps(Y_inter{i}(:,j)));
        Grps_tmp=Grps_tmp+Temp*3^(j-1);
    end
    Grps(:,i)=Grps2y(y2Grps(Grps_tmp));
    
end
adecomp_res.Grps =Grps; % stores factor/level classes

% Adds the residuals matrix as a last table
Xs{1,ntab} = adecomp_res.anovamat{1,end};

% Matrices with Means + Residues
adecomp_res.Xs=Xs;

% place 'R' for residuals into position of adecom_res.matid
[~, ntab]=size(adecomp_res.matid);
adecomp_res.matid{1,ntab+1}='R';


end