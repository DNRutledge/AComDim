function [A_ComDimIn, A_ComDimOut, ComDimResult, A_Collection]=A_ComDim_2020_LMS(X, lbl_var, GrpsCritNum, Options);
% A_ComDim - Performs ComDim analysis on a series of tables
% generated as in ASCA or APCA.
% function [A_ComDimIn, A_ComDimOut]=A_ComDim(X, GrpsCritNum, Options)
%
%
% INPUT :
%---------------
% X:  (nXrows x nXcols) matrix of initial data to be decomposed into ANOVA matrices
% lbl_var: variable names for loadings plots
% GrpsCritNum: (nXrows x nGrps) matrix of group numbers for each classification criterion
% Options.Plots: 0=No figures; 1=Last figures; 2=All figures; 3=Too many figures
% % % % % Options.Dims: Number of Common Components if <> Number of tables
% Options.ndim: Number of Common Components if <> Number of tables
% Options.Samples=0/1; Plot samples in ellipses or not
% Options.decomp: Method used to deflate X matrix : 'classical' (ANOVA-default)/'glm' (DoE)
% Options.MVA = 'PCA', 'ICA', 'CCA' (& other MVA methods ?)
%
% If Options.MVA = 'PLS'
% then Options.Y = data vector or matrix for ComDim_PLS
%
% OUTPUT :
%-----------------
% A_ComDimIn: Structure with fields:
% A_ComDimIn.X : (nXrows x nXcols) original data matrix (X)
% A_ComDimIn.Grps : (nXrows x nGrps) groups numbers for factors and interactions
% A_ComDimIn.X_i : (nTables+1 of nXrows x nXcols) Deflated Matrices before & after removing X_mean
% A_ComDimIn.X_r : (nTables of nXrows x nXcols) Matrices of group mean vectors + residual noise
% A_ComDimIn.Titre : (nTables x 1) list of strings with labels of Factors;
% A_ComDimIn.TotalVariance= : (nTables+1 x 1) % Extracted Variance
%
% A_ComDimOut: Cell array of (1 x nTables) Structures with fields:
% A_ComDimOut{i}.ComDimQ : (nXrows x 1) observations scores
% A_ComDimOut{i}.ComDimSalience : (1 x ndim) weight of the original tables in each dimension.
% A_ComDimOut{i}.ComDim_Lx : (1 x nXcols) weight of the original variables for the significant table.
% A_ComDimOut{i}.ComDimExpl : (1 x 1) percentage explanation given by dimension
% A_ComDimOut{i}.F1 : (1 x 1) Fisher's F based on saliences of CD1
%
% A_ComDimOut{i}.MaxTab : (1 x 1) Table with highest Salience for dimension i.
%
%
% ComDimResult: Structure with all the output of comdim_calib_V7
%
% A_Collection : (nGrps x 1) Structure of (nXrows x nXcols) Factor+Error tables
%
% CALLS :
% adecomp_DNR : to decompose X by ANOVA or GLM
% comdim_PCA_2020
% comdim_ICA_2020
% ANOVAPCAmegX_Ellipse : 2D ellipses of scores
% ANOVAPCAmegX_Ellipse_no_text_sample_mod
% Figure_DNR
%
% REFERENCES
%-----------------
% 'ComDim' method :
% E.M. Qannari, I. Wakeling, P. Courcoux and H. J. H. MacFie
% Food quality and Preference 11 (2000) 151-154
%
% 'APCA' methods :
% P. Harrington, N. Vieira, J. Espinoza, J. Nien, R. Romero, A. Yergey
% Analysis of variance-principal component analysis: A soft tool for proteomic discovery.
% Analytica Chimica Acta, 544, 118-127, 2005.
%
%
% TYPICAL EXAMPLE
%-----------------
% A (37x4600) matrix of spectra 'Wine_Data'
% and
% 4 grouping criteria indicated by a (37x4) matrix
% of non-zero groups numbers 'Wine_Grps'
% and
% 0 plots, 10 Common Components
% Opts.Plots=0;
% Opts.Dims=10;
% Opt.Samples='Samples';

%[A_ComDimIn, A_ComDimOut, A_Collection]=A_ComDim_2020(Wine_Data, Wine_Grps, Opts);
%
% LMS 27 Jan 21 updated F plot to include critical values
%
%%

[nXrows,nXcols]=size(X);
[nXrows, nGrpCrits]=size(GrpsCritNum);
% verify matrix sizes later !

GrpsCritNumTemp=GrpsCritNum;

if exist('Options')
    if isfield(Options, 'Plots')
        Plots=Options.Plots;
    else
        Plots=3;
    end
    
    %     if isfield(Options, 'Dims')
    %         nDims=Options.Dims;
    %     else
    %         nDims=2^(nGrpCrits);
    %     end
    
    if isfield(Options, 'ndim')
        ndim=Options.ndim;
    else
        ndim=2^(nGrpCrits);
    end
    
    if isfield(Options, 'Samples')
        PlotSamples=Options.Samples;
    else
        PlotSamples=0;
    end
    
    if isfield(Options, 'decomp')
        decomp=Options.decomp;
    else
        decomp='classic';
    end
    
    if ~isfield(Options,'MVA')
        Options.MVA='PCA';
        MVA=Options.MVA;
    else
        MVA=Options.MVA;
    end
    if MVA=='PLS'
        if isfield(Options, 'Y')
            Y=Options.Y;
        else
            Y=[];
            error('Y required for PLS');
        end
    end
    
else
    Plots=3;
    ndim=2^(nGrpCrits);
end


%% "ANOVA" or "GLM" decomposition of X
clear adecomp_res;

% Options.decomp='classical';
% Options.interactions=3;

[adecomp_res] = adecomp_DNR(X, GrpsCritNum, Options); % see adecomp.m for information

%% Number of tables
NumTab=size(adecomp_res.avgmat,2);
[nXrows,nXcols]=size(adecomp_res.avgmat{1});

%% Plot matrices before and after decomposing
if Plots==3
    for col=1:NumTab
        Figure_DNR(1);
        
        subplot(3,1,1);
        plot(adecomp_res.anovamat{1, col}'), axis tight;
        title(['Before ', num2str(adecomp_res.matid{col}),' removed']);
        
        subplot(3,1,2);
        plot(adecomp_res.avgmat{1, col}'), axis tight;
        title(adecomp_res.matid(col));
        title(num2str(adecomp_res.matid{col}));
        
        subplot(3,1,3);
        plot(adecomp_res.anovamat{1, col+1}'), axis tight;
        title(['After ', num2str(adecomp_res.matid{col}),' removed']);
    end
end


%%
clear A_ComDimIn;

A_ComDimIn.X=adecomp_res.anovamat{1, 1}; % Origiunal X
A_ComDimIn.Grps=adecomp_res.Grps; % Group values
A_ComDimIn.X_i=adecomp_res.anovamat; % X & Deflated matrices
A_ComDimIn.X_r=adecomp_res.Xs; % Factor means + Residues
A_ComDimIn.Titre=adecomp_res.matid; % Factor names
A_ComDimIn.TotalVariance=cell2mat(adecomp_res.ssqvarexp); % Extrcated Variance

Grps=adecomp_res.Grps;


%% Prepare dataset for ComDim
clear collection;

for i=1:NumTab
    collection(i).v=[1:nXcols];
    collection(1).i=[1:nXrows];
    collection(i).d=A_ComDimIn.X_r{1,i};
end

%% ComDim
clear ComDimResult A_Collection;

Options.ndim=ndim;
Options.Output='TPL';

switch MVA
    case 'PCA'
        %% Do PCA
        [ComDimResult]=comdim_PCA_2020(collection, Options);
        
    case 'CCA'
        %% Do CCA
        [ComDimResult]=comdim_CCA_2020(collection, Options);
        
    case 'ICA'
        %% Do ICA
        [ComDimResult]=comdim_ICA_2020(collection, Options);
        
    case 'PLS'
        %% Do PLS
        [ComDimResult]=comdim_PLS_2020(collection, Y, Options);
end

A_Collection=collection;

%%
for i=1:ndim
    A_ComDimOut{i}.ComDimQ=ComDimResult.Q.d(:,i);
    A_ComDimOut{i}.ComDimSalience=ComDimResult.saliences.d(:,i);
    
    %     [maxSalience,IndSalience]=max(ComDimResult.saliences.d(:,i));
    [val,MaxTab(i)]=max(ComDimResult.saliences.d(:,i));
    U=inv(A_ComDimOut{1,i}.ComDimQ(:,1)'*A_ComDimOut{1,i}.ComDimQ(:,1))*A_ComDimOut{1,i}.ComDimQ(:,1)';
    
    %     A_ComDimOut{i}.ComDim_Lx=U* A_ComDimIn.X_i{1,IndSalience+1};
    A_ComDimOut{i}.ComDim_Lx=U* A_ComDimIn.X_i{1,MaxTab(i)+1};
    
    A_ComDimOut{i}.ComDimExpl=ComDimResult.explained.d(:,i);
    
    A_ComDimOut{i}.MaxTab=MaxTab(i);
    
end


%%
% for i=1:ndim
%     [val,MaxTab(i)]=max(ComDimResult.saliences.d(:,i));
% end

%% Plot ComDim results
nCols=2;
if Plots>0
    Figs=double(floor((ndim+1)/nCols));
    
    Figure_DNR(1);
    for i=1:ndim
        subplot(Figs,nCols,i);
        bar(ComDimResult.saliences.d(:,i),'r','BarWidth',0.2), axis tight;
        set(gca,'XTickLabel',A_ComDimIn.Titre(2:end));
        set(gca,'XTickLabelRotation',90);
        temp=xlim;
        temp(1,1)=0;
        xlim(temp)
        ylabel(['CC', num2str(i)]);
        if MaxTab(i)<NumTab
%             title(['Salience CC', num2str(i) ' Factor ' ,A_ComDimIn.Titre{MaxTab(i)+1}]);
            title(['Factor ' ,A_ComDimIn.Titre{MaxTab(i)+1}]);
        end
    end
    suptitle('Saliences');
    
    Figure_DNR(1);
    for i=1:ndim
        subplot(Figs,nCols,i);
        plot(ComDimResult.Q.d(:,i),'b-'), axis tight;
        ylabel(['CC', num2str(i)]);
        if MaxTab(i)<NumTab
            hold on;
            scatter([1:nXrows],ComDimResult.Q.d(:,i),20,Grps(:,MaxTab(i)),'Filled');
            % text([1:nXrows],ComDimResult.Q.d(:,i),num2str(Grps(:,MaxTab(i))));
            
            %         end
            %         if MaxTab(i)<NumTab
            
            title(['Factor ' ,A_ComDimIn.Titre{MaxTab(i)+1}]);
        end
    end
    suptitle('Scores');
    
    Figure_DNR(1);
    for i=1:ndim
        subplot(Figs,nCols,i);
        plot(A_ComDimOut{i}.ComDim_Lx,'MarkerSize',2), axis tight;
        ylabel(['CC', num2str(i)]);
        if MaxTab(i)<NumTab
            title(['Factor ' ,A_ComDimIn.Titre{MaxTab(i)+1}]);
        end
    end
    suptitle('Loadings');
    
    Figure_DNR(1);
    subplot(1,3,1);
    plot(ComDimResult.explained.d,':s'), axis tight;
    xlabel('Dims');
    ylabel('Variance');
    subplot(1,3,2);
    plot(ComDimResult.Sum_saliences_Dim.d,':s'), axis tight;
    xlabel('Dims');
    ylabel('Sum Saliences');
    subplot(1,3,3);
    bar(ComDimResult.Sum_saliences_Tab.d'), axis tight;
    xlabel('Tables');
    ylabel('Sum Saliences');
    
    %% Calculate F from Saliences
    for i=ndim:-1:1
        Fup=ones(NumTab,1)*ComDimResult.saliences.d(NumTab,i);
        Fdown=ComDimResult.saliences.d(:,i);
        
        F1=Fup./Fdown;
        
        A_ComDimOut{i}.F1=F1;
    end
    
    %% Sum saliences of all CCs that have max salience for Residues
    [C,I]=max(ComDimResult.saliences.d);
%     SumSalResidues=sum(C(:,I==NumTab));
%     SumSalAll=sum(ComDimResult.saliences.d(:,I==NumTab)');
    MaxSalRes=find(I==NumTab);
    SumSalResidues=sum(C(:,MaxSalRes));
    
    if length(MaxSalRes)>1
        Fdown_All=sum(ComDimResult.saliences.d(:,MaxSalRes)')';
    else
        Fdown_All=(ComDimResult.saliences.d(:,MaxSalRes));
    end
    
    for i=ndim:-1:1
        Fup=ones(NumTab,1)*SumSalResidues;
        Fdown=ComDimResult.saliences.d(:,i);
        
        F1_All=Fup./Fdown;
        A_ComDimOut{i}.F1_All=F1_All;
        
        F_All=Fup'./Fdown_All';
        A_ComDimOut{i}.F_All=F_All;
    end
    
    %% Plot F1 & F1_All
    Figure_DNR(1);
    subplot(2,1,1);
    critF95=finv(0.95, nXrows-1, nXrows-1);
    critF99=finv(0.99, nXrows-1, nXrows-1);
    
    % axes('YScale','log','YMinorTick','on');
    % axes('YMinorTick','on');
    
    bar(F1(:,1),'r');
    set(gca,'XTickLabel',A_ComDimIn.Titre(2:end));
    set(gca,'XTickLabelRotation',90);
    hold on;
    line([0 size(F1,1)+1],[1 1],'LineWidth',2,'LineStyle',':', 'Color', 'black');
    line([0 size(F1,1)+1], [critF95 critF95], 'LineWidth',2,'LineStyle','--', 'Color', 'green');
    line([0 size(F1,1)+1], [critF99 critF99], 'LineWidth',2,'LineStyle','--', 'Color', 'blue');
    
    title('F* based on CC1 Salience ratios');
    ylabel('F*');
    xlabel('Factor');
    legend('Factor F value', 'F=1', 'Critical F p<0.05', 'Critical F p<0.01')
    
    subplot(2,1,2);
    critF95=finv(0.95, nXrows-1, nXrows-1);
    critF99=finv(0.99, nXrows-1, nXrows-1);
    
    % axes('YScale','log','YMinorTick','on');
    % axes('YMinorTick','on');
    
%     bar(F1_All(:,1),'r');
    bar(F_All,'r');
    set(gca,'XTickLabel',A_ComDimIn.Titre(2:end));
    set(gca,'XTickLabelRotation',90);
    hold on;
    line([0 size(F1_All,1)+1],[1 1],'LineWidth',2,'LineStyle',':', 'Color', 'black');
    line([0 size(F1_All,1)+1], [critF95 critF95], 'LineWidth',2,'LineStyle','--', 'Color', 'green');
    line([0 size(F1_All,1)+1], [critF99 critF99], 'LineWidth',2,'LineStyle','--', 'Color', 'blue');
    
    title(['F* based on Salience ratios for CCs ',num2str(MaxSalRes)]);
    ylabel('F*');
    xlabel('Factor');
    legend('Factor F value', 'F=1', 'Critical F p<0.05', 'Critical F p<0.01')
    
end


%% Plot ComDim
if Plots>1
    for i=1:ndim
        if MaxTab(i)<NumTab
            Var=floor(A_ComDimIn.TotalVariance(1,MaxTab(i)+1)*100)/100;
            
            Figure_DNR(1);
            
            if PlotSamples==0
                ANOVAPCAmegX_Ellipse_no_text_sample_mod([ComDimResult.Q.d(:,i),ComDimResult.Q.d(:,1)],Grps(:,MaxTab(i)),0.05,1,2); %'_mod' change LMS as this is faster
                %                 ANOVAPCAmegX_Ellipse_no_text_sample
            elseif PlotSamples==1
                ANOVAPCAmegX_Ellipse([ComDimResult.Q.d(:,i),ComDimResult.Q.d(:,1)],Grps(:,MaxTab(i)),0.05,1,2);
                %                 ANOVAPCAmegX_Ellipse
            end
            
            xlabel(['CC' num2str(i),' / SumSal= ',num2str(ComDimResult.Sum_saliences_Dim.d(1,i))]);
            ylabel(['CC' num2str(1),' / SumSal= ',num2str(ComDimResult.Sum_saliences_Dim.d(1,1))]);
            %             ylabel('CC1');
            %             title(['Factor ' A_ComDimIn.Titre{MaxTab(i)+1}]);
            title(['Factor ',num2str(A_ComDimIn.Titre{MaxTab(i)+1}), ' (Total Variance ',num2str(Var),' %)']);
        end
    end
end

%% Plot loadings as bar style
% need to add catches to stop errors if labels are incorrect length or if
% loadings are spectra
if Plots>2
    for i=1:ndim
        if MaxTab(i)<NumTab
            Var=floor(A_ComDimIn.TotalVariance(1,MaxTab(i)+1)*100)/100;
            
            Figure_DNR(1);
            
            h=bar(A_ComDimOut{i}.ComDim_Lx);
            [temp, nv]=size(A_ComDimOut{i}.ComDim_Lx);
            
            set(gca, 'XtickLabels', lbl_var);
            set(gca, 'YtickLabels', '');
            title(['Loadings CC', num2str(i)  ' Factor ' ,A_ComDimIn.Titre{MaxTab(i)+1}]);
            
            xlabel(['CC' num2str(i),' / SumSal= ',num2str(ComDimResult.Sum_saliences_Dim.d(1,i))]);
            %ylabel(['CC' num2str(1),' / SumSal= ',num2str(ComDimResult.Sum_saliences_Dim.d(1,1))]);
        end
    end
    
    
end


