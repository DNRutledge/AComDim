function ANOVAPCAmegX_Ellipse_no_text_sample_mod(score,IndClass,confi,CP1,CP2);
%==========================================================================
% Input  :
%---------
% score     : Les composantes principales soit de l'AFD, l'ACP, PLS, etc
% IndClasse : Variable qualitative qui indique l'indice de classe de chaque
%             observation
% confi     : le seuil de confiance de l'ellipse
%             exemple :pour avoir 95% de l'ellipse il suffit que confi=0.05
% CP1 CP2   : Les 2 Composantes à choisir pour la représentation graphiques
%
% on peut changer les lignes 156, 163 et 166 (pour dessiner des lignes ou des points.
% les points sont plus rapides)
% we can change the lines 156, 163 and 166 ( to draw lines or points. Points are faster )
%
%=========================== Plan Factoriel à 2 Dim =======================
% Nocairi Hicham
% 27-09-2006
%==========================================================================

%subplot(1,1,1,'Position',[0.08 0.08 0.87 0.87])
% ncize  = 6.5; % la taille des identifiants du nuages des individus
ncize  = 10; % la taille des identifiants du nuages des individus; the size of the identifies in the group clouds
ncize1 = 10; % la taille des labels des axe x et y et le titre; the size of the text for x/y labels title
csize  = 7.5;
ecrit  = 'Times';
ecrit  = 'Hevetica';
ecrit  = 'Comic Sans MS Gras';
Inter  = 'Tex';
if nargin==0
    error('il n''y a pas assez d''arguments en entrée');
end
% initialisation 
ellip  = []; iclass  = []; ikamsup = []; ikamQltsup = []; mx      = [];
bary   = []; cercle  = []; segm    = []; qlt        = []; droit   = [];
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% DEMARCHE :
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%- les données d'analyse :
% X       = pcaRes.Xdata;
% [n, p]  = size(pcaRes.score); % les composantes principales

%-les composantes principales  (CP1, CP2) :
Z       = score(:,[CP1, CP2]); % get scores for first 2 LV
QLT     = IndClass; %vector of class numbers

%  (Création de la couleur blue pour l'ensemble des observations :
CoulObs = ones(size(Z,1),1)*[0 0 1];
%
%--------------------------------------------------------------------------
% Création d'une couleur pour chaque classe :
Class     = IndClass;
m         = size(Class,1);
numColors = 1;
theGroups = 1;
groups    = 0;
cmap      = [0 0 1];
groups    = m;%max(Class.d);
if groups >= 1 & groups <= m
    theGroups = zeros(m,1);
    numColors = 0;
    for count = groups:-1:1
        if (theGroups(count) == 0)
            P = zeros(m,1);
            P(count) = 1;
            uu = find(Class==Class(count,1));
            P(uu,1) = ones(size(uu,1),1);            
            numColors = numColors + 1;
            theGroups(logical(P)) = numColors;
        end
    end
    cmap = hsv(numColors);
    cmap(end+1,:) = [0 0 0];
else
    groups = 1;
end

CoulIndCl = zeros(m,3);
h            = zeros(m,1);

for n = 1:(m)
    if n <= groups
        CoulIndCl(n,:) = cmap(theGroups(n),:);
    else
        CoulIndCl(n,:) = cmap(end,:);
    end
end
for n=1:m
    if CoulIndCl(n,:)==[1 1 0]
        CoulIndCl(n,:) = [0.53 0.32 0.20];
    end
end


% create a matrix (Coulg) of n classes x 3 for colours of classes

% use unique function

Coulg = unique(CoulIndCl, 'rows');
%================================(FIGURE 1)================================
% PLAN Factoriel A DEUX DIMENSIONS (CP1,CP2) :
%--------------------------------------------------------------------------
%plot(Z(:,1), Z(:,2), ' w') ;
%hold on

margin = 0.05;
minx   = min(Z(:,1));               miny   = min(Z(:,2));    
maxx   = max(Z(:,1));               maxy   = max(Z(:,2));
minx   = minx-(maxx-minx)*margin;   miny   = miny-(maxy-miny)*margin;
maxx   = maxx+(maxx-minx)*margin;   maxy   = maxy+(maxy-miny)*margin;

%axis([minx maxx miny maxy]);
%line([minx maxx],[0 0]); 
%line([0 0],[miny maxy]);

% Projection du nuage des individus sur le plan factoriel (CP1,CP2) :
%for i=1:size(Z,1)
%    ikam(i) = text(Z(i,1),Z(i,2),'\bullet','FontSize',ncize,'Fontweight','Demi',...
%        'FontName',ecrit,'Interprete',Inter,'Color',CoulObs(i,:));
  
%end
%set(gca,'Fontsize',ncize1,'Fontweight','Demi')
%labx   = [char('t') num2str(CP1)];
%laby   = [char('t') num2str(CP2)];
%xlabel(labx);
%ylabel(laby);

% title('Plan factoriel à 2 Dimension' )

%axis square;
%
%--------------------------------------------------------------------------

%----------------------------------------------------------
% AFFICHER L'ELLIPSE DE CONFIANCE DE CHAQUE GROUPE
%__________________________________________________________
%if ~isempty(ikam)       delete(ikam);       ikam       = []; end
[n, nc]        = size(IndClass);
% get number of classes in IndClass
no_class=size(getlabels(nominal(IndClass)));

[n1,p1]  = size(Z);
k0       = 0;
k3       = 0;



%for k=1:n %the number of samples ??? but should this be the number of classes
for k=1:no_class(1,2);

    
    l=1;
    
   % l    = find(IndClass==IndClass(k,1)); % gets sample placeholders for all samples in same class 

    if ~isempty(l)
       
        
        X         = Z((IndClass==k),:); % get scores for class objects (samples)
         nbrInd    = size(X,1); % number in class
         
         %Coulg     = CoulIndCl(l,:); % colour code for class
         % change Coulg to number fo classes not per sample
        nn1       = size(X,1); % number of class samples/objects
        %Cl2       = IndClass(l,1); 
        %[tf ind1] = ismember(k,IndClass,'rows');

        k0        = k0+1;
        G          = mean(X); %mean score for class 
        %Cl         = [char('g') num2str(k)]; % not sure what this does not used
        % une ACP du tableau centré X de la classe k
        XCR        = (X - ones(nn1,1)*G) ;  % mean centre scores      
        Q          = (XCR'*XCR)/(nn1-1); 
        Tr         = trace(Q);
        [V, S, H]  = svd(Q);
        [u1, idx]  = sort(diag(S)') ;
        u1         = fliplr(u1) ;
        idx        = fliplr(idx) ;
        eigenval   = u1(1:2);%%%%%%%%%%%%%%%%%%%%%%%%
        idx        = idx(1:2);
        vec        = V(:,idx); 
        l1         = eigenval(1);
        l2         = eigenval(2);
        %t          = 0:2*pi/128:2*pi;%CHANGE number in the middle TO 0.001
        t          = 0:0.001:2*pi;%CHANGE number in the middle TO 2*pi/128
        c          = tinv(1-(confi/2),nn1-1); % l'ellipse à confi% de confiance
        x          = c*sqrt(l1)*cos(t); % create x-vector
        y          = c*sqrt(l2)*sin(t); % create y-vector
        u          = (vec(1)*x-vec(2)*y)+G(1,1); % move vector to centre x-axis
        v          = (vec(2)*x+vec(1)*y)+G(1,2); % move vector to centre y-axis
       
       
        H=patch(u,v,1); %Make a patch object
%         set(H,'LineWidth', 2,'FaceColor',Coulg(k,:), 'EdgeColor' , Coulg(k,:), 'FaceAlpha', 0.1); %Set the color to grey
        set(H,'LineWidth', 2,'FaceColor',Coulg(k,:), 'EdgeColor' , Coulg(k,:), 'FaceAlpha', 0.05); %Set the fill color to lighter than outline
%         set(gca,'Fontsize',ncize1,'Fontweight','Demi');
        set(gca,'Fontsize',ncize1,'Fontweight','Normal');
        hold on;
       
        ikam(k,1)  = plot(X(:,1),X(:,2),'.','MarkerSize', 8, 'Color',Coulg(k,:)); % add labels for samples
        % could move above line out of immediate loop to facilitate correct
        % legend.  make second loop after inserting legend ie place legend
        % on after doing patch/loop and numbers
   
        % temp_txt={'S';  'NS'};
        
        % bary(k0)   = text(G(1,1),G(1,2),temp_txt(k),'FontSize',12, 'Fontweight', 'Bold','FontName',ecrit,'Color',[0 0 0]); %plots centre point
%         bary(k0)   = text(G(1,1),G(1,2),num2str(k),'FontSize',12,'Fontweight','Normal','FontName',ecrit,'Interprete',Inter,'Color',[0 0 0]); %plots centre point
         bary(k0)   = text(G(1,1),G(1,2),num2str(k),'FontSize',12,'Fontweight','Bold','FontName',ecrit,'Interprete',Inter,'Color',[0 0 0]); %plots centre point
        %ikam(l,1)  = text(X(:,1),X(:,2),num2str(k),'FontSize',ncize,'Fontweight','Demi','FontName',ecrit,'Interprete',Inter,'Color',Coulg(1,:)); % add labels for samples
        
        mx(k)      = min(u);Mx(k) = max(u);my(k) = min(v);My(k) = max(v);
        mx(k)      = mx(k)-(Mx(k)-mx(k))*margin;
        my(k)      = my(k)-(My(k)-my(k))*margin;
        Mx(k)      = Mx(k)+(Mx(k)-mx(k))*margin;
        My(k)      = My(k)+(My(k)-my(k))*margin;
       
    end
end;
mx = min([mx,minx]); Mx = max([Mx,maxx]);
my = min([my,miny]); My = max([My,maxy]);
axis([mx Mx my My]);
line([mx Mx],[0 0]);
line([0 0],[my My]);
grid on
%set(H, 'Legend', 'Addd1', 'Add2', 'Add3');
% add individual points
%for k=1:no_class(1,2);

%ikam(k,1)  = plot(X(:,1),X(:,2),'.','MarkerSize', 5, 'Color',Coulg(k,:)); % add labels for samples
%end

