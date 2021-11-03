function ANOVAPCAmegX_Ellipse(score,IndClass,confi,CP1,CP2);

%==========================================================================
% Input  :
%---------
% score     : Les composantes principales soit de l'AFD, l'ACP, PLS, etc
% IndClasse : Variabe qualitative qui indique l'indice de classe de chaque
%             observation
% confi     : le seuil de confiance de l'ellipse
%             exemple :pour avoir 95% de l'ellipse il suffit que confi=0.05
% CP1 CP2   : Les 2 Composantes à choisir pour la représentation graphiques

%=========================== Plan Factoriel à 2 Dim =======================
% Nocairi Hicham
% 27-09-2006
%==========================================================================
%figure
% subplot(1,1,1,'Position',[0.08 0.08 0.87 0.87]);

ncize  = 6.5; % la taille des identifiants du nuages des individus
ncize1 = 10; % la taille des labels des axe x et y et le titre
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
Z       = score(:,[CP1, CP2]);
QLT     = IndClass;

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
%================================(FIGURE 1)================================
% PLAN Factoriel A DEUX DIMENSIONS (CP1,CP2) :
%--------------------------------------------------------------------------
plot(Z(:,1), Z(:,2), ' w') ;
hold on

margin = 0.05;
minx   = min(Z(:,1));               miny   = min(Z(:,2));
maxx   = max(Z(:,1));               maxy   = max(Z(:,2));
minx   = minx-(maxx-minx)*margin;   miny   = miny-(maxy-miny)*margin;
maxx   = maxx+(maxx-minx)*margin;   maxy   = maxy+(maxy-miny)*margin;

axis([minx maxx miny maxy]);
line([minx maxx],[0 0]);
line([0 0],[miny maxy]);

% Projection du nuage des individus sur le plan factoriel (CP1,CP2) :
for i=1:size(Z,1)
    ikam(i) = text(Z(i,1),Z(i,2),'\bullet','FontSize',ncize,'Fontweight','Demi',...
        'FontName',ecrit,'Interprete',Inter,'Color',CoulObs(i,:));

end
set(gca,'Fontsize',ncize1,'Fontweight','Demi')
labx   = [char('t') num2str(CP1)];
laby   = [char('t') num2str(CP2)];
xlabel(labx);
ylabel(laby);
title('Plan factoriel à 2 Dimensions' )
%axis square;
%
%--------------------------------------------------------------------------
grid on
hold on
%----------------------------------------------------------
% AFFICHER L'ELLIPSE DE CONFIANCE DE CHAQUE GROUPE
%__________________________________________________________
if ~isempty(ikam)       delete(ikam);       ikam       = []; end
n        = size(IndClass,1);
[n1,p1]  = size(Z);
k0       = 0;
k3       = 0;

% for k=1:n
%     l    = find(IndClass==IndClass(k,1));

IndClass_Mat=y2Grps(IndClass);
n_Class_Mat=size(IndClass_Mat,2);
for n_Class=1:n_Class_Mat
    [l,temp]=find(IndClass_Mat(:,n_Class)==1);
    if ~isempty(l)
        nbrInd    = size(l,1);
        X         = Z(l,:);
        Coulg     = CoulIndCl(l,:);
        nn1       = size(X,1);
        Cl2       = IndClass(l,1);
        
        %         [tf ind1] = ismember(k,IndClass,'rows');
        
        k0        = k0+1;
        G          = mean(X);
        
        %         Cl         = [char('g') num2str(k)];
        Cl         = [char('g') num2str(n_Class)];
        
        %         % une ACP du tableau centré X de la classe k
        % une ACP du tableau centré X de la classe n_Class
        XCR        = (X - ones(nn1,1)*G) ;
        Q          = (XCR'*XCR)/(nn1-1);
        if ~isinf(Q)
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
            x          = c*sqrt(l1)*cos(t);
            y          = c*sqrt(l2)*sin(t);
            u          = (vec(1)*x-vec(2)*y)+G(1,1);
            v          = (vec(2)*x+vec(1)*y)+G(1,2);
            %ellip(k0)  = plot(u,v,'.','MarkerSize',10);
            ellip(k0)  = plot(u,v,'.','MarkerSize',4);
            set(ellip(k0),'Color',Coulg(1,:));
            %set(ellip(k0),'LineWidth',1.5)
            set(ellip(k0),'LineWidth',1)
            
            %         bary(k0)   = text(G(1,1),G(1,2),num2str(IndClass(k,1)),'FontSize',8.5,'Fontweight','Demi','FontName',ecrit,'Interprete',Inter,'Color',Coulg(1,:));
            bary(k0)   = text(G(1,1),G(1,2),num2str(IndClass(n_Class,1)),'FontSize',8.5,'Fontweight','Demi','FontName',ecrit,'Interprete',Inter,'Color',Coulg(1,:));
            
            ikam(l,1)  = text(X(:,1),X(:,2),num2str(Cl2),'FontSize',ncize,'Fontweight','Demi','FontName',ecrit,'Interprete',Inter,'Color',Coulg(1,:));
            
            %         mx(k)      = min(u);Mx(k) = max(u);my(k) = min(v);My(k) = max(v);
            %         mx(k)      = mx(k)-(Mx(k)-mx(k))*margin;
            %         my(k)      = my(k)-(My(k)-my(k))*margin;
            %         Mx(k)      = Mx(k)+(Mx(k)-mx(k))*margin;
            %         My(k)      = My(k)+(My(k)-my(k))*margin;
            
            mx(n_Class)      = min(u);
            Mx(n_Class) = max(u);
            my(n_Class) = min(v);
            My(n_Class) = max(v);
            
            mx(n_Class)      = mx(n_Class)-(Mx(n_Class)-mx(n_Class))*margin;
            my(n_Class)      = my(n_Class)-(My(n_Class)-my(n_Class))*margin;
            Mx(n_Class)      = Mx(n_Class)+(Mx(n_Class)-mx(n_Class))*margin;
            My(n_Class)      = My(n_Class)+(My(n_Class)-my(n_Class))*margin;
        end
        
    end
end

% end;

mx = min([mx,minx]);
Mx = max([Mx,maxx]);
my = min([my,miny]);
My = max([My,maxy]);

axis([mx Mx my My]);
line([mx Mx],[0 0]);
line([0 0],[my My]);
