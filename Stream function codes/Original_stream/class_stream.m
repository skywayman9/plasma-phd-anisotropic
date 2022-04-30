 e = 2; g = 1;
 [x,y] = meshgrid(-20:20,-20:20);  % This makes regular grid
 u = e*x-g*y;                  % Linear velocity field
 v = g*x-e*y;
 
%  u = 0*x-2.5*cos(90-26.5651);
%  v = 0*x+2.5*sin(26.5651);
%  
%   u = 0*x-2.5*cos(0);
%  v = 0*x+2.5*sin(90);
 
 
 N=200;
dx=1/(N);
dy=1/(N);
NX=7;
NY=NX;
% NX=9;
% NY=NX-2;
% NX=12;
% NY=NX-4;

% NX=16;
% NY=NX-2;


[x,y] = meshgrid(-2:dx:2,-2:dy:2);

% Perpendicular lines
u=y;
v=x*1+0.0;

%Source
u=x;
v=y*1.0+0.0;

[x,y] = meshgrid(1:dx:2,1:dy:2);
%

u= 5*cos(45)+y*0;
v= 5*sin(45)+y*0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 45 degree case
% N=400;
% dx=1/(N);
% dy=1/(N);
% [x,y] = meshgrid(-2:dx:2,-1:dy:1);
% [x,y] = meshgrid(0:dx:2,0:dy:1);
% NX=14; %5, 8  and 11
% NY=NX;
% u = 1/sqrt(2)+y*0;
% v = 1/sqrt(2)+x*0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 26 degree case
% NX=3;
% NY=NX+1;
% 
% NX=7;
% NY=NX+2;
% 
% NX=11;
% NY=NX+3;
% 
% NX=15;
% NY=NX+4;
% 
% [x,y] = meshgrid(0:dx:4,0:dy:2);
% u = 0.8944+y*0;
% v = 0.4472+x*0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NX=9;
% NY=NX-2;
% [x,y] = meshgrid(-0.8944:dx:0.8944,-0.4472:dy:0.4472);
% u = 0.4472+x*0;
% v = 0.8944+y*0;



% % About 30 degree
% NX=6;
% NY=NX+3;
% % NX=10;
% [x,y] = meshgrid(0:dx:2,0:dy:1);
% % [x,y] = meshgrid(-2:dx:2,-1:dy:1);
% u = sqrt(1/4)+y*0;
% v = sqrt(3/2)+x*0;


 [phi,psi] = flowfun(u,v);  % Here comes the potential and streamfun.

minx = round(min(psi(:)),1);
maxx = round(max(psi(:)),1);
levels =  minx:200.0:maxx;
% NX=levels;
contour(x,y,psi,NX,'-r','LineWidth',1.5); % Contours of streamfunction
hold on;
minx = round(min(phi(:)),1);
maxx = round(max(phi(:)),1);
levels =  minx:200.0:maxx;
% NY=levels;
contour(x,y,phi,NY,'-b','LineWidth',1.5); % Contours of potential
hold on

% meshc(x,y,psi,phi,10s)/
xlabel('$x$','FontSize',22,'Interpreter','latex','LineWidth',1.5,'fontname','Times')
ylabel('$y$','FontSize',22,'Interpreter','latex','LineWidth',1.5,'fontname','Times')
set(gca,'FontSize',20)
%  quiver(u,v,'w')         % Now superimpose the velocity field
 
 
ax = gca; 
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickDir','out');
% yticks([0 25 50 75 100])
% xticks([0 50 100 150 200])
set(gca,'linewidth',2)
print('Streamfunction_exp','-depsc2','-r600');