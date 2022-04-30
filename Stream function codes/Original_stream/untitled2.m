mu=[0 0];
sigma=[5 -2 ;-2 2];
x = -10:0.1:10; y = x;
[X,Y] = meshgrid(x,y);
F = mvnpdf([X(:) Y(:)], mu,sigma);
F = reshape(F,length(x),length(x));
contour(x,y,F);
xlabel('X'); ylabel('Y');
grid on
axis square
%%DRawing line 
yi=-4:2:4 ;
xi = interp1(y,x,yi) ;
hold on
plot(xi,yi,'r')