N=2; %order of the Gaussian quadrature
[w,ptGaussRef]=gaussValues2DQuad(N);
% Draw Gauss points in the reference quadrilateral
plotRectangle([-1 -1],[1,-1],[1,1],[-1,1]);
hold on;
plot(ptGaussRef(:,1),ptGaussRef(:,2),'ro');
hold off;
%