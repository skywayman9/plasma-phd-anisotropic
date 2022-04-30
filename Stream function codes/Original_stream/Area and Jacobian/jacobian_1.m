F=@(x,y) 2*exp(x+y.^2).*y; %present function
xmin=0; % rectangle dimensions
xmax=1;
ymin=0;
ymax=1;
% -----------------------------------------------
% In order to create a rectangular mesh in Matlab
% we can make the following procedure
% -----------------------------------------------

numDiv=4; %number of subdivisions
hx=(xmax-xmin)/numDiv;
hy=(ymax-ymin)/numDiv;
xi=xmin:hx:xmax; %all points in the x-axis
eta=ymin:hy:ymax; %all points in the y-axis
[X,Y] = meshgrid(xi,eta); % create a grid of points
Z = F(X,Y); % function values
surf(X,Y,Z); % optional: is just to visualize the result
[elem,vertex] = surf2patch(X,Y,Z); % the main variables
numElem=size(elem,1) %total number of elements
numVert= size(vertex,1) % total number of vertices

N=2; %number of Gauss points = NxN
integralTot=0;
for i=1:numElem %compute the integral on each element
    v1=[vertex(elem(i,1),1),vertex(elem(i,1),2)];
    v2=[vertex(elem(i,2),1),vertex(elem(i,2),2)];
    v3=[vertex(elem(i,3),1),vertex(elem(i,3),2)];
    v4=[vertex(elem(i,4),1),vertex(elem(i,4),2)];
    vertices=[v1;v2;v3;v4];
    elemInteg=integQuad(F,vertices,N);
    integralTot=integralTot+elemInteg;
end
actualIntegVal= (exp(1)-1)^2 %exact known value
errorInt=abs(actualIntegVal-integralTot) %absolute error


function valInteg = integQuad(F,vertices,N)
    [w,ptGaussRef]=gaussValues2DQuad(N);
    % Shape functions
    Psi1=@(x,y)(1-x).*(1-y)/4;
    Psi2=@(x,y)(1+x).*(1-y)/4;
    Psi3=@(x,y)(1+x).*(1+y)/4;
    Psi4=@(x,y)(1-x).*(1+y)/4;
    % Shape function derivatives
    dPsi11=@(x,y) -(1-y)/4;
    dPsi21=@(x,y) (1-y)/4;
    dPsi31=@(x,y) (1+y)/4;
    dPsi41=@(x,y) -(1+y)/4;
    dPsi12=@(x,y) -(1-x)/4;
    dPsi22=@(x,y) -(1+x)/4;
    dPsi32=@(x,y) (1+x)/4;
    dPsi42=@(x,y) (1-x)/4;
    % Gradient matrix
    Jacb =@(x,y) [dPsi11(x,y), dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
                  dPsi12(x,y), dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];
    % evaluate Shape functions on Gaussian reference points
    xi = ptGaussRef(:,1);
    eta = ptGaussRef(:,2);
    evalPsi1 = Psi1(xi,eta);
    evalPsi2 = Psi2(xi,eta);
    evalPsi3 = Psi3(xi,eta);
    evalPsi4 = Psi4(xi,eta);
    % from the change of variables function
    ptGaussDomain = [evalPsi1,evalPsi2,evalPsi3,evalPsi4]*vertices;
    % evaluate Jacobian contribution for each point
    for i=1:size(xi,1)
        evalDetJacb(i) = abs(det(Jacb(xi(i),eta(i))*vertices));
    end
    %evaluate the function on the domain points
    evalF=F(ptGaussDomain(:,1),ptGaussDomain(:,2));
    % Finally, apply Gauss formula
    suma=0;
    for i=1:size(ptGaussDomain,1)
        suma=suma+w(i)*evalF(i)*evalDetJacb(i);
    end
    valInteg = suma;
end