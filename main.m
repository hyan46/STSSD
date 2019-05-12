load data.mat
Y = data;
allgamma = [1e-2 1 2 3 5 9 11 30 50 90 120 160 200 240 280 320];
%linspace(10,200,20);
lambda = [0.1 0.1 0.6];
nx = size(Y,1);
ny = size(Y,2);
nt = size(Y,3);
kx = 50; ky = 50; kt = 3;
sdx = 3; sdy = 3; sdt = 3;

%%
% B{1} = [];
% B{2} = [];
% You need specify the number of knots based on how smooth is the
% background (kx,ky) and how large is your anomaly (skx,sky)
B{1} = bsplineBasis(nx,kx,sdx);
B{2} = bsplineBasis(ny,ky,sdy);
skx = round(nx/2);
sky = round(nx/2);
Bs{1} = bsplineBasis(nx,skx,1);
Bs{2} = bsplineBasis(ny,sky,1);


[T2,Snow,Yhat,t,Itr,defect,Tte] =ewmamonit(Y,B,[],lambda,allgamma,'maxIter',3,'issave',1,'type','h');
%%
T2tr=T2;
[ mT2,sd,Ttr,Itr] = chartIC( T2tr(:,1:100));
[ Ttr,Itr] = chartOC( T2tr,mT2,sd)
L = max(Ttr(1:150));
[ Ttr,Itr] = chartOC( T2tr,mT2,sd)

%%

odx = find(Ttr>L);
n = size(Y,3)
f = @(x) log(x-min(Ttr)+1e-6);

plot(2:n,f(Ttr(2:n)),'k.-',2:n,f(L)*ones(n-1,1),'k--',odx,f(Ttr(odx)),'ro','MarkerSize',5)
%  plot(1:n,log(Tte),'k*-',1:n,log(TH)*ones(n,1),'k--')
 %save resulthard.mat
set(gca,'FontSize',14)
% title('Control Chart')
ylabel('log(Testing Statistics)')
xlabel('Time')
