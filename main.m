%% Plot densities
clear all; close all;

%generate sample
n=1e3;
a=randn(n,1);
b=randn(n,1);
points=[a b];
w=ones(n,1)./n;

m=2;
delta=20;
sample = example(delta,a,b);
id = kmeans(sample,m);
sample_1 = sample(id==1); 
sample_2 = sample(id==2);
points1=[a(id==1) b(id==1)];
points2=[a(id==2) b(id==2)];
sample_separated={sample_1,sample_2};

%plot
if 0
f1=figure;
plot3(a,b,sample,'r+');
hold on
h(1)=plot3(a(id==1),b(id==1),sample(id==1),'bo');
h(2)=plot3(a(id==2),b(id==2),sample(id==2),'go');

f2=figure;
xcoord=linspace(-50,50,1000);
[f1,xi1] = ksdensity(sample,xcoord,'Weights',w);
plot(xi1,f1,'LineWidth',2,'DisplayName','Reference');
end
%% Build metamodel of indicator using kriging

new_points=randn(1e2,2);
[m,std]=krigeage(new_points,points,id);
new_id=round(m);

f3=figure;
plot3(points(:,1),points(:,2),id,'.k','LineWidth',2,'DisplayName','Original Points')
hold on
plot3(new_points(:,1),new_points(:,2),m,'r+','LineWidth',2,'DisplayName','mean')
plot3(new_points(:,1),new_points(:,2),m+std,'g^','LineWidth',2,'DisplayName','mean+std')
plot3(new_points(:,1),new_points(:,2),m-std,'g^','LineWidth',2,'DisplayName','mean-std')
plot3(new_points(:,1),new_points(:,2),new_id,'bo','LineWidth',2,'DisplayName','Interpoled indicatrix')
legend show
xlabel('$\xi_1$','FontSize',16,'Interpreter','latex')
ylabel('$\xi_2$','FontSize',16,'Interpreter','latex')
zlabel('$Ind(\xi_1,\xi_2)$','FontSize',16,'Interpreter','latex')
title('Kriging interpolation of the indicator function')


%% Simple gPC estimation

% generate quadrature
order_max = 100;
dim=2;
[ points, weights ] = make_quadrature( dim, order_max );
Nq=size(points,1);

% estimate variable at quadrature points
sample_quad=example(delta,points(:,1),points(:,2));

% generate & evaluate polynomials
He=poly1D(order_max,'hermite-prob-norm');
alpha=multi_index(dim,order_max);
P_max=size(alpha,1);
pol=ones(Nq,P_max);
for i=1:P_max
    %polynomial evaluation
    for j=1:dim
        pol(:,i)=polyval(He{alpha(i,j)+1},points(:,j)).*pol(:,i);
    end
end

% Compute the gPC coefficients
coeffs=NaN(P_max,1);
for i=1:P_max
    coeffs(i,1)=sum(sample_quad.*pol(:,i).*weights);
end

% Use the gPC expansion to approx density
points_gPC=[a b];

for order=[5 15 40 80]
    temp=(order+1):(order+dim);
    P=prod(temp)/factorial(dim);
    pol_MC=ones(n,P);
    for i=1:P
        %polynomial evaluation
        for j=1:dim
            pol_MC(:,i)=polyval(He{alpha(i,j)+1},points_gPC(:,j)).*pol_MC(:,i);
        end
    end
    
    sample_gPC=sum(repmat(coeffs(1:P),[1 n]).*pol_MC.',1);
    [f2,xi2] = ksdensity(sample_gPC);
    plot(xi2,f2,'--','LineWidth',2,'DisplayName',['p = ' num2str(order)]);
end

legend show
xlim([-6 6])
ylim([0 0.35])