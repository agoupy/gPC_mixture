%% Plot densities
clear all; 
close all;
opt_plot=1;

delta = 20;
m = 2;


switch 'cloud'
    
    case 'cloud'
        %generate sample
        n = 1e2;
        a = randn(n,1);
        b = randn(n,1);
        w=ones(n,1)./n;
        p=[a b];
        
        sample = example(delta,a,b);
        id = kmeans(sample,m);
        
    case 'quad'
        order = 10;
        dim=2;
        [p,w]= make_quadrature(dim, order);
        n=size(p,1);
        a=p(:,1);
        b=p(:,2);
        sample = example(delta,a,b);
        id = kmeans(sample,m);
end

%plot
n_new=1e3;

if opt_plot
    f1=figure(1);
    
    xcoord=linspace(-50,50,round(n_new/5));
    % [f1,xi1] = ksdensity(sample,xcoord,'Weights',w);
    
    [f1,xi1] = ksdensity(example(delta,randn(n_new,1),randn(n_new,1)),xcoord);    
    p1(1) = plot(xi1,f1,'LineWidth',2,'DisplayName','Original Points');
    
    f2=figure(2);
    p2(1)=plot3(a,b,sample,'r+','DisplayName','Original Points');
    %     hold on
    %     plot3(a(id==1),b(id==1),sample(id==1),'bo');
    %     plot3(a(id==2),b(id==2),sample(id==2),'go');
    
end
%% Build metamodel of indicator using kriging

new_p=randn(n_new,2);
[mean_id,std_id]=krigeage(new_p,p,id);
new_id=round(mean_id);

f3=figure(3);
p3(1)=plot3(p(:,1),p(:,2),id,'sc','LineWidth',2,'DisplayName','Indicator on original points');
hold on
% p3(3)=plot3(new_p(:,1),new_p(:,2),mean_id,'r+','LineWidth',2,'DisplayName','mean');
% p3(4)=plot3(new_p(:,1),new_p(:,2),mean_id+std_id,'g^','LineWidth',2,'DisplayName','mean+std');
% p3(5)=plot3(new_p(:,1),new_p(:,2),mean_id-std_id,'g^','LineWidth',2,'DisplayName','mean-std');
p3(2)=plot3(new_p(:,1),new_p(:,2),new_id,'bo','LineWidth',2,'DisplayName','Indicator on new points');
xlabel('$\xi_1$','FontSize',16,'Interpreter','latex')
ylabel('$\xi_2$','FontSize',16,'Interpreter','latex')
zlabel('$Ind(\xi_1,\xi_2)$','FontSize',16,'Interpreter','latex')
l3=legend(p3);
x=sort(p(:,1));
xlim([-5 5])
ylim([-5 5])
line(x,1/6-0.5*x.^2,2+0.*x,'LineWidth',2,'Color','r')
l3.Location='northwest';
title('Kriging interpolation of the indicator function')
set(gcf, 'Renderer', 'painters');

%% gPC estimation on each sample

% polynomials
order_max=5;
dim=2;
He=poly1D(order_max,'hermite-prob-norm');
alpha=multi_index(dim,order_max);
P_max=size(alpha,1);
pol=ones(n,P_max);

% polynomial evaluation
for i=1:P_max
    for j=1:dim
        pol(:,i)=polyval(He{alpha(i,j)+1},p(:,j)).*pol(:,i);
    end
end

[a_mod1,a_mod1_std,mse1] = gPC_mod_LS(sample(id==1),pol(id==1,:),w(id==1));
[a_mod2,a_mod2_std,mse2] = gPC_mod_LS(sample(id==2),pol(id==2,:),w(id==2));


%% Use indicator and gPC to generate a sample and plot density

%Compute densities
pol_gPC=ones(n_new,P_max);
%polynomial evaluation
for i=1:P_max
    for j=1:dim
        pol_gPC(:,i) = polyval(He{alpha(i,j)+1},new_p(:,j)).*pol_gPC(:,i);
    end
end

% sample_gPC1=pol_gPC(new_id==1,:)*a(:,1);
% sample_gPC2=pol_gPC(new_id==2,:)*a(:,2);

sample_gPC1=pol_gPC(new_id==1,:)*a_mod1;
sample_gPC2=pol_gPC(new_id==2,:)*a_mod2;

if opt_plot
    figure(1); hold on;
    %     xcoord1=linspace(0,80,500);
    %     xcoord2=linspace(-80,0,500);
    %     m1=mean(sample_gPC1);m2=mean(sample_gPC2);
    %     [f1,xi1] = ksdensity(sample_gPC1,xcoord1.*(m1>0)+xcoord2.*(m1<0));
    %     [f2,xi2] = ksdensity(sample_gPC2,xcoord1.*(m2>0)+xcoord2.*(m2<0));
    [f1,xi1] = ksdensity(cat(1,sample_gPC1,sample_gPC2),xcoord);
    p1(2)=plot(xi1,f1,'LineWidth',2,'DisplayName','gPC');
    %     p1(3)=plot(xi2,f2,'LineWidth',2,'DisplayName','gPC - sample 2');
    l1=legend(p1);
    l1.Location='northwest';
    title('Densities comparison')
    set(gcf, 'Renderer', 'painters');
    
    figure(2);hold on;
    p2(2)=plot3(new_p(new_id==1,1),new_p(new_id==1,2),sample_gPC1,'bo','DisplayName','gPC - sample 1');
    p2(3)=plot3(new_p(new_id==2,1),new_p(new_id==2,2),sample_gPC2,'go','DisplayName','gPC - sample 2');
    l2=legend(p2);
    l2.Location='northwest';
    title('Response Surface')
    set(gcf, 'Renderer', 'painters');
    
end