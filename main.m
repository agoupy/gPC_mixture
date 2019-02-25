%% Plot densities
%generate sample
n=100000;
a=randn(n,1);
b=randn(n,1);

delta=2;
sample=example(delta,a,b);

%plot
close(figure(1))
fi=figure(1);
fi.Position=[1336 417 707 552];
[f1,xi1] = ksdensity(sample);
plot(xi1,f1,'LineWidth',2,'DisplayName','Reference');
hold on
ax=gca;
ax.FontSize=16;
set(gcf, 'Renderer', 'painters');
title(['Theoretical density vs gPC reconstruction for $\delta$ = ' num2str(delta)], ...
    'Interpreter', 'latex');

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