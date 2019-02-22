%% Plot densities
%generate sample
n=10000;
a=randn(n,1);
b=randn(n,1);

delta=[1 1.5 2 2.5];
n_sample=size(delta,2);
sample=nan(n_sample,n);
for i_sample=1:n_sample
    sample(i_sample,:)=example(delta(i_sample),a,b);
end

%plot
f1=figure(1);
for i_sample=1:n_sample
[f1,xi1] = ksdensity(sample(i_sample,:));
plot(xi1,f1,'DisplayName',['delta = ' num2str(delta(i_sample))]);
hold on
end
legend show
