%% simulate Poisson Jump Diffusion
mu=0.5; sigma=0.3; lambda=2; S0=100;

a=0; b=0.5;

T=1; n=1000; dt=T/n;

X=log(S0)*ones(1,n);
J=ones(1,n);
for t=2:n
    Z=normrnd(0,1);
    Nt=pois_process(lambda,dt);
    
    logY=randn(1,Nt);
    logY=a+b*logY;
    M=sum(logY(1:Nt));
    
    X(t)=X(t-1)+(mu-.5*sigma^2)*dt+sigma*sqrt(dt)*Z+M;
    J(t)=J(t-1)*exp(M);
end
S=exp(X);

%%
close all
figure

subplot(2,1,1)
t = linspace(0,T,n);
plot(t,S);

subplot(2,1,2)
plot(t,J);
ylim([0,1.2])

saveas(gcf,'jumpDiffusionSim','epsc')


%% generate Poisson process
function Nt= pois_process(lambda,t)
    p=exp(-lambda*t);
    F=p;
    n=0;

    U=rand(1,1);

    while U>F
        n=n+1;
        p=p*lambda*t/n;
        F=F+p; 
    end

    Nt=n;
end