%plotting feasible region
x=linspace(0,10,100); 
y=linspace(0,10,100);
f1=6-x;
f2=3+.5.*x

for i=1:length(x)
    r(i)=min([f1(i),f2(i)]);
    disp('ri');
    r(i)
end

figure(1)
plot(x,f1,'-.b','DisplayName',' x(1)+x(2) \leq 6');
axis([0 10 0 10])
hold on
plot(x,f2,'--','DisplayName','-x(1)+2x(2) \leq 6');
length(x)
length(r)
h3=area(x,r,'DisplayName','Feasible Region');
h3.FaceColor = [0 0.25 0.25];

% constraints setup
A = [1,1;-1,2];
b = [6,6];
ub = [inf inf];
lb = [0,0];
Aeq=[];
beq=[];

% different values of lambda
lambda = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 4 5 6 7];

c = [-1 -3]; % cost vector
cp = [2 1]; % cost pertubation direction
solx=[];
soly=[];

% minimizing objective function based on lambda
for i=1:length(lambda)
    l = lambda(i);
    cf = [c(1)+cp(1)*l c(2)+cp(2)*l]; %coefficients of the objective function
    xmin = linprog(cf,A,b,Aeq,beq,lb,ub);
    z = cf*xmin;
    plot(xmin(1),xmin(2),'-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6],'HandleVisibility','off')
    
    % plotting
    xp = linspace(0,10,100);
    yp = (-cf(1)/cf(2))*xp + z/cf(2);
    axis([0 7 0 7]);
    hold on
    txt=strcat('{\lambda}=  ', num2str(l));
    plot(xp,yp,'LineWidth',2,'DisplayName',txt);
end

% graph labels
title('Perturbation of the Cost Vector');
xlabel('x(1)');
ylabel('x(2)');

hold off
legend show




    
    
    
    
    
    
