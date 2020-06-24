% different values of lambda
lambda = [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7];

c = [-1 -3]; %coefficients of the objective function
bi = [6 6]; 
bp = [-1 1]; %rhs perturbation direction

% minimizing objective function based on lambda
for i=1:length(lambda)
    l = lambda(i)
    
    %plotting feasible region
    x=linspace(0,10,100); 
    y=linspace(0,10,100);
    f1=(bi(1)+(l*bp(1)))-x;
    f2=((bi(2)+(l*bp(2)))/2)+.5.*x;

    for i=1:length(x)
        r(i)=min([f1(i),f2(i)]);
    end

    figure(1)
    plot(x,f1,'-.b','HandleVisibility','off');
    axis([0 10 0 10])
    hold on
    plot(x,f2,'--','HandleVisibility','off');
    rtxt=strcat('{\lambda}=  ', num2str(l));
    h3=area(x,r,'DisplayName',rtxt);
    p=l/10;
    h3.FaceColor = [0 0.25+p 0.5];

     % constraints setup
    A = [1,1;-1,2];
    b = [bi(1)+(l*bp(1)),bi(2)+(l*bp(2))];
    ub = [inf inf];
    lb = [0,0];
    Aeq=[];
    beq=[];
    
    xmin = linprog(c,A,b,Aeq,beq,lb,ub);
    z = c*xmin
    
    plot(xmin(1),xmin(2),'-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6],'HandleVisibility','off')
    
    % plotting
    xp = linspace(0,10,100);
    yp = (-c(1)/c(2))*xp + z/c(2);
    axis([0 7 0 7]);
    hold on
    txt=strcat('{\lambda}=  ', num2str(l));
    plot(xp,yp,'LineWidth',1,'HandleVisibility','off');
end

% graph labels
title('Perturbation of the Right-hand-side Vector');
xlabel('x(1)');
ylabel('x(2)');

hold off
legend show




    
    
    
    
    
    
