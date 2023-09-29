clear
clf

hList = [0.1,0.25,0.5,0.8,1];
hold on
for i=1:5
    h = hList(i);
    x=-4:h:4;
    m=length(x);
    k = 1;
   
    D=-30*eye(m);
    D=D+16*diag(ones(m-1,1),+1);
    D=D+16*diag(ones(m-1,1),-1);
    D=D-diag(ones(m-2,1), +2);
    D=D-diag(ones(m-2,1), -2);
    D=D/(12*h^2);
    D=-0.5*D; % factor for schrodinger equation
    V=0.5*k*diag(x.^2);
        
    A=D+V;
    
    [vmat, emat] = eig(A);
    plot(x, vmat(:,1), 'r-')
    error = abs(emat(1,1)-0.5);
    E(i) = error;
end
hold off
clf
figure(2)
plot(hList,E, 'ro-')



