clear

hList = [0.1,0.25,0.5,0.8,1];
for i=1:5
    h = hList(i);
    x=-4:h:4;
    m=length(x);
    k = 1;
   
    D=2*eye(m);
    D=D-diag(ones(m-1,1),+1);
    D=D-diag(ones(m-1,1),-1);
    D=D/(2*h^2);
    V=0.5*k*diag(x.^2);
        
    A=D+V;
    
    [vmat, emat] = eig(A);
    plot(x, vmat(:,1), 'r-')
    error = abs(emat(1,1)-0.5);
    E(i) = error;
end
plot(hList,E, 'ro-')



