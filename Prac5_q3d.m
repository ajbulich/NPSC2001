clear
clf

h=0.25;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

potential = 10*ones(m,1);
potential(abs(x)<1)=0;
V=diag(potential);

A=D+V;

[vmat,emat] = eig(A);
B = inv(A);

%bar(eig(A))


i = 1;
for omega=1:0.01:1.2
    %disp(omega)
    psi = ones(m,1);
    for j=1:10
        %plot(x,psi,'ro-')
        %title(j-1)
    
        old = psi;
        psi=psi*(1-omega) + omega*B*psi;
        psi = psi/max(abs(psi));
        error(j) = sum(abs(old./psi-1));
    end
    mError(i) = min(error);
    i = i + 1;
end

semilogy(1:0.01:1.2,mError, 'ro-')
ylabel("Error")
xlabel("Omega")


