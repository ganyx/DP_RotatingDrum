function [CC,Radius]=CircleThru3Dots(A,B,C)
Ah=A*A';
Bh=B*B';
Ch=C*C';
CC=zeros(size(A));
G=(C(2)-B(2))*A(1)+(A(2)-C(2))*B(1)+(B(2)-A(2))*C(1);
CC(1)=((Bh-Ch)*A(2)+(Ch-Ah)*B(2)+(Ah-Bh)*C(2))/(2*G);
CC(2)=-((Bh-Ch)*A(1)+(Ch-Ah)*B(1)+(Ah-Bh)*C(1))/(2*G);
Radius=sqrt((A-CC)*(A-CC)');
theta=linspace(0,2*pi,101);
x=CC(1)+Radius*cos(theta);
y=CC(2)+Radius*sin(theta);
plot(x,y,'r-')
ABC=[A;B;C];
hold on
plot(ABC(:,1),ABC(:,2),'b.','markersize',20)
plot(CC(1),CC(2),'r.','markersize',20)
grid on
box off
axis equal
end