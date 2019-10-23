%Euler Method:  Apoptosis Intrinsic Pathway

%reaction rates (1/s)
k1 = 2*10^(-4) ; 
k2 = 2*10^(-4);
k3= 2.8*10^(-7);
k4= 2.8*10^(-4);
k5 =2*10^(-5);
k6 =2*10^(-4); %Bcl + tbid forward rate
k8 =2*10^(-3); %Bcl + Bax forward rate


%initial conditions (nM)
tbid= [30,30,30,30] ;
Bax= [80,80,80,80];
tbidB= [0,0,0,0];
Bax2= [0,0,0,0] ;
cytc= [100,100,100,100];
apaf= [100,100,100,100];
apop = [0,0,0,0];
P9 = [20,20,20,20];
C9 = [0,0,0,0];
P3 = [200,200,200,200];
C3 = [0,0,0,0];
Bcl2 = [30,45,59,90]; %inhibitor
tBcl = [0,0,0,0];
BBcl = [0,0,0,0] ;


dt = 1;
tgraph = zeros(1,20001);
C9graph = zeros(1,20001);
P9graph = zeros(1,20001);
P3graph = zeros(1,20001);
C3graph = zeros(1,20001);
Bax2graph = zeros(1,20001);
cytgraph = zeros(1,20001);
Bclgraph = zeros(1,20001);
for i=1:4

for t = 0:20000
    if Bax2(i)<=20
    
        tbid(i) = tbid(i) + (dt*(-k1*tbid(i)*Bax(i) + k1*tbidB(i)*Bax(i) -k6*tbid(i)*Bcl2(i)));
        Bax(i) = Bax(i) + (dt*(-k1*tbid(i)*Bax(i) -k1*tbidB(i)*Bax(i) -k8*Bcl2(i)*Bax(i) ));
        tbidB(i) = tbidB(i) +(dt*(k1*tbid(i)*Bax(i) -k1*tbid(i)*Bax(i) ));
        Bax2(i) = Bax2(i) + (dt*(k1*tbid(i)*Bax(i) ));
        Bcl2(i) = Bcl2(i) + (dt*(-k6*Bcl2(i)*tbid(i) -k8*Bcl2(i)*Bax(i)));
        tBcl(i) =  tBcl(i) + (dt*(k6*Bcl2(i)*tbid(i)));
        BBcl(i) = BBcl(i) + (dt*(k8*Bcl2(i)*Bax(i)));
        
    else
        tbid(i) = tbid(i) + (dt*(-k1*tbid(i)*Bax(i) + k1*tbidB(i)*Bax(i) -k6*tbid(i)*Bcl2(i)));
        Bax(i) = Bax(i) + (dt*(-k1*tbid(i)*Bax(i) -k1*tbidB(i)*Bax(i) -k8*Bcl2(i)*Bax(i) ));
        tbidB(i) = tbidB(i) +(dt*(k1*tbid(i)*Bax(i) -k1*tbid(i)*Bax(i) ));
        Bax2(i) = Bax2(i) + (dt*(k1*tbid(i)*Bax(i) ));
        Bcl2(i) = Bcl2(i) + (dt*(-k6*Bcl2(i)*tbid(i) -k8*Bcl2(i)*Bax(i)));
        tBcl(i) =  tBcl(i) + (dt*(k6*Bcl2(i)*tbid(i)));
        BBcl(i) = BBcl(i) + (dt*(k8*Bcl2(i)*Bax(i)));
        
        cytc(i) = cytc(i) + (dt*(-k3*cytc(i)*apaf(i)));
        apaf(i) = apaf(i) + (dt*(-k3*cytc(i)*apaf(i)));
        apop(i) = apop(i) + (dt*(k3*cytc(i)*apaf(i)));
        P9(i) = P9(i) + (dt*(-k4*apop(i)*P9(i)));
        C9(i) = C9(i) + (dt*(k4*apop(i)*P9(i)));
        C3(i) = C3(i) + (dt*(k5*C9(i)*P3(i)));
        P3(i) = P3(i) + (dt*(-k5*C9(i)*P3(i)));    
    end
        
        
C9graph(i,t+1) = C9(i);
P9graph(i,t+1) = P9(i);
P3graph(i,t+1) = P3(i);
C3graph(i,t+1) = C3(i);
cytgraph(i,t+1) = cytc(i);
Bclgraph(i,t+1) =Bcl2(i);
Bax2graph(i,t+1) = Bax2(i);
tgraph(i,t+1) = t;

end
end
figure
x=plot(C3graph(1,:))
hold on;
y=plot(C3graph(2,:))
hold on;
z=plot(C3graph(3,:))
hold on;
m=plot(C3graph(4,:))
xlabel('time(sec)')
ylabel('Caspase 3 concentration (nM)')
set([x y z m],'LineWidth',1.5)
xlim([0 14000])
ylim([-10 200])
legend('[Bcl2] = baseline','[Bcl2] = + 50%','[Bcl2] = + 97%','[Bcl2] = + 200%')


