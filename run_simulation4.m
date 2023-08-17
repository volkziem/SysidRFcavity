% run_sumulation.m, V. Ziemann, 230619
% x=[Vreal;Vimag];   % actually the differences to steady-state
% u=[Ireal,Iimag];
clear all;  close all
Niter=1000000;
Nforget=20e4;
alpha=1-1/Nforget;
%alpha=1;
R=1;            % shunt impedance
Z=0.5;          % weight of current
sig=0.01;       % process noise level
dt=1e-7;        % sample time at rate 10 MHz
omega12=2e5;    %=1/tau=filling time=1/10 microseconds
domega=-1e5;    %=cavity detuning =10kHz
xset=0;         % Voltage set point
q0=[omega12*dt,domega*dt];        % bandwidth and detuning
F0=[-q0(1),-q0(2);q0(2),-q0(1)];  % eq. 17
Kp=(1-sqrt(1+(R/Z)^2))/R;  % feedback gain (from Riccati equation)
P=eye(2);                   % initial value of PT
Areal=eye(2)+F0;            % eq. 15
Breal=R*omega12*dt*eye(2);  % defined just after eq. eq. 15
Bplus=Kp*eye(2);            % feedback matrix for P-regulator
qhat=zeros(2,1);            % inintial parameter estimate

% convergence analysis for P-controller
Csys=Areal+Breal*Bplus             % A+Kp*B 
[V,Lam]=eig(Csys*Csys');           % eq. 35
lambda=Lam(1,1)
h=zeros(2,2,2);
h(1,:,:)=[-1+R*Kp,0;0,1];          % eq. 29
h(2,:,:)=[0,-1;-1+R*Kp,0];

hb=zeros(2,2,2); 
hb(:,:,1)=[-1+R*Kp,0;0,1];         % eq. 31
hb(:,:,2)=[0,-1+R*Kp;-1,0];

HH=tensorprod(hb,h,[3,2],[1,2])    % eq. 38
Q_first=HH/(1-lambda)

mu1=HH(1,1)/(1-lambda);            % eq. 39
mu2=HH(2,2)/(1-lambda);
mu_P=[mu1,mu2]  % mu for proportional controller only!

Kd=0;   
Ki=1;

%calculate_mu_PD, mu_PD=[mu1,mu2]              % appendix A
%calculate_mu_PI, mu_PI=[mu1,mu2]              % appendix B
calculate_mu_PID, mu_PID=real([mu1,mu2])      % appendix C

x=sig*randn(2,1);     % initialize cavity voltage reading, process noise
data=zeros(Niter,8);  % storage for later plotting
xprev=zeros(size(x)); % previous voltage for PD regulator
xint=zeros(size(x));  % initialzie integrator for PI regulator
tic
for iter=1:Niter      % main iteration loop
  if 0   % pulsed operation
    if iter==Niter/4, xset=[5;0]; end
    if iter==Niter/2, xset=[0;0]; end
  end
  if 0   % step reponse
    if iter==Niter/4
      q0=2*q0;
      Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
      Breal=Breal*2;
    end
    if iter==Niter/2
      q0=0.5*q0;
      Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
      Breal=Breal*0.5;
    end
  end
  if 0 % periodic perturbation
    q0(2)=0.01*sin(10*6.2832e-6*iter);  
    Areal=eye(2)+[-q0(1),-q0(2);q0(2),-q0(1)];
  end
  x=x-xset;                                 % adjust feedback set point
  xint=xint+x;                              % update integrator for PI
  u=Bplus*x-Kd*(x-xprev)-Ki*xint;           % feedback
  xnew=Areal*x+Breal*u+sig*randn(size(x));  % eq. 15
  %...................................system identification
  y=xnew-x;                                       % eq. 20, right
  G=[-x(1)+u(1)*R,-x(2);-x(2)+u(2)*R,x(1)];       % eq. 20, left
  tmp2=eye(2)-P*G'*inv(alpha*eye(2)+G*P*G')*G;    % eq. 53 (or 25)
  Pnew=tmp2*P/alpha;                              % rest of eq. 53
  qhat=tmp2*(qhat+P*G'*y/alpha);                  % eq. 54 (or 26)
  xprev=x;
  x=xnew;
  P=Pnew;
  %................................save for later plotting
  data(iter,1)=x(1);
  data(iter,2)=x(2);
  data(iter,3)=u(1);
  data(iter,4)=u(2);
  data(iter,5)=qhat(1);
  data(iter,6)=qhat(2);
  data(iter,7)=P(1,1);
  data(iter,8)=P(2,2);
end
toc
mm=1:Niter;
mm2=Niter/2:Niter;   % to calculate averages over second half
sigma=[std(data(:,1)),std(data(:,3))]   % rms voltage and current, real part 

if 0    % plot of voltages and currents
  figure(1)
  subplot(2,1,1);
  plot(mm,data(:,1),'k',mm,data(:,2),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('v_r, v_i');
  ylim([-4,10.2])
  set(gca,'FontSize',16);
  subplot(2,1,2);
  plot(mm,data(:,3),'k',mm,data(:,4),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('i_r, i_i');
  set(gca,'FontSize',16);
end

if 0   % evolution of fit parameters
  figure(2); clf
%  subplot(2,1,1);
  plot(mm,data(:,5),'k',mm,data(:,6),'r','LineWidth',2);
  xlabel('Iterations'); ylabel('q(1), q(2)');
  legend('q(1)=\omega_{12}dt','q(2)=\Delta\omega dt');
  set(gca,'FontSize',16);
  pavg=[rms(data(mm2,5)),rms(data(mm2,6))] 
end

if 1   % convergence plots of estimation error and PT
  if alpha==1    % without forgetting, assumes p0=1
    pt=@(mu,T)1./(1+sig^2*mu*T);
    gt1=@(mu,T)1./(1+sig^2*mu*T);
  else           % with forgetting
    b=(1-alpha)/alpha;
    pt=@(mu,T) b./(mu*sig^2/alpha^2+(b-mu*sig^2/alpha^2)*exp(-b*T)); 
    gamma1=mu1*sig^2/(b*alpha^2);
    c1=1-1/gamma1;
    gt11=@(T)(((1-c1)*exp(-b*T))./(1-c1*exp(-b*T))).^alpha;
    gamma2=mu2*sig^2/(b*alpha^2);
    c2=1-1/gamma2;
    gt12=@(T)(((1-c2)*exp(-b*T))./(1-c2*exp(-b*T))).^alpha;
  end
  mmm=logspace(0,log10(Niter),30); % for plotting asterisks

  figure(3); clf                  % Plot estimation error
  eb1=sig*sqrt(data(:,7));        % empirical error bar of q(1)
  data(:,1)=abs(data(:,5)-q0(1)); % estimation error of q(1)
  if alpha==1
    loglog(mm,data(:,1),'k',mm,eb1,'b-.',mmm,data(1,1)*gt1(mu1,mmm),'r*','LineWidth',2);
  else
    loglog(mm,data(:,1),'k',mm,eb1,'b-.',mmm,data(1,1)*gt11(mmm),'r*','LineWidth',2);
  end
  xlim([1,Niter]); ylim([9e-7,1.3*data(1,1)])
  xlabel('Iterations'); ylabel('Estimation error |a_T(1)|')
  legend('Simulation','Errorbars','Model');
  set(gca,'FontSize',16);
  set(gcf,'Position',[30,500,900,500])
  
  figure(4); clf %   plot PT
  loglog(mm,data(:,7),'k',mmm,pt(mu1,mmm),'r*',  ...
         mm,data(:,8),'k--',mmm,pt(mu2,mmm),'b*','LineWidth',2);
  xlim([1,Niter]); ylim([5e-4,1.3])
  xlabel('Iterations'); ylabel('P_T')
  legend('Simulation P_T(1,1)','Model with \mu_1', ...
    'Simulation P_T(2,2)','Model with \mu_2');
  set(gca,'FontSize',16);
  set(gcf,'Position',[30,30,900,500])
end

%save('tmp.mat',"data")
