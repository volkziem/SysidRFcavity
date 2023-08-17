% calculate_mu_PD.m, V. Ziemann, 230810
  Csys=Areal+Breal*Bplus-Breal*Kd;
  Dsys=Breal*Kd;
  Z1=0.5*Csys+sqrtm(0.25*Csys^2+Dsys);
  Z2=0.5*Csys-sqrtm(0.25*Csys^2+Dsys);
  RR=inv(Z1-Z2); 
  A11=inv(eye(2)-Z1*Z1');
  A12=inv(eye(2)-Z1*Z2');
  A21=inv(eye(2)-Z2*Z1');
  A22=inv(eye(2)-Z2*Z2');

  Evpvp=RR*(Z1*A11*Z1'-Z1*A12*Z2'-Z2*A21*Z1'+Z2*A22*Z2')*RR';
  Evpv=RR*(Z1*A11-Z1*A12-Z2*A21+Z2*A22)*RR';
  Evvp=RR*(A11*Z1'-A12*Z2'-A21*Z1'+A22*Z2')*RR';
  Evv=RR*(A11-A12-A21+A22)*RR';

  h=zeros(2,2,2);
  h(1,:,:)=[-1+R*(Kp-Kd),0;0,1];
  h(2,:,:)=[0,-1;-1+R*(Kp-Kd),0];
  hb=zeros(2,2,2); 
  hb(:,:,1)=[-1+R*(Kp-Kd),0;0,1];
  hb(:,:,2)=[0,-1+R*(Kp-Kd);-1,0];

  jj=zeros(2,2,2);
  jj(1,:,:)=R*Kd*[1,0;0,0];
  jj(2,:,:)=R*Kd*[0,0;1,0];
  jjb=zeros(2,2,2);
  jjb(:,:,1)=R*Kd*[1,0;0,0];
  jjb(:,:,2)=R*Kd*[0,1;0,0];

  tmp=tensorprod(hb,Evpvp,3,1);
  hvvh=tensorprod(tmp,h,[3,2],[1,2]);

  tmp=tensorprod(hb,Evpv,3,1);
  hvvj=tensorprod(tmp,jj,[3,2],[1,2]);
  tmp=tensorprod(jjb,Evvp,3,1);
  jvvh=tensorprod(tmp,h,[3,2],[1,2]);

  tmp=tensorprod(jjb,Evv,3,1);
  jvvj=tensorprod(tmp,jj,[3,2],[1,2]);
  Q=hvvh+hvvj+jvvh+jvvj
%   mu1=Q(1,1);
%   mu2=Q(2,2);
  ev=eig(Q);
  mu1=max(ev);
  mu2=min(ev);