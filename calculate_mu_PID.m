% calculate_mu_PID.m, V. Ziemann, 230512
  r=-(eye(2)+Areal+Breal*Bplus-Breal*Kd-Breal*Ki);  % r, Bronstein
  s=Areal+Breal*Bplus-2*Breal*Kd;                   % s
  t=Breal*Kd;                                       % t
  [V,Lam]=eig(r); 
  rr1=Lam(1,1);
  rr2=Lam(2,2);
  tmp=V'*s*V;
  ss1=tmp(1,1);
  ss2=tmp(2,2);
  tt1=t(1,1);
  tt2=t(2,2);
  rts1=roots([1,rr1,ss1,tt1]);
  rts2=conj(rts1);
  %rts2=roots([1,rr2,ss2,tt2])

  Z1=V*[rts1(1),0;0,rts2(1)]*V';
  Z2=V*[rts1(2),0;0,rts2(2)]*V';
  Z3=V*[rts1(3),0;0,rts2(3)]*V';
  z1=rts1(1); z2=rts1(2); z3=rts1(3);
  c1=[1,1,1;z1,z2,z3;z1^2,z2^2,z3^2]\[0;0;1];
  absz=abs([z1,z2,z3])  % for stability must be < 1
  C1=V*[c1(1),0;0,conj(c1(1))]*V';
  C2=V*[c1(2),0;0,conj(c1(2))]*V';
  C3=V*[c1(3),0;0,conj(c1(3))]*V';
  
  R1=C1*(Z1-eye(2));
  R2=C2*(Z2-eye(2));
  R3=C3*(Z3-eye(2));
  
  A11=inv(eye(2)-Z1*Z1');
  A12=inv(eye(2)-Z1*Z2');
  A13=inv(eye(2)-Z1*Z3');
  A21=inv(eye(2)-Z2*Z1');
  A22=inv(eye(2)-Z2*Z2');
  A23=inv(eye(2)-Z2*Z3');
  A31=inv(eye(2)-Z3*Z1');
  A32=inv(eye(2)-Z3*Z2');
  A33=inv(eye(2)-Z3*Z3');

  Evv=R1*A11*R1'+R1*A12*R2'+R1*A13*R3' ...
     +R2*A21*R1'+R2*A22*R2'+R2*A23*R3' ...
     +R3*A31*R1'+R3*A32*R2'+R3*A33*R3';

  Evpvp=R1*Z1*A11*Z1'*R1'+R1*Z1*A12*Z2'*R2'+R1*Z1*A13*Z3'*R3' ...
       +R2*Z2*A21*Z1'*R1'+R2*Z2*A22*Z2'*R2'+R2*Z2*A23*Z3'*R3' ...
       +R3*Z3*A31*Z1'*R1'+R3*Z3*A32*Z2'*R2'+R3*Z3*A33*Z3'*R3';

  Evpv=R1*Z1*A11*R1'+R1*Z1*A12*R2'+R1*Z1*A13*R3' ...
      +R2*Z2*A21*R1'+R2*Z2*A22*R2'+R2*Z2*A23*R3' ...
      +R3*Z3*A31*R1'+R3*Z3*A32*R2'+R3*Z3*A33*R3';

  Evvp=R1*A11*Z1'*R1'+R1*A12*Z2'*R2'+R1*A13*Z3'*R3' ...
         +R2*A21*Z1'*R1'+R2*A22*Z2'*R2'+R2*A23*Z3'*R3' ...
         +R3*A31*Z1'*R1'+R3*A32*Z2'*R2'+R3*A33*Z3'*R3';

  Evpy=R1*Z1*A11*Z1'^2*C1'+R1*Z1*A12*Z2'^2*C2'+R1*Z1*A13*Z3'^2*C3' ...
      +R2*Z2*A21*Z1'^2*C1'+R2*Z2*A22*Z2'^2*C2'+R2*Z2*A23*Z3'^2*C3' ...
      +R3*Z3*A31*Z1'^2*C1'+R3*Z3*A32*Z2'^2*C2'+R3*Z3*A33*Z3'^2*C3';
  Eyvp=Evpy';

  Evy=R1*A11*Z1'^2*C1'+R1*A12*Z2'^2*C2'+R1*A13*Z3'^2*C3' ...
        +R2*A21*Z1'^2*C1'+R2*A22*Z2'^2*C2'+R2*A23*Z3'^2*C3' ...
        +R3*A31*Z1'^2*C1'+R3*A32*Z2'^2*C2'+R3*A33*Z3'^2*C3';
  Eyv=Evy';

  Eyy=C1*Z1^2*A11*Z1'^2*C1'+C1*Z1^2*A12*Z2'^2*C2'+C1*Z1^2*A13*Z3'^2*C3' ...
        +C2*Z2^2*A21*Z1'^2*C1'+C2*Z2^2*A22*Z2'^2*C2'+C2*Z2^2*A23*Z3'^2*C3' ...
        +C3*Z3^2*A31*Z1'^2*C1'+C3*Z3^2*A32*Z2'^2*C2'+C3*Z3^2*A33*Z3'^2*C3';

%   EE=@(n,m)C1*Z1^n*A11*(Z1')^m*C1'+C1*Z1^n*A12*(Z2')^m*C2'+C1*Z1^n*A13*(Z3')^m*C3'+ ...
%      +C2*Z2^n*A21*(Z1')^m*C1'+C2*Z2^n*A22*(Z2')^m*C2'+C2*Z2^n*A23*(Z3')^m*C3'+ ...
%      +C3*Z3^n*A31*(Z1')^m*C1'+C3*Z3^n*A32*(Z2')^m*C2'+C3*Z3^n*A33*(Z3')^m*C3';
%   Evpvp=EE(2,2)-EE(1,2)-EE(2,1)+EE(1,1)
%   Evpv=EE(2,1)-EE(1,1)-EE(2,0)+EE(1,0)
%   Evpy=EE(2,2)-EE(1,2)
%   Evvp=Evpv'
%   Evv=EE(1,1)-EE(0,1)-EE(1,0)+EE(0,0)
%   Evy=EE(1,2)-EE(0,2)
%   Eyvp=Evpy'
%   Eyv=Evy'
%   Eyy=EE(2,2) 
  
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

  kk=zeros(2,2,2);
  kk(1,:,:)=R*Ki*[1,0;0,0];
  kk(2,:,:)=R*Ki*[0,0;1,0];
  kkb=zeros(2,2,2);
  kkb(:,:,1)=R*Ki*[1,0;0,0];
  kkb(:,:,2)=R*Ki*[0,1;0,0]; 

  tmp=tensorprod(hb,Evpvp,3,1);
  hvpvph=tensorprod(tmp,h,[3,2],[1,2]);

  tmp=tensorprod(hb,Evpv,3,1);
  hvpvj=tensorprod(tmp,jj,[3,2],[1,2]);

  tmp=tensorprod(hb,Evpy,3,1);
  hvpyk=tensorprod(tmp,kk,[3,2],[1,2]);

  tmp=tensorprod(jjb,Evvp,3,1);
  jvvph=tensorprod(tmp,h,[3,2],[1,2]);

  tmp=tensorprod(jjb,Evv,3,1);
  jvvj=tensorprod(tmp,jj,[3,2],[1,2]);

  tmp=tensorprod(jjb,Evy,3,1);
  jvyk=tensorprod(tmp,kk,[3,2],[1,2]);
 
  tmp=tensorprod(kkb,Eyvp,3,1);
  kyvph=tensorprod(tmp,h,[3,2],[1,2]);
 
  tmp=tensorprod(kkb,Eyv,3,1);
  kyvj=tensorprod(tmp,jj,[3,2],[1,2]);

  tmp=tensorprod(kkb,Eyy,3,1);
  kyyk=tensorprod(tmp,kk,[3,2],[1,2]);

  Q=hvpvph+hvpvj-hvpyk ...
   +jvvph+jvvj-jvyk ...
   -kyvph-kyvj+kyyk;

  ev=eig(Q);
  mu1=real(max(ev));
  mu2=real(min(ev));



