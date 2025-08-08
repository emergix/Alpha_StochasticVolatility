(* ::Package:: *)

<<"MultivariateStatistics`"


NormalDis[x_]:=(Erf[x/Sqrt[2]]+1)/2


NormalDis[x0_,\[Sigma]_,x_]:=(Erf[((x-x0)/\[Sigma])/Sqrt[2]]+1)/2


NormalDens[x0_,\[Sigma]_,x_]:=E^(-((x-x0)^2/(2 \[Sigma]^2)))/(Sqrt[2 \[Pi]] \[Sigma])


NormalDis[x_]:=(Erf[x/Sqrt[2]]+1)/2


NormDens[x0_,\[Sigma]_,x_]:=E^(-((x-x0)^2/(2 \[Sigma]^2)))/(Sqrt[2 \[Pi]] \[Sigma])


NormDens[x_]:=E^(-((x)^2/2 ))/Sqrt[2 \[Pi]]


NormDens2[x_]:=If[Abs[x]>38100,0,E^(-((x)^2/2 ))/Sqrt[2 \[Pi]]]


NormalDens[x_]:=E^(-((x)^2/2 ))/Sqrt[2 \[Pi]]


NormDis[x_]:=(Erf[x/Sqrt[2]]+1)/2


NormaDisInv[x_]:=Sqrt[2] InverseErf[0,-1+2 x]


Nd[x_]:=N[(Erf[x/Sqrt[2]]+1)/2,40]


BS[s_,k_,t_,r_,d_,v_]:= If[k==0,s,s Exp[-d t]  NormalDis[(Log[(s Exp[r t-d t])/k]+(v^2 t)/2)/(v Sqrt[t])]-k Exp[-r t] NormalDis[(Log[(s Exp[r t-d t])/k]-(v^2 t)/2)/(v Sqrt[t])]]


BS[s_,k_,t_,r_,v_]:=If[k==0,s, s NormalDis[(Log[(s Exp[r t])/k]+(v^2 t)/2)/(v Sqrt[t])]-k Exp[-r t] NormalDis[(Log[(s Exp[r t])/k]-(v^2 t)/2)/(v Sqrt[t])]]


BS[f_,k_,t_,v_]:= If[k==0,f,f  NormalDis[(Log[f/k]+(v^2 t)/2)/(v Sqrt[t])]-k  NormalDis[(Log[f/k]-(v^2 t)/2)/(v Sqrt[t])]]


BSPut[f_,k_,v_]:=BS[f,k,v]-(f-k)


BSPut[f_,k_,t_,v_]:= BS[f,k,t,v]-(f-k)


BSDelta[f_,k_,t_,v_]:=  NormalDis[(Log[f/k]+(v^2 t)/2)/(v Sqrt[t])]


BSPutDelta[f_,k_,t_,v_]:=-NormalDis[-((Log[f/k]+(v^2 t)/2)/(v Sqrt[t]))]


BSVega[f_,k_,t_,v_]:= f  NormalDens[(Log[f/k]+(v^2 t)/2)/(v Sqrt[t])] Sqrt[t]


BSGamma[f_,k_,t_,v_]:= 1/(f v Sqrt[t])  NormalDens[(Log[f/k]+(v^2 t)/2)/(v Sqrt[t])]


BSDigitalPut[f_,k_,v_]:=1-NormalDis[(Log[f/k]-v^2/2)/v ]


BSDigitalCall[f_,k_,v_]:=NormalDis[(Log[f/k]-v^2/2)/v ]


BSDigitalCall[f_,k_,t_,v_]:=NormalDis[(Log[f/k]-(v^2 t)/2)/(v Sqrt[t])]


BS1[f_,k_,v_]:= f  NormalDis[(Log[f/k]+v^2/2)/v ]-k  NormalDis[(Log[f/k]-v^2/2)/v ]


PhiVol[m_,u_]:=Exp[m] Nd[m/u+u/2]-Nd[m/u-u/2]


ImpVolBS[F_,K_,t_,discount_,opt_]:=Module[{o=opt/(K  discount),m=Log[F/K],ss},
ss=FindRoot[PhiVol[m,u]==o,{u,0.2,0.0001,7.},AccuracyGoal->8,WorkingPrecision->30,MaxIterations->200];ss[[1,2]]/Sqrt[t]]


ImpVolBS[F_,K_,t_,opt_]:=Module[{o=opt/K ,m=Log[F/K],ss},
If[Im[F]!=0  ,Print["ImpVolBS: Im[F]\[NotEqual]0 ! F=",F ]];
If[ F<0 ,Print["ImpVolBS: F<0! F=",F ]];
If[Im[K]!=0 ,Print["ImpVolBS: Im[K]\[NotEqual]0 ! K=",K ]];
If[K<0,Print["ImpVolBS: K<0! K=",K ]];
If[o>=Exp[m],Print["ImpVolBS: opt too big! opt=",opt, " limit= ",K Exp[m]];7,
ss=FindRoot[PhiVol[m,u]==o,{u,0.2,0.0001,7.},AccuracyGoal->8,WorkingPrecision->30,MaxIterations->200];ss[[1,2]]/Sqrt[t]]]


ImpVolBSPut[F_,K_,t_,opt_]:=Module[{o=(opt+F-K)/K ,m=Log[F/K],ss},
If[Im[F]!=0  ,Print["ImpVolBS: Im[F]\[NotEqual]0 ! F=",F ]];
If[ F<0 ,Print["ImpVolBS: F<0! F=",F ]];
If[Im[K]!=0 ,Print["ImpVolBS: Im[K]\[NotEqual]0 ! K=",K ]];
If[K<0,Print["ImpVolBS: K<0! K=",K ]];
ss=FindRoot[PhiVol[m,u]==o,{u,0.2,0.0001,7.},AccuracyGoal->8,WorkingPrecision->30,MaxIterations->200];ss[[1,2]]/Sqrt[t]]


ImpVolBS[F_,List[A__],t_,discount_,opt_]:=Table[ImpVolBS[F,List[A][[i]],t,discount,opt[[i]]],{i,1,Length[List[A]]}]



ImpVolBS[F_,List[A__],t_,opt_]:=Table[ImpVolBS[F,List[A][[i]],t,opt[[i]]],{i,1,Length[List[A]]}]



ImpVolBSPut[F_,List[A__],t_,opt_]:=Table[ImpVolBSPut[F,List[A][[i]],t,opt[[i]]],{i,1,Length[List[A]]}]


AsymptoticNormDisDiff[x_,y_]:=(-Sign[x]2NormDens2[x])/(Abs[x]+Sqrt[x^2+4])+(Sign[y]2NormDens2[y])/(Abs[y]+Sqrt[y^2+4])


AsymptoticNd[x_]:=1-E^(-x^2/2)/(Sqrt[\[Pi]](x/Sqrt[2]+Sqrt[x^2/2+2]))


AsymptoticN[x_]:=1-E^(-x^2/2)/(Sqrt[\[Pi]]((2x)/Sqrt[2]))


AsymptoticNdm[x_]:=E^(-x^2/2)/(Sqrt[\[Pi]](-x/Sqrt[2]+Sqrt[x^2/2+2]))


AsymptoticNd2[x_]:=If[x>0,AsymptoticNd[x],AsymptoticNdm[x]]


AsymptoticPhiVol[m_,u_]:=If[m>=0,Exp[m] AsymptoticNd[m/u+u/2]-AsymptoticNd[m/u-u/2],Exp[m] AsymptoticNdm[m/u+u/2]-AsymptoticNdm[m/u-u/2]]


PhiVolInitial[m_,u_]:=Exp[m] NormalDis[m/u+u/2]-NormalDis[m/u-u/2]


PhiVol[m_,u_,p_]:=If[Abs[m/u]>5,AsymptoticPhiVol2[SetPrecision[m,p],SetPrecision[u,p]],PhiVolInitial[SetPrecision[m,p],SetPrecision[u,p]]]


AsymptoticPhiVol2[m_,u_]:=If[m>=0,-1+E^m+2 E^(-((-2 m+u^2)^2/(8 u^2))) Sqrt[2/\[Pi]] u* (1/(2 m+u *(-u+Sqrt[16-4 m+(4 m^2)/u^2+u^2]))-1/(2 m+u *(u+Sqrt[16+4 m+(4 m^2)/u^2+u^2]))),2 E^(-((-2 m+u^2)^2/(8 u^2))) Sqrt[2/\[Pi]] u *(1/(2 m-u* (u+Sqrt[16-4 m+(4 m^2)/u^2+u^2]))-1/(2 m+u* (u-Sqrt[16+4 m+(4 m^2)/u^2+u^2])))]


AsymptoticPhiVol3[m_,u_]:=If[m>=0,-1+E^m+2 E^(m/2-m^2/(2 u^2)) Sqrt[2/\[Pi]]  (u^3)/(4 m^2), 2 E^(m/2-m^2/(2 u^2)) Sqrt[2/\[Pi]]  u^3/(4 m^2)]


BS[f_,k_,t_,v_]:= If[k==0,f,f  NormalDis[(Log[f/k]+(v^2 t)/2)/(v Sqrt[t])]-k  NormalDis[(Log[f/k]-(v^2 t)/2)/(v Sqrt[t])]]


BS[f_,k_,t_,v_]:=If[k==0,f,Module[{m=Log[f/k],u=v Sqrt[t]},k PhiVol[m,u]]]


PhiVol[m_,u_]:=Exp[m] Nd[m/u+u/2]-Nd[m/u-u/2]


AsymptoticBS[f_,k_,T_,\[Sigma]_]:=If[f>k,(f-k)+k (Sqrt[f/(2\[Pi] k)] E^(-(Log[f/k]^2/(2 \[Sigma]^2 T))) ( \[Sigma]^3 T^(3/2))/Log[f/k]^2),k (2Sqrt[f/(2\[Pi] k)] E^(-(Log[f/k]^2/(2 \[Sigma]^2 T))) ( \[Sigma]^3 T^(3/2))/Log[f/k]^2)]


AsymptoticBSTimeValue[f_,k_,T_,\[Sigma]_]:=k (Sqrt[f/(2\[Pi] k)] E^(-(Log[f/k]^2/(2 \[Sigma]^2 T))) ( \[Sigma]^3 T^(3/2))/Log[f/k]^2)


AsymptoticBSVolImplicit[f_,k_,T_,C_]:=If[f>k,Log[f/k]/(Sqrt[3T] Sqrt[ProductLog[Log[f/k]^(2/3)/(3 ((C+k-f)^(2/3) ((2\[Pi])/(f k))^(1/3)))]]),Log[k/f]/(Sqrt[3T] Sqrt[ProductLog[Log[k/f]^(2/3)/(3 ((C)^(2/3) ((2\[Pi])/(f k))^(1/3)))]])]


AsymptoticPhiVolInverse[m_,y_]:=Module[{x},If[m>=0,x=((y+1-E^m)/((2 E^(m/2))/(4 m^2) Sqrt[2/\[Pi]]));m/(Sqrt[3] Sqrt[ProductLog[m^2/(3 x^(2/3))]]),
x=y/((2 E^(m/2))/(4 m^2) Sqrt[2/\[Pi]]);m/(Sqrt[3] Sqrt[ProductLog[m^2/(3 x^(2/3))]])]]


AsymptoticPhiVolSumToDo[m_,u_]:={-1+E^m,2 E^(m/2-m^2/(2 u^2)) Sqrt[2/\[Pi]]  (u^3)/(4 m^2)}


LogNormalDens[x0_,\[Sigma]0_,x_]:=If[x<=0,0,Exp[-1/2 (Log[x/x0]/\[Sigma]0)^2]/(x \[Sigma]0 Sqrt[2\[Pi]])]


LogNormalDistrib[x0_,\[Sigma]0_,x_]:=If[x<=0,0,(1-Erf[-Log[x/x0]/(Sqrt[2] \[Sigma]0)])/2]


LogNormalQuantile[x0_,\[Sigma]0_,y_]:=E^(-Sqrt[2] \[Sigma]0 InverseErf[0,1-2 y]) x0


BSp[f_,k_,t_,v_,bond_,p_]:= Module[{m=SetPrecision[f,p]/SetPrecision[k,p],V=SetPrecision[v,p] Sqrt[SetPrecision[t,p]]},(SetPrecision[f,p] NormalDis[(Log[m]+V^2/2)/V]-SetPrecision[k,p]  NormalDis[(Log[m]-V^2/2)/V])SetPrecision[bond,p]]


MertonCall[F_,K_,t_,\[Sigma]_,\[Lambda]_,muJ_,\[Sigma]J_,Nbterms_]:=Module[
{\[Lambda]p,i,r},
\[Lambda]p=\[Lambda] Exp[muJ+\[Sigma]J^2/2];
Sum[
r=Exp[i Log[\[Lambda] t]]/Gamma[i+1];
r  BS[F Exp[(\[Lambda] -\[Lambda]p)t +i (muJ+(\[Sigma]J^2 )/2)],K,t,Sqrt[\[Sigma]^2+i \[Sigma]J^2/t]],{i,0,Nbterms}] E^(-\[Lambda] t)]
		


MertonPut[F_,K_,t_,\[Sigma]_,\[Lambda]_,muJ_,\[Sigma]J_,Nbterms_]:=Module[
{\[Lambda]p,i,r},
\[Lambda]p=\[Lambda] Exp[muJ+\[Sigma]J^2/2];
Sum[
r=Exp[i Log[\[Lambda] t]]/Gamma[i+1];
r  (BS[F Exp[(\[Lambda] -\[Lambda]p)t +i (muJ+(\[Sigma]J^2 )/2)],K,t,Sqrt[\[Sigma]^2+i \[Sigma]J^2/t]]-(F-K)),{i,0,Nbterms}] E^(-\[Lambda] t)]
		


MertonCallPercentJ[F_,K_,t_,\[Sigma]T_,\[Lambda]_,\[Mu]J_,pJ_,Nbterms_]:=Module[
{\[Lambda]p,i,\[Sigma],\[Sigma]p},
\[Sigma]=Sqrt[(1-pJ)]\[Sigma]T;
\[Sigma]J=Sqrt[pJ/\[Lambda]]\[Sigma]T;
\[Lambda]p=E^(\[Mu]J+\[Sigma]J^2/2) \[Lambda];
Sum[(t \[Lambda])^i/Gamma[1+i]  BS[E^(t (\[Lambda]-\[Lambda]p)+i (\[Mu]J+\[Sigma]J^2/2)) F,K,t,Sqrt[\[Sigma]^2+(i \[Sigma]J^2)/t]],{i,0,Nbterms}] E^(-\[Lambda] t)]
		


Erf1[x_]:=Module[{a1=0.0705230784,a2=0.0422820123,a3=0.0092705272,a4=0.001520143,a5=0.0002765672,a6=0.000430638},1-1/(1+a1 x+a2 x^2+a3 x^3+a4 x^4+a5 x^5+a6 x^6)^16]


Erf2[x_]:=If[x>0,Erf1[x],-Erf1[-x]]


Erf3[x_]:=2/Sqrt[\[Pi]] (x-x^3/3+x^5/(5 2!)-x^7/(7 3!)+x^9/(9 4!)-x^11/(11 5!)+x^13/(13 6!)-x^15/(15 7!)+x^17/(17 8!)-x^19/(19 9!)+x^21/(21 10!)-x^23/(23 11!)+x^25/(25 12!)-x^27/(27 13!)+x^29/(29 14!)-x^31/(31 15!)+x^33/(33 16!)-x^35/(35 17!))


Erf4[x_,n_]:=2/Sqrt[\[Pi]] (x-x^3/3+Sum[(-1)^i x^(2i+1)/((2i+1)i!),{i,2,n}])


ErfDefinitif[x_]:=If[Abs[x]<4,Erf4[x,62],Erf2[x]]


BivariateDensity[x_,y_,\[Rho]_]:=E^((-x^2-y^2+2 x y \[Rho])/(2 (1-\[Rho]^2)))/(2 \[Pi] Sqrt[1- \[Rho]^2])


BivariateDensity[x0_,\[Sigma]0_,x1_,\[Sigma]1_,\[Rho]_,x_,y_]:=BivariateDensity[(x-x0)/\[Sigma]0,(y-x1)/\[Sigma]1,\[Rho]]


BivariateCummulative2D[x_,Y_,\[Rho]_]:=Module[{y},NIntegrate[ (E^(-(y^2/2))  (1+ Erf[(x-y \[Rho])/(Sqrt[2] Sqrt[1-\[Rho]^2])]))/(2 Sqrt[2\[Pi]]),{y,-Infinity ,Y}]]


BivariateCummulative2D[S1_,\[Sigma]1_,S2_,\[Sigma]2_,\[Rho]_,x_,Y_]:=Module[{y},NIntegrate[ (E^(-(y^2/2))  (1+ Erf[(((x-S1)/\[Sigma]1)-y \[Rho])/(Sqrt[2] Sqrt[1-\[Rho]^2])]))/(2 Sqrt[2\[Pi]]),{y,-Infinity ,((Y-S2)/\[Sigma]2)}]]


DyBivariateCummulative2D[x_,y_,\[Rho]_]:= (E^(-(y^2/2))  (1+ Erf[(x-y \[Rho])/(Sqrt[2] Sqrt[1-\[Rho]^2])]))/(2 Sqrt[2\[Pi]])


DyBivariateCummulative2D[S1_,\[Sigma]1_,S2_,\[Sigma]2_,\[Rho]_,K1_,K2_]:= 1/\[Sigma]2 DyBivariateCummulative2D[(K1-S1)/\[Sigma]1,(K2-S2)/\[Sigma]2,\[Rho]]


DxBivariateCummulative2D[x_,y_,\[Rho]_]:=DyBivariateCummulative2D[y,x,\[Rho]]


DxBivariateCummulative2D[S1_,\[Sigma]1_,S2_,\[Sigma]2_,\[Rho]_,K1_,K2_]:= 1/\[Sigma]1 DxBivariateCummulative2D[(K1-S1)/\[Sigma]1,(K2-S2)/\[Sigma]2,\[Rho]]


BivariateCummulative2D[x_,Y_,\[Rho]_]:=Module[{y},NIntegrate[ (E^(-(y^2/2))  NormalDis[(x-y \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]],{y,-Infinity ,Y}]]


BiLogCummulative2D[x1_,\[Sigma]1_,x2_,\[Sigma]2_,\[Rho]_,X1_,X2_]:=If[(X1<=0)||(X2<=0),0,BivariateCummulative2D[Log[X1/x1]/\[Sigma]1,Log[X2/x2]/\[Sigma]2,\[Rho]]]


DyBiLogCummulative2D[x1_,\[Sigma]1_,x2_,\[Sigma]2_,\[Rho]_,X1_,X2_]:=If[(X1<=0)||(X2<=0),0,DyBivariateCummulative2D[Log[X1/x1]/\[Sigma]1,Log[X2/x2]/\[Sigma]2,\[Rho]]/(\[Sigma]2 X2)]


BivariateCummulative2DY[x_,Y_,\[Rho]_]:=Module[{y},NIntegrate[ (y E^(-(y^2/2))  NormalDis[(x-y \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]],{y,-Infinity ,Y}]]


BivariateCummulative2DX[X_,y_,\[Rho]_]:=Module[{x},NIntegrate[ (x E^(-(x^2/2)) NormalDis[(y-x \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]] ,{x,-Infinity ,X}]]


BivariateCummulative2DYY[x_,Y_,\[Rho]_]:=Module[{y},NIntegrate[ (y^2 E^(-(y^2/2))  NormalDis[(x-y \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]],{y,-Infinity ,Y}]]


BivariateCummulative2DXX[X_,y_,\[Rho]_]:=Module[{x},NIntegrate[ (x^2 E^(-(x^2/2)) NormalDis[(y-x \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]] ,{x,-Infinity ,X}]]


BivariateCummulative2DXY[x_,y_,\[Rho]_]:=(E^(-(x^2/2)) (1-\[Rho]^2))/Sqrt[2 \[Pi]] ( NormalDens[(y-x \[Rho])/Sqrt[1- \[Rho]^2]] Sqrt[1-\[Rho]^2]- x \[Rho]  NormalDis[(y-x \[Rho])/Sqrt[1- \[Rho]^2]])+\[Rho] BivariateCummulative2DYY[x,y,\[Rho]]


BivariateCummulative2DXYapp[X_,Y_,\[Rho]_]:=NIntegrate[(x -\[Rho] y)y  BivariateDensity[x,y,\[Rho]],{x,-Infinity,X},{y,-Infinity,Y}]+\[Rho] BivariateCummulative2DYY[X,Y,\[Rho]]


DxBivariateCummulative2D[x_,y_,\[Rho]_]:=(E^(-(x^2/2)) NormalDis[(y-x \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]] 


DyBivariateCummulative2D[x_,y_,\[Rho]_]:=(E^(-(y^2/2)) NormalDis[(x-y \[Rho])/ Sqrt[1-\[Rho]^2]])/Sqrt[2\[Pi]] 


BivariateCopula[x_,y_,\[Rho]_]:=BivariateCummulative2D[NormaDisInv[x],NormaDisInv[y],\[Rho]]


DxBivariateCopula[x_,y_,\[Rho]_]:=If[x<0.00001,0,DxBivariateCummulative2D[NormaDisInv[x],NormaDisInv[y],\[Rho]]/NormalDens[NormaDisInv[x]]]


DyBivariateCopula[x_,y_,\[Rho]_]:=If[y<0.00001,0,DyBivariateCummulative2D[NormaDisInv[x],NormaDisInv[y],\[Rho]]/NormalDens[NormaDisInv[y]]]


DxyBivariateCopula[x_,y_,\[Rho]_]:=BivariateDensity[NormaDisInv[x],NormaDisInv[y],\[Rho]]/(NormalDens[NormaDisInv[x]]NormalDens[NormaDisInv[y]])


BivariateCummulative2DApprox[x_,Y_,\[Rho]_,legendreCoefs_]:=Module[{ccc=BetweenNormalCoeffs[legendreCoefs,-6,Y,4]},CoeffBasedIntegrate[(   (1+ Erf[(x-# \[Rho])/(Sqrt[2] Sqrt[1-\[Rho]^2])])/2 )&,ccc]];


NormalDis[x_]:=(Erf[x/Sqrt[2]]+1)/2


NormalCall[F_,\[Sigma]_,K_,T_]:=(F-K)NormDis[(F-K)/(\[Sigma] Sqrt[T])]+\[Sigma] Sqrt[T] NormDens[(F-K)/(\[Sigma] Sqrt[T])]


NormalCallVega[F_,K_,T_,\[Sigma]_]:=(E^(-((F-K)^2/(2 T \[Sigma]^2))) Sqrt[T])/Sqrt[2 \[Pi]]


NormalImplicitVol[F_,K_,T_,o_]:=
Module[{ss,u,msg},If[o<F-K,msg=Join["Intrinseq vale>option! Val. Intr.",ToString[F-K]," opt.=",ToString[o]];
Message[NormalImplicitVol::OPTLESSVI,msg];Throw[msg],
ss=FindRoot[NormalCall[F,u,K,T]==o,{u,0.2,0.0001,15},AccuracyGoal->8,MaxIterations->200];u /. ss[[1]]]]


ImpliedNormalCorrelation[Volnormal1_,volonormal2_,volnormalspread_]:=(volnormalspread^2-Volnormal1^2-volonormal2^2)/(-2Volnormal1 volonormal2)


NormalCallStd[\[Sigma]inv_]:=NormDis[\[Sigma]inv]+ NormDens[\[Sigma]inv]/\[Sigma]inv


NormalImplicitVol1[F_,K_,T_,o_]:=(F-K)/(NormalImplicitVolStd[o/(F-K)] Sqrt[T])


NormalImplicitCorrelation[\[Sigma]1_,\[Sigma]2_,vol_]:=(\[Sigma]1^2+\[Sigma]2^2-vol^2)/(2\[Sigma]1 \[Sigma]2)


NewtonFindRoot[func_,initial_,shift_,accuracy_,nbloopmax_]:=
Module[{derf,x=initial,i=0,fx},
fx=func[x];
While[(Abs[fx]>accuracy)&&(i<nbloopmax),
derf=(func[x+shift]-fx)/shift;
x-=fx/derf;
i++;
fx=func[x];
];
{x,i,Abs[fx]}]


NewtonFindMinimum[func_,initial_,shift_,accuracy_,nbloopmax_,printflag_]:=
Module[{der2f,derf,x=initial,i=0,fx,fx1,fx2},
fx2=func[x-shift];
fx1=func[x+shift];
derf=Re[(fx1-fx2)/(2shift)];
If[printflag>0,Print["NewtonFindMinimum initial: x=",x," derf=",derf]];
While[(Abs[derf]>accuracy)&&(i<nbloopmax),
fx=func[x];
der2f=(fx1+fx2-2fx)/(shift*shift);
x-=Re[derf/der2f];
If[printflag>0,Print["NewtonFindMinimum: x=",x," fx=",fx," derf=",derf," der2f=",der2f]];
i++;
fx2=func[x-shift];
fx1=func[x+shift];
derf=(fx1-fx2)/(2shift);
];
{x,i,Abs[derf]}]



