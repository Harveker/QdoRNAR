/* 
Franco Endrigo
Contact: franco.endrigo.r@gmail.com
Date of creation: 12/10/2021,2021
*/

//data import

//description of system
/*
dMi/dt = function of volume
d[R]total/dt= MIr-Kdr*[R]total-[R]total/V*dv/dt
d[F]/dt = kf[DNA]-Kdf[F]-[F]/F*dV/dt
[DNA]=[DNA]tot -[R][DNA]tot/Kd+[R]=[DNA]tot(1-[R]/Kd+[R])

For a NAR system the Concentration of the repressor R is:
a*[R]^3+b*[R]^2+c*[R]+d = 0
where the constants are:
a=𝐾_𝑑𝑄+[𝑄]
b=𝐾_𝑑𝑄(𝐾_𝑑𝑁𝐴𝑅+𝐾_𝑑𝐵𝑀𝑔+[𝑅_𝑔]_𝑇𝑂𝑇+[𝐵𝑀_𝑔]_𝑇𝑂𝑇+[𝑅_𝑝]_𝑇𝑂𝑇)+[𝑄]*(𝐾_𝑑𝑁𝐴𝑅+𝐾_𝑑𝐵𝑀𝑔)
c=𝐾_𝑑𝑄(𝐾_𝑑𝑁𝐴𝑅*(1+[𝑅]_𝑇𝑂𝑇+[𝐵𝑀_𝑔]_𝑇𝑂𝑇))+(𝐾_𝑑𝐺𝑀𝑔*(1+[𝑅_𝑔]_𝑇𝑂𝑇+(𝐾_𝑑𝑁𝐴𝑅*[𝑄])/𝐾_𝑑𝑄+[𝑅_𝑝]_𝑇𝑂𝑇))
d=−[𝑅_𝑝]_𝑇𝑂𝑇*𝐾_𝑑𝑄*𝐾_𝑑𝑁𝐴𝑅*𝐾_𝑑𝐵𝑀𝑔

Since it is a cubic function, we need the Cardano-Tartaglia Method, and therefore:
𝑥_1=(−𝑏)/3𝑎+∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
𝑥_2= (−𝑏)/3𝑎+(−1/2+ⅈ √3/2) ∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +(−1/2+ⅈ √3/2) ∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
𝑥_3=(−𝑏)/3𝑎+(−1/2−ⅈ √3/2) ∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +(−1/2+ⅈ √3/2) ∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
Where 𝑝=𝑐/𝑎−𝑏^3/(3𝑎^2) and 𝑞=𝑑/𝑎−𝑏𝑐/(3𝑎^2 )+(2𝑏^3)/(27𝑎^3 )

Parameters	
v_r	Velocity of production of repressor(__strength of the promoter__) (mol repressor L-1 min-1)
KDR	First order rate constant for the degradation of repressor (mol L-1) (mol L-1) min-1 
KDgfp	First order rate constant for the degradation of GFP (mol L-1) (mol L-1) min-1
[Q]	concentration of free quercetin in the cell (mol L-1)
[R_g]_TOT	Total concentration of Repressor gene in the cell (mol L-1)
[B_Mg]_TOT	Total concentration of Biomarker gene in the cell (mol L-1)
[R_p]_TOT	Total concentration of Repressor gene in the cell (mol L-1)
K_dQ	Dissociation constant for the “Q~R” complex
K_dNAR	Dissociation constant for the “Promoter~R” complex
K_dBMg	Dissociation constant for the "Biomarker~R" complex
mi	constant for volume increase
kfDNA	constant of formation of dna

Differential variables	
V	cell volume (L)
[R]TOT	total concentration of repressor protein in the cell (mol L-1) i.e. in all forms
[F]	concentration of GFP in the cell (mol L-1)
t	time (min)

*/
//functions
function dy = fct(t,y,cte,VrSTR)

    //redundant to simplify programming
    V = y(1)
    Rtot = y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    cte=abs(cte)
    KfDNA=cte(1)
    kdr=cte(2)
    kdgfp=cte(3)
    Q=cte(4)
    DNAtot=cte(5)
    KQ=cte(6)
    KDNA=cte(7)
    mi=cte(8) //mi(v)
    vr=VrSTR(:) //promoter strength
    
    k= 0.25    //carrying capacity for the ricker model.
    IGR= 2 //Intrinsic Growth Rate, still testing this one.
    //Ext= 1.2    //extinction rate. Not on use yet.
    a =1+Q/KQ
    b= KDNA*a+DNAtot-Rtot
    c=-Rtot*KDNA
    d=
    Del=((b^2)-4*a*c)^0.5
    R = (-b+Del)/(2*a)    //GFP DNA concentration
    DNA = DNAtot*(1-R/(KDNA+R))    //freeDNA concentration at time t
    disp(R)
    disp(DNA, "DNA")
    //dy(1)= mi*V    // Volume over time eq
    //dy(1)= mi*V*(1-V/2)    //Logistic map
    //dy(1)= V*%e^(IGR*(1-V/k))    //Ricker model
    dy(1)=V*%e^(1-(V/k))    //Ricker model modified works to show a "overshoot" pattern but still
    //doesn't simulate the stabilization of the system well
    dy(2)= vr-kdr*Rtot-(Rtot/V)*dy(1)    //total repressor over time eq
    // this if statement is to make sure the y(3) can be call either as a matrix or
    //as a single variable.
    if length(y)==3 then
        dy(3)= KfDNA*DNA - kdgfp*F-F/V*dy(1)    //fluorescence over time
    else
        for i = 1:length(F)
            //here we call F as an vector, since we accept that F is a submatrix of y(),
            //containing multiple values of Fluorescence.
            dy(2+i)= KfDNA*DNA - kdgfp*F(i)-F(i)/V*dy(1)    //fluorescence over time
        end
    end

    disp("dy s",dy(:)) //when we call fminsearch it minimizes resulting in complex numbers, 
    //which the ODE doesn't handle well. Trying to fix that.
endfunction

function par = APR(pz,Yzero,time,VrSTR)
    par=0
    ycalc=[]
    for i = 1:6
        y = ode(Yzero,0,time,list(fct,pz,VrSTR(i)))
        ycalc= [ycalc; +y(:,:)]
    end
    yfull=[]
    j=1
    for i = 1:6
        y1=(ycalc(j+1,:)./ycalc(j,:))
        y2=(ycalc(j+2,:)./ycalc(j,:))
        yfull= [yfull; [+y1; +y2]]
        j=j+3
    end
    
    for i=1:length(time)
        for j= 1:6
        param(j)=(yfull(j*2-1,i)./4)-(yfull(j*2,i))
        par= par + sum(param)
        end
    end
    disp(param,"param")
endfunction

function SDQ =obj(pz,Yzero,time,FluorExp,VrSTR)
    SDQ=0
    y_calc= ode(Yzero,0,time,list(fct,pz,VrSTR)) //generate 6 row of results
    for i = 1:length(time)
        SSDQ(1) = (y_calc(3,i) - FluorExp(1,i))^2
        //SSDQ(2) = (y_calc(4,i) - FluorExp(2,i))^2
        //SSDQ(3) = (y_calc(5,i) - FluorExp(3,i))^2
        //SSDQ(4) = (y_calc(6,i) - FluorExp(4,i))^2
        //SSDQ(5) = (y_calc(7,i) - FluorExp(5,i))^2
        //SSDQ(6) = (y_calc(8,i) - FluorExp(6,i))^2
        SDQ = SDQ + sum(SSDQ)
    end
    disp("final sum",SDQ)
endfunction

//return
disp("function working!")



