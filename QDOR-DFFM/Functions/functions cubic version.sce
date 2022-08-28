/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 12/10/2021,2021
*/

//data import

//description of system
/*
dMi/dt = function of volume
d[R]total/dt= MIr-Kdr*[R]total-[R]total/V*dv/dt
d[F]/dt = kf[DNA]-Kdf[F]-[F]/F*dV/dt
[DNA]=[DNA]tot -[R][DNA]tot/Kd+[R]=[DNA]tot(1-[R]/Kd+[R])

Now [R_p]_free and [BM_p]_free will have a MM equation that describes them:


For a NAR system the Concentration of the repressor R is:
a*[R]^3+b*[R]^2+c*[R]+d = 0
where the constants are:
a=𝐾_𝑑𝑄+[𝑄]
b=𝐾_𝑑𝑄(𝐾_𝑑𝑁𝐴𝑅+𝐾_𝑑𝐵𝑀𝑔+[𝑅_𝑔]_𝑇𝑂𝑇+[𝐵𝑀_𝑔]_𝑇𝑂𝑇+[𝑅_𝑝]_𝑇𝑂𝑇)+[𝑄]*(𝐾_𝑑𝑁𝐴𝑅+𝐾_𝑑𝐵𝑀𝑔)
c=𝐾_𝑑𝑄(𝐾_𝑑𝑁𝐴𝑅*(1+[𝑅]_𝑇𝑂𝑇+[𝐵𝑀_𝑔]_𝑇𝑂𝑇))+(𝐾_𝑑𝐺𝑀𝑔*(1+[𝑅_𝑔]_𝑇𝑂𝑇+(𝐾_𝑑𝑁𝐴𝑅*[𝑄])/𝐾_𝑑𝑄+[𝑅_𝑝]_𝑇𝑂𝑇))
d=−[𝑅_𝑝]_𝑇𝑂𝑇*𝐾_𝑑𝑄*𝐾_𝑑𝑁𝐴𝑅*𝐾_𝑑𝐵𝑀𝑔

Since it is a cubic function, we need the Cardano-Tartaglia Method, and therefore:
𝑥_1=(−𝑏)/3𝑎+∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
𝑥_2= (−𝑏)/3𝑎+(−1/2+ⅈ √3/2) ∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +(−1/2-ⅈ √3/2) ∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
𝑥_3=(−𝑏)/3𝑎+(−1/2−ⅈ √3/2) ∛(−𝑞/2+√(𝑞^2/4+𝑝^3/27)) +(−1/2+ⅈ √3/2) ∛((−𝑞)/2−√(𝑞^2/4+𝑝^3/27)) 
Where 𝑝=𝑐/𝑎−𝑏^3/(3𝑎^2) and 𝑞=𝑑/𝑎−𝑏𝑐/(3𝑎^2 )+(2𝑏^3)/(27𝑎^3 )

Parameters
(ceasse to exist??) v_r »       Velocity of production of repressor(__strength of the promoter__) (mol repressor L-1 min-1)
KDR»        First order rate constant for the degradation of repressor (min-1) 
KDgfp»      First order rate constant for the degradation of GFP (min-1)
Q»        Concentration of free quercetin in the cell (mol L-1)
R_gTOT»      Total concentration of Repressor gene in the cell (mol L-1)
R_gfree»      Concentration of free Repressor gene in the cell (mol L-1)
BM_gTOT»         Total concentration of Biomarker gene in the cell (mol L-1)
BM_gfree»         Concentration of free Biomarker gene in the cell (mol L-1)
R_pTOT»      Total concentration of Repressor gene in the cell (mol L-1)
R_pfree»      Concentration of free Repressor protein in the cell (mol L-1)
K_dQ»       Dissociation constant for the “Q~R” complex
K_dNAR»         Dissociation constant for the “Promoter~R” complex
K_dBMg»         Dissociation constant for the "Biomarker~R" complex
mi»         constant for volume increase
kfQdoR»         constant of formation of the qdoR protein
kfGFP»      constant of formation of GFProtein


Differential variables
V»      cell volume (L)
[R]TOT»     total concentration of repressor protein in the cell (mol L-1) i.e. in all forms
[F]»        concentration of GFP in the cell (mol L-1)
t»      time (min)

*/

//functions
function R_pfree = newtonBissection(cte,y)
    V = y(1)
    R_pTOT= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    
    cte=abs(cte)
    k_fqdoR=cte(1)
    k_fgfp=cte(2)
    K_dR=cte(3)
    K_dgfp=cte(4)
    Q = cte(5)
    R_gTOT = cte(6)
    BM_gTOT=cte(7)
    K_dQ=cte(8)
    K_dNAR=cte(9)
    K_dBMg=cte(10)
    mi=cte(11) //mi(v)
    
    //________________________________________________________   
Rt = y(2)    // total repressor
Dt = R_gTOT    // total of gene for repressor
Ft = BM_gTOT    // total of gene for GFP
Qf = Q    // free quercetin concentration - remains constant
Kd = K_dNAR     // repressor
Kf = K_dBMg     // GFP
Kq = K_dQ    // quercetina
tolerance = 0.000001


upper = Rt
lower = 0
Rb = (upper-lower)/2
del = 1
iteration = 0

while abs(del) > tolerance
    iteration = iteration + 1
    Rf = Rt-Rb
    del = Rt - Rf*(Dt/(Kd+Rf) + Ft/(Kf+Rf) + 1 + Qf/Kq)
    
//    disp(iteration)
//disp(Rb, "rb")
//    disp(del)
    
    if del < 0
        then lower = Rb
            Rb = lower + (upper-lower)/2

        else 
            upper = Rb
            Rb = Rb/2

        end
        
        
        if iteration>1000 
            then 
            disp("more than 1000 iterations")
            break
            end
end
R_pfree = Rt-Rb

//disp('final value of Rp_free = Rtot-Rbound')
//disp(R_pfree)

//__________________________________________________________
    
    
 // DAVID: Personally, I would test the program first with a very simple growth kinetic equation - maybe linear growth ... dy(3) = constant)
endfunction

function R_pfree =  cubicRoot (cte, y)
    
    V = y(1)
    R_pTOT= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    
    cte=abs(cte)
    k_fqdoR=cte(1)
    k_fgfp=cte(2)
    K_dR=cte(3)
    K_dgfp=cte(4)
    Q = cte(5)
    R_gTOT = cte(6)
    BM_gTOT=cte(7)
    K_dQ=cte(8)
    K_dNAR=cte(9)
    K_dBMg=cte(10)
    mi=cte(11) //mi(v)
    
    k= 0.25     //carrying capacity for the ricker model.
    IGR= 2      //Intrinsic Growth Rate, still testing this one.
    //Ext= 1.2    //extinction rate. Not on use yet.
    
    //description of the cubic function
    a = K_dQ+Q
    b = K_dQ*(K_dNAR+K_dBMg+R_gTOT+R_pTOT)+Q*(K_dNAR+K_dBMg)
    c = K_dQ*(K_dNAR*(1+R_gTOT+BM_gTOT))+(K_dBMg*(1+R_gTOT+(K_dNAR*Q/K_dQ)+R_pTOT))
    d = R_pTOT*K_dQ*K_dNAR*K_dBMg
   // disp ('a,b,c,d', a,b,c,d)
    //The Cardano-Tartaglia(CT) method
    p=(c/a)-(b^3)/(3*a^2)
    q=(d/a)-(b*c/3*a^2)+(2*b^3/27*a^3)
    //Some common operations in the root of the CT method
    G=((q^2)/4+(p^3)/27)^(1/3)
    H(1)=(-q/2+G)^(1/3)
    H(2)=(-q/2-G)^(1/3)
    L(1)=(-1/2+(%i*sqrt(3)/2))
    L(2)=(-1/2-(%i*sqrt(3)/2))
    //disp('G''s',G,'H''s:', H,'l''s',L)
    //Roots of the CT method
    R_pfree(1) =-b/(3*a)+H(1)+H(2)
    R_pfree(2) =-b/(3*a)+L(1)*H(1)+L(2)*H(2)
    R_pfree(3) =-b/(3*a)+L(2)*H(1)+L(1)*H(2)
    
endfunction

function dy = diferentialSolver(t,y,cte,VrSTR)
    [R_pfree]= newtonBissection(cte,y)
    //[R_pfree] = cubicRoot(cte,y)
    //disp(cubicRoot(cte,y))
    //disp( newtonBissection(cte,y))
    //redundant to simplify programming
    V = y(1)
    R_pTOT= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    //cte=abs(cte)
    k_fqdoR=cte(1)
    k_fgfp=cte(2)
    K_dR=cte(3)
    K_dgfp=cte(4)
    Q = cte(5)
    R_gTOT = cte(6)
    BM_gTOT=cte(7)
    K_dQ=cte(8)
    K_dNAR=cte(9)
    K_dBMg=cte(10)
    mi=cte(11) //mi(v)
    vr=VrSTR(:) //promoter strength
    
    k= 0.33    //carrying capacity for the ricker model.
    IGR= 2 //Intrinsic Growth Rate, still testing this one.
    //Ext= 1.2    //extinction rate. Not on use yet.
    
    R_gfree =R_gTOT/(K_dNAR+R_pfree(1))
    disp(R_gfree)    //Free repressor gene concentration at time t
    BM_gfree=BM_gTOT/(K_dBMg +R_pfree(1))
    disp(BM_gfree)
   // disp(list(["Rpfree: " + string(R_pfree)],["Rgfree: " + string(R_gfree)], ["BM_gfree: " + string(BM_gfree)]))
    //dy(1)= mi*V    // Volume over time eq
    //dy(1)= mi*V*(1-V/2)    //Logistic map
    //dy(1)= V*%e^(IGR*(1-V/k))    //Ricker model
    //dy(1)=V*%e^(IGR*(1-(V/k)))    //Ricker model modified works to show a "overshoot" pattern but still
    dy(1)= mi
    //doesn't simulate the stabilization of the system well
     dy(2)= (vr*R_gfree)-(K_dR*R_pfree(1))-((y(2)/V)*dy(1))    //total repressor over time eq
    // this if statement is to make sure the y(3) can be call either as a matrix or
    //as a single variable.
    if length(y)==3 then
        dy(3)= k_fgfp*BM_gfree - K_dgfp*F-F/V*dy(1)    //fluorescence over time
    else
        for i = 1:length(F)
            //here we call F as an vector, since we accept that F is a submatrix of y(),
            //containing multiple values of Fluorescence.
            dy(2+i)= k_fgfp*BM_gfree - K_dgfp*F(i)-F(i)/V*dy(1)    //fluorescence over time
        end
    end

    disp("dy s",dy(:)) //when we call fminsearch it minimizes resulting in complex numbers, 
    disp("first interation go as plan")//which the ODE doesn't handle well. Trying to fix that.
    
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


