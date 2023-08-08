/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 12/10/2021,2021
*/
//functions
global output;
function [status, yout, nextstep, proteinFree] = Munnin(t, y, flag, proteinFree,cte)
    k_fR=cte(1)
    k_fB=cte(2)
    k_dR=cte(3)
    k_dB=cte(4)
    Q = cte(5)
    r_T = cte(6)
    b_T=cte(7)
    K_dQ=cte(8)
    K_dQ2=cte(9)
    K_dR=cte(10)
    K_dB=cte(11)
    r=cte(12) //mi(v)
    k= cte(13)
    if flag == 'cont'
        // Calculate proteinFree at this time step (assuming you have the appropriate code to do this)
        R_gfree =r_T*K_dR/(K_dR+R_pf) // Calculate R_gfree
        BM_gfree=b_T*K_dB/(K_dB +R_pf) // Calculate BM_gfree
        R_pf = proteinFree(1) // Calculate R_pf
        
        // Append the new data to the proteinFree matrix
        proteinFree = [proteinFree, [R_gfree; BM_gfree; R_pf]];
    end
    
    // Set status and other output arguments
    status = 0;
    yout = [];
    nextstep = [];
endfunction

function R_pfree = newtonBissection(cte,y)
    V = y(1)
    R_T= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    
    //cte=abs(cte)
    k_fR=cte(1)
    k_fB=cte(2)
    k_dR=cte(3)
    k_dB=cte(4)
    Q = cte(5)
    r_T = cte(6)
    b_T=cte(7)
    K_dQ=cte(8)
    K_dQ2=cte(9)
    K_dR=cte(10)
    K_dB=cte(11)
    r=cte(12) //mi(v)
    k= cte(13)
    vr=VrSTR(:) //promoter strength 
    
    //________________________________________________________   
Rt = y(2)    // total repressor
Dt = r_T    // total of gene for repressor
Ft = b_T    // total of gene for GFP
Qf = Q    // free quercetin concentration - remains constant
Kd = k_dR     // repressor
Kf = K_dB     // GFP
Kq = K_dQ    // quercetina1
//Kq2 = k_dQ2    // quercetina2
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

function [TD, R_pfree] = cubicRoot (cte, y)
    
    V = y(1)
    R_T= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    
    //cte=abs(cte)
    k_fR=cte(1)
    k_fB=cte(2)
    k_dR=cte(3)
    k_dB=cte(4)
    Q = cte(5)
    r_T = cte(6)
    b_T=cte(7)
    K_dQ=cte(8)
    K_dQ2=cte(9)
    K_dR=cte(10)
    K_dB=cte(11)
    r=cte(12) //mi(v)
    k= cte(13)
    vr=VrSTR(:) //promoter strength 
    
    
    k= 0.25     //carrying capacity for the ricker model.
    IGR= 2      //Intrinsic Growth Rate, still testing this one.
    //Ext= 1.2    //extinction rate. Not on use yet.
    
    //description of the cubic function
    a = K_dQ+Q+(Q/K_dQ2)
    b = K_dQ*(R_T-K_dR-K_dB+r_T+b_T+Q*(K_dR+K_dB)/K_dQ+Q*(K_dR+K_dB)/K_dQ*K_dB)
    c = K_dQ*(R_T*K_dR+R_T*K_dB-K_dR*K_dB+K_dB*r_T+K_dR*b_T+(Q*K_dR*K_dB)/K_dQ +(Q*K_dR*K_dB)/(K_dQ^2*K_dQ2 ))
    d = K_dQ*K_dR*K_dB*R_T
    //Dp= disp ('a,b,c,d', a,b,c,d)
    //Dp= return(Dp)
    //The Cardano-Tartaglia(CT) method
    p=(c/a)-(b^2)/(3*a^2)
    q=(d/a)-(b*c/3*a^2)+(2*b^3/27*a^3)
    //Some common operations in the root of the CT method
    G=((q^2)/4+(p^3)/27)^(1/3)
    H(1)=(-q/2+G)^(1/3)
    H(2)=(-q/2-G)^(1/3)
    L(1)=(-1/2+(%i*sqrt(3)/2))
    L(2)=(-1/2-(%i*sqrt(3)/2))
    //disp('G''s',G,'H''s:', H,'l''s',L)
    //discriminant to analise the cubic equation
    Delta= 18*a*b*c*d-4*(b^3)*d+(b^2)*(c^2)-4*a*(c^3)-27*(a^2)*(d^2)
    //Roots of the CT method
    R_pfree(1) =-b/(3*a)+H(1)+H(2)
    R_pfree(2) =-b/(3*a)+L(1)*H(1)+L(2)*H(2)
    R_pfree(3) =-b/(3*a)+L(2)*H(1)+L(1)*H(2)
    
    [TD]= return [Delta;a;b;c;d;[R_pfree]]

endfunction


function [dy] = diferentialSolver(t,y,cte,VrSTR)
    //[R_pfree]= newtonBissection(cte,y)
    n=1
    [R_pfree] =(cubicRoot(cte,y))
    [cData] = [R_pfree]
    //disp(cubicRoot(cte,y))
    R_pf=[]
    
    // this will select the root in that is a real number, meaning that will always have two complex roots in this system
    if abs(imag(R_pfree(6)))<1D-10 then
        R_pf = abs(real(R_pfree(6)))
        //disp("here1")
    elseif (abs(imag(R_pfree(7))))<1.D-10 then
        R_pf = abs(real(R_pfree(7)))
        //disp("here2")
    elseif abs(imag(R_pfree(8)))<1.D-10
        R_pf = abs(real(R_pfree(8)))
        //disp("here3")
    end

    //disp( newtonBissection(cte,y))
    //redundant to simplify programming
    //[diffdata] = return [R_pfree]
    //R_pf= R_pfree
    V = y(1)
    R_T= y(2)
    // this if statement is to make sure the y(3) can be call either as a matrix or as a single variable.
    if length(y)==3 then
        F = y(3)
    else
        [F]=y(3:$)
    end
    //cte=abs(cte)
    k_fR=cte(1)
    k_fB=cte(2)
    k_dR=cte(3)
    k_dB=cte(4)
    Q = cte(5)
    r_T = cte(6)
    b_T=cte(7)
    K_dQ=cte(8)
    K_dQ2=cte(9)
    K_dR=cte(10)
    K_dB=cte(11)
    r=cte(12) //mi(v)
    k= cte(13)
    vr=VrSTR(:) //promoter strength 
    

    //IGR= 2 Intrinsic Growth Rate, still testing this one.
    //Ext= 1.2    //extinction rate. Not on use yet.
    
    R_gfree =r_T*K_dR/(K_dR+R_pf)
    //disp(R_gfree)    //Free repressor gene concentration at time t
    BM_gfree=b_T*K_dB/(K_dB +R_pf)
    
    //disp(BM_gfree)
   // disp(list(["Rpfree: " + string(R_pfree)],["Rgfree: " + string(R_gfree)], ["BM_gfree: " + string(BM_gfree)]))
    //dy(1)= V*%e^r*(1-V/k)        // Volume over time eq Ricker Model
    
    //the differential equations begins here
    dy(1)= r*V*(1-(V/k))    // Volume over time eq

    //does simulate the stabilization of the system well
     dy(2)= (k_fR*vr*R_gfree)-(k_dR*R_pf)-((y(2)/V)*dy(1))    //total repressor over time eq
    // this if statement is to make sure the y(3) can be call either as a matrix or
    //as a single variable.
    if length(y)==3 then
        dy(3)= k_fB*BM_gfree - k_dB*F-F/V*dy(1)    //fluorescence over time
        dy(4)= R_pf
        dy(5)=r_T*K_dR/(K_dR+R_pf)
        dy(6)=b_T*K_dB/(K_dB +R_pf)
        dy(7)= (Q^2)*R_pf/K_dQ2*K_dQ
        dy(8)= Q*R_pf/K_dQ
        dy(9)= R_pf*BM_gfree/K_dB
        dy(10)= R_pf*R_gfree/K_dR
    else
        for i = 1:length(F)
            //here we call F as an vector, since we accept that F is a submatrix of y(),
            //containing multiple values of Fluorescence.
            dy(2+i)= k_fB*BM_gfree - k_dB*F(i)-(F(i)/V*dy(1))    //fluorescence over time
            dy(10)=r_T*K_dR/(K_dR+R_pf)
            dy(11)=b_T*K_dB/(K_dB +R_pf)
        end
    end
    disp(VrSTR,dy(4:10))

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


