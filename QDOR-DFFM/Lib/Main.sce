/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 30/09/2021-2021
*/
clc
clear
cd home
tic()
//data import

cd(get_absolute_file_path("Main.sce"))
cd ..
disp(["The present directory is: ",pwd()])       //debug info
cd "Input"
data = csvRead("input.csv",';',',','string')
datadouble = strtod(data(3:31,3:$),'.') //declaring a variable containing only the doubles of the archive.
n=5
for i = 1:5
    datadouble(n,:)= [] //excluding some lines in the csv archive, to simplify further data processing.
    n=(n+5)-1
end

clear n i

//description of system
/*
dMi/dt = function of volume
d[R]total/dt= MIr-Kdr*[R]total-[R]total/V*dv/dt
d[F]/dt = kf[DNA]-Kdf[F]-[F]/F*dV/dt
[DNA]=[DNA]tot -[R][DNA]tot/Kd+[R]=[DNA]tot(1-[R]/Kd+[R])
[R]= (-(Kd(1+[Q]/Kq)+[DNA]tot-[R]tot)+sqrroot((Kd(1+[Q]/Kq)+[DNA]tot-[R]tot)²-4(1+[Q]/Kq)[R]tot*Kd)/2(1+[Q]/Kq)

Parameters
vr = velocity of production of repressor (mol repressor L-1 min-1)
k_dr = first order rate constant for the degradation of repressor (mol L-1) (mol L-1) min-1 
k_dr = first order rate constant for the degradation of GFP (mol L-1) (mol L-1) min-1
[Q] = concentration of free quercetin in the cell (mol L-1)
[DNA]TOT = concentration of DNA in the cell (mol L-1)
KQ = dissociation constant for the complex “Q~R”
KD = dissociation constant for the complex “DNA~R”

Differential variables	
V = cell volume (L)
[R]TOT = total concentration of repressor protein in the cell (mol L-1) i.e. in all forms
[F] = concentration of GFP in the cell (mol L-1)
t = t (min)
Are the various assumptions about mathematical relationships correct?

GFP needs oxygen in order for the fluorescence to occur – so you can penalize GFP fluorescence with some function related to the total volume of biomass

*/
//functions
cd ..
cd "Functions"
exec("functions cubic version.sce",[-1])
//Parameters
//pFinal of 25/11
//pz1=[5.856452627863837, 0.9910359014797071, 0.015122725353707884, 5.133119899353982E-4, 18.82280793257661, 3.0727121407212827, 4.390882279546192, 0.04169975746516301]

pz(1) = 0.005    //k_fqdoR»     constant of formation of the qdoR protein
pz(2)= 0.012    //k_fgfp»      constant of formation of GFProtein(mol/min-1)
pz(3) = 0.1    //K_dR»     First order rate constant for the degradation of repressor (min-1) 
pz(4) = 0.00035    //K_dgfp»        First order rate constant for the degradation of GFProtein (min-1)
pz(5) = 100D-9    //[Q]»        concentration of free quercetin in the cell (mol L-1)
pz(6) = 12D-9    //[R_g]_TOT»      Total concentration of Repressor gene in the cell (mol L-1)
pz(7)= 12D-9    //[BM_g]_TOT»     Total concentration of Biomarker gene in the cell (mol L-1)
pz(8) = 0.7    //K_dQ»      Dissociation constant for the “Q~R” complex
pz(9) = 0.04    //K_dNAR»      Dissociation constant for the “Promoter~R” complex
pz(10)= 0.8    //K_dBMg»     Dissociation constant for the "Biomarker~R" complex
pz(11) = 0.001    //mi»       Constant for volume increase
[VrSTR]= [0.1; 0.15; 0.16; 0.24; 0.33; 0.45]*0.2    //vr = velocity of production of repressor (mol repressor L-1 min-1)

//Graph Parameters
//here it will open more graph windows
f1=scf(0)
//f2=scf(1)
f3=scf(2) //not using yet

//Data
t = (strtod(data(2,3:$),'.')).*60
QPG114 = datadouble(1:4,:)
QPG115 = datadouble(5:8,:)
QPG116 = datadouble(9:12,:)
QPG105 = datadouble(13:16,:)
QPG110 = datadouble(17:20,:)
QPG = datadouble(21:24,:)
FluorExp = datadouble(3:4:$,:)

//initial conditions
Yz(1)=2
Yz(2)=0
Yz(3)=0

Yzero= [Yz(1);Yz(2);FluorExp(1:6,1);]
//test= ode(Yz,0,t,list(fct,pz,VrSTR(1)))


//main program
ycalc=[]
for i = 1:6
    y = ode(Yz,0,t,list(diferentialSolver,pz,VrSTR(i)))
    ycalc= [ycalc; +y(:,:)]
end

//this was already working at 26/10/2021

///////////////////////////////////
/////until here is works just fine,//////
/////if you uncomment the ode///////
////and comment the fminsearch/////
//////////////////////////////////
//pFinal = fminsearch(list(APR,Yzero,t,VrSTR),pz) 
/*
PFF=[] //parfinalfull
for i = 1:6
    pFinal = fminsearch(list(obj,Yzero,t,FluorExp,VrSTR(i)),pz) 
    //here i'm testing to see if i can adjust first the parameters(pz) for each experiment, 
    //without the complex numbers issue, after fixing that, the next step is to adjust the pz's to all
    //experiments simultaneously.
    PFF = [PFF; +pFinal(:,:)]
end
*/


//plot of the math model of the experiment

scf(f1)
yfull=[]
j=1
for i = 1:6
y=(ycalc(j+1,:)./ycalc(j,:))
y2=(ycalc(j+2,:)./ycalc(j,:))
yfull= [yfull; [+y; +y2]]
j=j+3
end
j=1
titulos=["QPG114","QPG115","QPG116","QPG105","QPG110","QPG100"]

da=gda() // obtendo manipulador dos eixo modelos para ver e editar campos
// título padrão
da.title.text="QdoR NAR"
da.title.foreground = 12;
da.title.font_size = 4;
// rótulos x padrões
da.x_label.text="Tempo(mins)";
da.x_label.font_style = 3;
da.x_label.font_size = 3;
da.x_label.foreground = 5;
da.x_location = "bottom";
// rótulos y padrões
da.y_label.text="M (mol/L)";
da.y_label.font_style = 3;
da.y_label.font_size = 3;
da.y_label.font_angle = 0;
da.y_label.foreground = 3;
da.y_location = "left";
da.thickness = 2;
da.foreground = 1;

if yfull(2,:)==yfull(4,:) then
    plot(t',[yfull(1,:);yfull(2,:)]'), title("QdoR Nar")
else
    for i=1:6
        subplot(2,3,i), plot(t',[yfull(j,:);yfull(j+1,:)]'), title(titulos(i))
        j=j+2
    end
end

/*

*/

toc()
/*plot of experimental data
scf(f2)
subplot(3,2,1), plot(t',[QPG114([1,3],:)]'), title("QPG114")
subplot(3,2,2), plot(t',[QPG115([1,3],:)]'), title("QPG115")
subplot(3,2,3), plot(t',[QPG116([1,3],:)]'), title("QPG116")
subplot(3,2,4), plot(t',[QPG105([1,3],:)]'), title("QPG105")
subplot(3,2,5), plot(t',[QPG110([1,3],:)]'), title("QPG110")
subplot(3,2,6), plot(t',[QPG([1,3],:)]'), title("QPG")
*/


//plot of the behavior of the volume based on population growth
scf(f3)
plot(t',ycalc(1,:)', 'k-'), title("Volume Behavior") 


/*[10/11] David Alexander Mitchell
//Inputting values of the parameters
vr = 1    // velocity of production of repressor (mol repressor L-1 min-1)
kf = 0.05 // rate constant for GFP production from free DNA
kdr = 0.1 // first order rate constant for the degradation of repressor (mol L-1) (mol L-1) min-1 
kdf = 0.1 // first order rate constant for the degradation of GFP (mol L-1) (mol L-1) min-1 
concQ = 5
concDNAtot = 20
KQ = 1
KD = 1
dVdt = 0.00001 // for the moment, let V increase linearly
*/
