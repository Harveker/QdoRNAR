/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 30/09/2021-2021
*/
clf([0,1,2])
clear
cd home
tic()
realtimeinit(1)
//data import

cd(get_absolute_file_path("Main.sce"))
cd ..
disp("The present directory is: ", pwd())       //debug info
cd "Input"
data = csvRead("input.csv",';',',','string')
dataOD105 = csvRead("105 qdoR PqdoI gfp.csv",";",".","string")
dataOD110 = csvRead("110 qdoR PqdoI gfp.csv",";",".","string")
dataOD114 = csvRead("114 qdoR PqdoI gfp.csv",";",".","string")
dataOD115 = csvRead("115 qdoR PqdoI gfp.csv",";",".","string")
dataOD116 = csvRead("116 qdoR PqdoI gfp.csv",";",".","string")
dataODPqdoiqdoR = csvRead("PqdoI qdoR PqdoI gfp.csv",";",".","string")  
dataODPqdorqdoR = csvRead("PqdoR qdoR PqdoI gfp.csv",";",".","string")
datadouble = strtod(data(3:31,3:$),'.') //declaring a variable containing only the doubles of the archive.


fullODData = [[dataOD105];[dataOD110];[dataOD114];[dataOD115];[dataOD116];[dataODPqdoiqdoR];[dataODPqdorqdoR]]
fullODData = strtod(fullODData(2:$,3:$),'.')

for i =1:143
    if isnan(fullODData(i,$))  then
        fullODData(i,:)=[]
    end
end

for i =2:137
    if fullODData(i,:)==fullODData(1,:) then
        fullODData(i,:)= []
    end
end

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

pz(1) = 6.60    //k_fqdoR»     constant of formation of the qdoR protein(min-1) (Ross & Orlowski, 1982)
pz(2)= 11.04    //k_fgfp»      constant of formation of GFProtein(min-1) (Berset et al., 2017)
pz(3) = 3D-5    //k_dR»     First order rate constant for the degradation of repressor (min-1) (choosen)
pz(4) = 5D-4    //k_dgfp»        First order rate constant for the degradation of GFProtein (min-1) (Berset et al., 2017)
pz(5) = 2.5D-12    //[Q]»        concentration of total quercetin in the cell (mol L-1)(Said et al., 2016)
pz(6) = 166D-11    //[r]_T»      Total concentration of Repressor gene in the cell (mol L-1) (Conversion in mol/L of 1 copy/cell  )(Berset et al., 2017)
pz(7)= 166D-11    //[b]_T»     Total concentration of Biomarker gene in the cell (mol L-1)(Conversion in mol/L of 1 copy/cell  )(Berset et al., 2017)
pz(8) = 1D-4 
//pz(8) = 1D-4    //K_dQ»      Dissociation constant for the “Q~R” complex (M)(Hirooka & Fujita, 2011)
pz(9)= 1D-4       //K_dQ2»      Dissociation constant for the “Q²~R” complex (M)(choosen)
//pz(10) = 4D-9    //K_dR»      Dissociation constant for the “Promoter~R” complex (M)(Hirooka & Fujita, 2011)
pz(10) = 9D-9 
pz(11)= 1D-9    //K_dB»     Dissociation constant for the "Biomarker~R" complex (M)(Hirooka & Fujita, 2011)
pz(12) = 0.0033    //r»       Constant for volume increase.(choosen)
pz(13) =3.5D-15     // k»     Carrying capacity. (µL) (Choosen)
pz(14) = 7.20D-2       //GFPm»       Maturation constant of gfp(min-1) (Berset et al., 2017) (not used yet)
[VrSTR]= [0.1; 0.15; 0.16; 0.24; 0.33; 0.45]    //vr = velocity of production of repressor (mol repressor L-1 min-1) (anderson promoter)

//Graph Parameters
//here it will open more graph windows
f1=scf(0) //Volume behavior
f2=scf(1) //Fluorescence Graphs
f3=scf(2) //Repressor graphs
//Data
//t = (strtod(data(2,3:$),'.')).*60
t=1:900/91:900
QPG114 = datadouble(1:4,:)
QPG115 = datadouble(5:8,:)
QPG116 = datadouble(9:12,:)
QPG105 = datadouble(13:16,:)
QPG110 = datadouble(17:20,:)
QPG = datadouble(21:24,:)
FluorExp = datadouble(3:4:$,:)

//initial conditions
Yz(1)=1D-15     //volume (Kubitschek & Friske, 1986) and (Prats & De Pedro, 1989)
Yz(2)= 72.D-12//repressor
Yz(3)=0 //GFP

Yzero= [Yz(1);Yz(2);FluorExp(1:6,1);]
//test= ode(Yz,0,t,list(fct,pz,VrSTR(1)))


//main program
ycalc=[]
for i = 1:6
    y = ode(Yz,0,t,list(diferentialSolver,pz,VrSTR(i)))
    ycalc= [ycalc; +y(:,:)]
end

ycalcExtract = ycalc'
//plot of the math model of the experiment

scf(f1)
yfull=[]
j=1
for i = 1:6
y=(ycalc(j+1,:)./ycalc(j,:))
//y= y./max(y)      //this is a normalized version, so every concentration is divided by the maximal 
y2=(ycalc(j+2,:)./ycalc(j,:))
//y2=y2./max(y2)      //this is a normalized version, so every concentration is divided by the maximal 
yfull= [yfull; [+y; +y2]]          // that is, in the matrix odd line #'s are repressor and even line #'s are fluorescence concentrations

j=j+3
end
j=1
titulos=["QPG114","QPG115","QPG116","QPG105","QPG110","QPG100"]

da=gda() // Getting axis manipulator
// Standart Title
da.title.text="modelFDFM-QdoR-GFP"
da.title.foreground = 12;
da.title.font_size = 4;
// standart x labels
da.x_label.text="Tempo(mins)";
da.x_label.font_style = 3;
da.x_label.font_size = 3;
da.x_label.foreground = 5;
da.x_location = "bottom";
// standart y labels
da.y_label.text="µM";
da.y_label.font_style = 3;
da.y_label.font_size = 3;
da.y_label.font_angle = 0;
da.y_label.foreground = 3;
da.y_location = "left";
da.thickness = 2;
da.foreground = 1;
//plot(t',[yfull(1,:);yfull(2,:)]'), title("QdoR Nar")

if yfull(2,:)==yfull(4,:) then
    plot(t',[yfull(1,:);yfull(2,:)]'), title("QdoR Nar")
    else
        for i=1:6
        scf(f1)
        //subplot(2,3,i), plot(t',yfull(j+1,:)','g-'), title(titulos(i))       //fluorescence plot
        plot(t',yfull(j+1,:)','g-')     //this will plot every biomarker curve in the same graph window
        scf(f2)
        //subplot(2,3,i), plot(t',yfull(j,:)','b-'), title(titulos(i))         //QdoR plot
        plot(t',yfull(j,:)')    //this will plot every repressor curve in the same graph window
        j=j+2
        //scf(f1)
        //subplot(2,3,i), plot(t',yfull(j+1,:),'g-', t',yfull(j,:),'b-'), title(titulos(i))
        //j=j+2
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
