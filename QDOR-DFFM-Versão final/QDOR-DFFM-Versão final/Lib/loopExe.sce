/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 08/08/2023
*/

clf([0,1,2])
clear
cd home
tic()
realtimeinit(1)
global pz
//data import

cd(get_absolute_file_path("loopExe.sce"))
cd ..
//disp("The present directory is: ", pwd())       //debug info
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
cd ..
cd "Lib"

//Parameters

pz(1) = 6.60    //k_fqdoR»     constant of formation of the qdoR protein(min-1) (Ross & Orlowski, 1982)
pz(2)= 11.04    //k_fgfp»      constant of formation of GFProtein(min-1) (Berset et al., 2017)
pz(3) = 3D-5    //k_dR»     First order rate constant for the degradation of repressor (min-1) (choosen)
pz(4) = 5D-4    //k_dgfp»        First order rate constant for the degradation of GFProtein (min-1) (Berset et al., 2017)
pz(5) = 2.5D-12    //[Q]»        concentration of total quercetin in the cell (mol L-1)(Said et al., 2016)
pz(6) = 166D-11    //[r]_T»      Total concentration of Repressor gene in the cell (mol L-1) (Conversion in mol/L of 1 copy/cell  )(Berset et al., 2017)
pz(7)= 166D-11    //[b]_T»     Total concentration of Biomarker gene in the cell (mol L-1)(Conversion in mol/L of 1 copy/cell  )(Berset et al., 2017)
pz(8) = 1D-4    //K_dQ»      Dissociation constant for the “Q~R” complex (M)(Hirooka & Fujita, 2011)
pz(9)= 1D-4       //K_dQ2»      Dissociation constant for the “Q²~R” complex (M)(choosen)
pz(10) = 4D-9    //K_dR»      Dissociation constant for the “Promoter~R” complex (M)(Hirooka & Fujita, 2011)
pz(11)= 2D-8    //K_dB»     Dissociation constant for the "Biomarker~R" complex (M)(choosen)
pz(12) = 0.0030    //r»       Constant for volume increase.(choosen)
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
Yz(1)=1D-15     //volume (Liters)(Kubitschek & Friske, 1986) and (Prats & De Pedro, 1989)
Yz(2)= 0//repressor
Yz(3)=0 //GFP

Yzero= [Yz(1);Yz(2);FluorExp(1:6,1);]
//y= ode(Yz,0,t,list(diferentialSolver,pz,VrSTR(1)))


//Data collection and Main execution

dataSave = struct ('Repressor',[], 'Fluorescence',[], 'Volume', [])
k=0
l=0
for i = 0.0023:0.0005:0.0103    //Changing this you will be able to change the range of the variations and loops
    l= l+1
    pz(12)= i   //by change this you can choose which variable you want to see variations
    exec ("Main.sce", [-1])
    for j=1:6
    dataSave(l).Repressor(j)= {yfull(k+j,:)}
    k=k+1
    dataSave(l).Fluorescence(j)={yfull(k+j,:)}
    end
    k=0
    dataSave(l).Volume = {ycalc(1,:)}
end



//maintaining this so i may use to understand how to concat strings to create multiple files
for i=1:11
repFile(i)='rep'+string(i)+'.csv'
end
for i=1:11
fluFile(i)='flu'+string(i)+'.csv'
end
for i=1:11
volFile(i)= 'vol'+string(i)+'.csv'
end

repMat=[]
fluMat=[]
volMat=[]
for i=1:11
    for j=1:6
    repMat=[repMat, cell2mat(dataSave(i).Repressor(j))';]
    fluMat=[fluMat, cell2mat(dataSave(i).Fluorescence(j))';]
    end
    volMat=[volMat, cell2mat(dataSave(i).Volume(1))';]
end
j=6
k=0
for i=1:11
    disp(k,j)
    csvWrite(repMat(:,k+1:j),repFile(i),';',',')
    csvWrite(fluMat(:,k+1:j),fluFile(i),';',',')
    csvWrite(volMat(:,i),volFile(i),';',',')
    j=j+6
    k=k+6

end
