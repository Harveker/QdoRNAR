/* 
Franco Endrigo
Contact: 
franco.endrigo.r@gmail.com
francoendrigo@ufpr.br
Date of creation: 30/09/2021-2021
*/
//main program
proteinFree =[]
ycalc=[]
for i = 1:6
    [y] = ode(Yz,0,t,list(diferentialSolver,pz,VrSTR(i)))
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
yfull= [yfull; [+y; +y2]]          // that is, in the matrix odd line #'s are repressor and even line #'s are fluorescence concentration

j=j+3
end

//csvWrite(yfull', "data.csv", ";",",")
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
