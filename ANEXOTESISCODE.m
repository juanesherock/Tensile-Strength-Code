clc, clear, close all
testID = 'XXX';
sect1 = (0.02*0.000366);   %meters^2

ws = xlsread('UTM-Tests.xlsx',testID);
l = length(ws);     %SAVE LENGTH OF WORKSHEET
[m, tp] = deal(zeros(l(1),1));  %CREATES A ZEROS VECTOR FOR STORAGE SLOPE
e = true;

for i=10:l    %CREATES A CYCLE FOR READ ALL POSITIONS OF WS
    xo = ws(1,1);   
    xi = ws(i,1);
    yo = ws(1,2);
    yi = ws(i,2);
    m(i) = (yi - yo)/(xi - xo); %SLOPE EQUATION
    if m(i) < m(i-1)            %WHEN SLOPE STARTS TO DECREASE
        while e
            k = m(i-1);         %SAVE INITIAL M FOR COMPARISON ONE TIME
            e = false;
            c = i-1;
        end        
        tp(i) = m(i)/k;         %CALCULATES THE PERCENTAGE OF VARIATION WITH RESPECT TO K POINT
        if tp(i) < 0.99          %WHEN VARIATION WILL BE MINOR THAN SOME PERCENTAGE FINISH 
            d = i;
            break
        end
    end
end

%%REGRESSION EQUATION ESTIMATION

x = ws(c:i,1); %STROKE (mm)
y = ws(c:i,2); %FORCE (N)
format long
b1 = x\y;      %CALCULATING AVERAGED SLOPE

yCalc1 = b1*x; %INITIAL CALCULATION
X = [ones(length(x),1) x];  %GENERATES A VECTOR WITH THE SAME POSITIONS OF SLOPE EQUATION
b = X\y;                    %DIVIDES X BETWEEN y TO FIND THE INCERTECPT

yCalc2 = X*b;               %NOW MULTIPLY "ORIGINAL" X WITH b 

Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);

%%FIND THE AVERAGE STEP BETWEEN X POSITIONS

stepx = 0;
for i = 2:length(x)
    stepx = (x(i)-x(i-1)) + stepx;
end
stepx = stepx / length(x);

%%RECALCULATES TOE UNTIL Y <= 0 

toesup = zeros(c,2);        %TOE SUPRESS VECTOR
execyc = 0;
if y(1) > 0 
    xx = [1 x(1)];      %SET INITIAL POINT FOR REGRESSION    
    i = 1;              %SET COUNTER INDEX
    while true
        if i == 13
            execyc = execyc + 1; %DOES THAT REALLY WORK???----
            toesup((i-1)*(execyc+1),2) = 0;     %INCREASES THE SIZE OF TOESUP OTHER TWELVE POS
        end
        xx(2) = xx(2) - stepx;    %CALCULATE PREVIOUS X POSITION
        yy = xx * b;              %CALCULATE PREVIOUS Y POSITION
        toesup(i,1) = xx(2);
        toesup(i,2) = yy;
        if yy <= 0            
            break
        end
        %disp(i)
        %pause(0.5)
        i = i+1;            %INCREMENT COUNTER
    end
end
toesup = flip(toesup);  %MIRROR ALONG COLUMNS. BOTTOM UP

curve = [toesup; [x yCalc2]; ws(d+1:end,1:2)]; %RECONSTRUCTION OF THE WHOLE CURVE
disx_ws = dist(min(ws(:,2)),max(ws(:,2)));  %DISTANCE BETWEEN MIN AND MAX VALUE AT X (ORIGINAL DATA)
disx_fc = dist(min(curve(:,1)),max(curve(:,1)));  %DISTANCE BETWEEN MIN AND MAX VALUE AT X (FITTED DATA)

fc = disx_fc/disx_ws;   %DETERMINE THE INCREASE RATIO

%%DELETES REMAINING ZEROS AT THE START

while true
    if curve(1,2) <= 0
        curve(1,:) = []; %IF FIRST Y POS <= 0 DELETE THE WHOLE ROW, SO THE NEXT IS NOW THE 1ST ROW
    else
        break %UNTIL 
    end
end
curve(:,1) = curve(:,1) + abs(min(curve(:,1)));
curve(:,1) = curve(:,1) / fc;
curve = [0 0; curve];

%%STRESS STRAIN CURVE

%SPACE FOR SECT1 WHICH MOVES UP TO FAST THIS PROCESS
testpiecelg = 0.1;  %meters
elastic_region = curve(1:d,2)/sect1;

elongation = curve(:,1)/1000; %mm to m
sect2 = (sect1 * testpiecelg)./(testpiecelg + elongation(d+1:end,1));
plastic_region = curve(d+1:end,2)./sect2;

stress = [elastic_region ; plastic_region];

lgth = testpiecelg + elongation;
strain = (lgth - testpiecelg)/testpiecelg;

E = stress(1:d)./strain(1:d);
while true
    if isnan(E(1)) || isinf(E(1))
        E(1) = []; %DELETES INCOHERENT VALUES 
    else
        break
    end
end
E = mean(E); %AVERAGES ELASTIC MODULUS VECTOR WHICH EACH ELEMENT IS THE SAME AS OTHERS
[TS,I] = max(stress); %FINDS THE TENSILE STRENGTH

fp = TS/max(elastic_region);    %FINDS A FACTOR TO PROYECTATE THE ELASTIC REGION FOR A BETTER COMPREHENSION
proy_x = strain(d,1) * fp;      %MULTIPLIES STRAIN AT THE END OF ELASTIC REGION BETWEEN THE FACTOR

ggg = gradient(stress,strain);
for i = length(ggg):-1:1
    if ggg(i) > 0
        break
    end
end

graph = figure(1);
plot(strain, stress,'k')
hold on
plot([0,proy_x],[0,TS],'--b')  
scatter(strain(I),TS,'filled','b')
stem(strain(i),stress(i),'-.pr','filled','LineWidth',1)

set(gca,'FontName','Times New Roman')     %ADJUST FONT SIZE OF AXIS TO 14
xlabel('Strain (\epsilon: %)','FontName','Times New Roman'), %PERSONALIZATION AXIS X
ylabel('Stress (\sigma: Pa)','FontName','Times New Roman')   %PERSONALIZATION AXIS Y

legend('Stress-Strain Curve',strcat('Elastic Modulus: ',num2str(E,'%.3E')),strcat('Tensile Strength: ',num2str(TS,'%.3E')),strcat('Elongation at Break: ',num2str(strain(i),'%.3f')),'Location','best');

saveas(graph,strcat(testID,'.svg'))