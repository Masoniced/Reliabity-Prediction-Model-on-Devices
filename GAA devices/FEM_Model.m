function out = model
%
% matlab_comsol.m
%
% Model exported on Sep 6 2015, 14:25 by COMSOL 4.4.0.248.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('GAA');

model.modelPath('C:\Users\Sen\Documents\Comsol Multiphysics Projects');



%% 3D

model.modelNode.create('comp2');

model.geom.create('geom2', 3);

model.mesh.create('mesh2', 'geom2');

model.geom('geom2').lengthUnit('nm');

model.param.set('HKT', '2');HKT=2;
model.param.set('ILT', '2');ILT=2;
model.param.set('MGT', '5');MGT=5;
model.param.set('SBT', '10');SBT=10;
model.param.set('PST', '3');PST=3;
model.param.set('GAL', '60');PST=60;



model.geom('geom2').create('cyl1', 'Cylinder');
model.geom('geom2').feature('cyl1').set('r', 'HKT+ILT+MGT+SBT/2+PST');
model.geom('geom2').feature('cyl1').set('h', 'GAL');
model.geom('geom2').feature('cyl1').setIndex('layer', 'PST', 0);
model.geom('geom2').feature('cyl1').setIndex('layer', 'MGT', 1);
model.geom('geom2').feature('cyl1').setIndex('layer', 'HKT', 2);
model.geom('geom2').feature('cyl1').setIndex('layer', 'ILT', 3);
model.geom('geom2').feature('cyl1').set('type', 'solid');
model.geom('geom2').feature('cyl1').set('pos', {'0' '-GAL/2' '0'});
model.geom('geom2').feature('cyl1').set('axistype', 'y');
model.geom('geom2').feature('cyl1').set('createselection', 'on');


model.geom('geom2').geomRep('comsol');
model.geom('geom2').runAll;
model.geom('geom2').run;



%% Mat Setting

model.material.create('mat1', 'Common', 'comp2');
model.material('mat1').name('HfO2');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'25' '0' '0' '0' '25' '0' '0' '0' '25'});
model.material('mat1').propertyGroup('def').set('electricconductivity', {'1e-15' '0' '0' '0' '1e-15' '0' '0' '0' '1e-15'});
model.material('mat1').propertyGroup('def').set('thermalconductivity', {'0.8' '0' '0' '0' '0.8' '0' '0' '0' '0.8'});
model.material('mat1').propertyGroup('def').set('density', '9.68e3');
model.material('mat1').propertyGroup('def').set('heatcapacity', '120');
model.material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'6e-6' '0' '0' '0' '6e-6' '0' '0' '0' '6e-6'});
model.material('mat1').set('family', 'plastic');

model.material.create('mat2', 'Common', 'comp2');
model.material('mat2').name('SiO2');
model.material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'12.3e-6' '0' '0' '0' '12.3e-6' '0' '0' '0' '12.3e-6'});
model.material('mat2').propertyGroup('def').set('density', '2.648e3');
model.material('mat2').propertyGroup('def').set('electricconductivity', {'1e-14' '0' '0' '0' '1e-14' '0' '0' '0' '1e-14'});
model.material('mat2').propertyGroup('def').set('heatcapacity', '700');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'6' '0' '0' '0' '6' '0' '0' '0' '6'});
model.material('mat2').propertyGroup('def').set('resistivity', {'1e17' '0' '0' '0' '1e17' '0' '0' '0' '1e17'});
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat2').propertyGroup('def').set('poissonsratio', '0.17');
model.material('mat2').set('family', 'plastic');

model.material.create('mat3', 'Common', 'comp2');
model.material('mat3').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat3').name('Highly Doped Silicon');
model.material('mat3').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').propertyGroup('def').set('electricconductivity', {'7.69[S/m]' '0' '0' '0' '7.69[S/m]' '0' '0' '0' '7.69[S/m]'});
model.material('mat3').propertyGroup('def').set('thermalexpansioncoefficient', {'2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]'});
model.material('mat3').propertyGroup('def').set('heatcapacity', '700[J/(kg*K)]');
model.material('mat3').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.material('mat3').propertyGroup('def').set('density', '2329[kg/m^3]');
model.material('mat3').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]'});
model.material('mat3').propertyGroup('Enu').set('youngsmodulus', '170e9[Pa]');
model.material('mat3').propertyGroup('Enu').set('poissonsratio', '0.28');
model.material('mat3').propertyGroup('RefractiveIndex').set('n', '');
model.material('mat3').propertyGroup('RefractiveIndex').set('ki', '');
model.material('mat3').propertyGroup('RefractiveIndex').set('n', {'3.48' '0' '0' '0' '3.48' '0' '0' '0' '3.48'});
model.material('mat3').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.material('mat3').set('family', 'plastic');

model.material.create('mat4', 'Common', 'comp2');
model.material('mat4').name('TiN');
model.material('mat4').propertyGroup('def').set('thermalexpansioncoefficient', {'9.3e-6' '0' '0' '0' '9.3e-6' '0' '0' '0' '9.3e-6'});
model.material('mat4').propertyGroup('def').set('density', '5.22e3');
model.material('mat4').propertyGroup('def').set('electricconductivity', {'1.3e6' '0' '0' '0' '1.3e6' '0' '0' '0' '1.3e6'});
model.material('mat4').propertyGroup('def').set('heatcapacity', '750');
model.material('mat4').propertyGroup('def').set('relpermittivity', {'1200' '0' '0' '0' '1200' '0' '0' '0' '1200'});
model.material('mat4').propertyGroup('def').set('resistivity', {'7.69e-7' '0' '0' '0' '7.69e-7' '0' '0' '0' '7.69e-7'});
model.material('mat4').propertyGroup('def').set('thermalconductivity', {'19.2' '0' '0' '0' '19.2' '0' '0' '0' '19.2'});
model.material('mat4').propertyGroup('def').set('poissonsratio', '0.25');
model.material('mat4').set('family', 'plastic');

model.material.create('mat5', 'Common', 'comp2');
model.material('mat5').name('PIT HfO2');
model.material('mat5').propertyGroup('def').set('electricconductivity', {'1.42857e-4' '0' '0' '0' '1.42857e-4' '0' '0' '0' '1.42857e-4'});
model.material('mat5').propertyGroup('def').set('relpermittivity', {'25' '0' '0' '0' '25' '0' '0' '0' '25'});
model.material('mat5').propertyGroup('def').set('thermalconductivity', {'0.8' '0' '0' '0' '0.8' '0' '0' '0' '0.8'});
model.material('mat5').propertyGroup('def').set('density', '8.68e3');
model.material('mat5').propertyGroup('def').set('heatcapacity', '120');
model.material('mat5').set('family', 'plastic');

model.material.create('mat6', 'Common', 'comp2');
model.material('mat6').name('PIT SiO2');
model.material('mat6').propertyGroup('def').set('relpermittivity', {'8.4' '0' '0' '0' '8.4' '0' '0' '0' '8.4'});
model.material('mat6').propertyGroup('def').set('density', '2.648e3');
model.material('mat6').propertyGroup('def').set('heatcapacity', '700');
model.material('mat6').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat6').propertyGroup('def').set('electricconductivity', {'1.2e-3' '0' '0' '0' '1.2e-3' '0' '0' '0' '1.2e-3'});
model.material('mat6').set('family', 'plastic');

model.material.create('mat7', 'Common', 'comp2');
model.material('mat7').name('BD HfO2');
model.material('mat7').propertyGroup('def').set('electricconductivity', {'1.42857' '0' '0' '0' '1.42857' '0' '0' '0' '1.42857'});
model.material('mat7').propertyGroup('def').set('relpermittivity', {'25' '0' '0' '0' '25' '0' '0' '0' '25'});
model.material('mat7').propertyGroup('def').set('thermalconductivity', {'0.8' '0' '0' '0' '0.8' '0' '0' '0' '0.8'});
model.material('mat7').propertyGroup('def').set('density', '8.68e3');
model.material('mat7').propertyGroup('def').set('heatcapacity', '120');
model.material('mat7').set('family', 'plastic');

model.material.create('mat8', 'Common', 'comp2');
model.material('mat8').name('BD SiO2');
model.material('mat8').propertyGroup('def').set('relpermittivity', {'8.4' '0' '0' '0' '8.4' '0' '0' '0' '8.4'});
model.material('mat8').propertyGroup('def').set('density', '2.648e3');
model.material('mat8').propertyGroup('def').set('heatcapacity', '700');
model.material('mat8').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat8').propertyGroup('def').set('electricconductivity', {'1.2e1' '0' '0' '0' '1.2e1' '0' '0' '0' '1.2e1'});
model.material('mat8').set('family', 'plastic');

model.material.create('mat9', 'Common', 'comp2');
model.material('mat9').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat9').name('PS Gate');
model.material('mat9').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat9').propertyGroup('def').set('electricconductivity', {'1e5[S/m]' '0' '0' '0' '1e5[S/m]' '0' '0' '0' '1e5[S/m]'});
model.material('mat9').propertyGroup('def').set('thermalexpansioncoefficient', {'4.4e-6[1/K]' '0' '0' '0' '4.4e-6[1/K]' '0' '0' '0' '4.4e-6[1/K]'});
model.material('mat9').propertyGroup('def').set('heatcapacity', '678[J/(kg*K)]');
model.material('mat9').propertyGroup('def').set('relpermittivity', {'4.5' '0' '0' '0' '4.5' '0' '0' '0' '4.5'});
model.material('mat9').propertyGroup('def').set('density', '2648[kg/m^3]');
model.material('mat9').propertyGroup('def').set('thermalconductivity', {'22.4[W/(m*K)]' '0' '0' '0' '22.4[W/(m*K)]' '0' '0' '0' '22.4[W/(m*K)]'});
model.material('mat9').propertyGroup('Enu').set('youngsmodulus', '160e9[Pa]');
model.material('mat9').propertyGroup('Enu').set('poissonsratio', '0.22');
model.material('mat9').set('family', 'plastic');


model.geom('geom2').geomRep('comsol');
model.geom('geom2').runAll;
model.geom('geom2').run;

% model.material('mat1').selection.set([5 6 12 15]);
% model.material('mat2').selection.set([7 8 13 14]);
% model.material('mat3').selection.set([9]);
% model.material('mat4').selection.set([3 4 11 16]);
% model.material('mat9').selection.set([1 2 10 17]);


%% Defect

S=load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\Results\2.1V_6_GAA.mat');
M=S.M;
R=S.R;
T_hk=S.T_hk;
MM=M(1:100,:,:);
RR=R(1:100,:,:);
[a1,b1,c1]=size(MM);
% l_x=(2*FinH+GateW-2*ILT)/a1;
% l_y=(GateL-MGT)/b1;
% l_z=(HKT+ILT)/c1;

hf_name=1;
il_name=1;
hf_bd_cell={};
il_bd_cell={};
hf_pit_cell={};
il_pit_cell={};
con1=1;
con2=1;
con3=1;
con4=1;

for i=1:a1*b1*c1
   [x,y,z]=ind2sub([a1,b1,c1],i);
   
   if MM(i)==1
%% BD Defect       
       if i<=a1*b1*T_hk
           
           pos = [cos(y*2*pi/b1)*(SBT/2+2+0.4*(4-z+0.5)),(x-50)*0.4,sin(y*2*pi/b1)*(SBT/2+2+0.4*(4-z+0.5))];
           model.geom('geom2').feature.create(num2str(hf_name), 'Block');
           model.geom('geom2').lengthUnit('nm');
           model.geom('geom2').feature(num2str(hf_name)).set('size', [0.4,0.4,0.4]);
           model.geom('geom2').feature(num2str(hf_name)).set('axistype', 'y');
           model.geom('geom2').feature(num2str(hf_name)).set('rot', y*360/b1);
           model.geom('geom2').feature(num2str(hf_name)).set('base', 'center');
           model.geom('geom2').feature(num2str(hf_name)).set('pos', [pos(1),pos(2),pos(3)]);  
           model.geom('geom2').feature(num2str(hf_name)).set('type', 'solid');
           model.geom('geom2').run(num2str(hf_name));          
           
           hf_bd_cell(con1)={num2str(hf_name)}; 
           con1=con1+1;
           hf_name=hf_name+1;
       else

           pos = [cos(y*2*pi/b1)*(SBT/2+0.4*(8-z+0.5)),(x-50)*0.4,sin(y*2*pi/b1)*(SBT/2+0.4*(8-z+0.5))];
           model.geom('geom2').feature.create(num2str(-il_name), 'Block');
           model.geom('geom2').lengthUnit('nm');
           model.geom('geom2').feature(num2str(-il_name)).set('size', [0.4,0.4,0.4]);
           model.geom('geom2').feature(num2str(-il_name)).set('axistype', 'y');
           model.geom('geom2').feature(num2str(-il_name)).set('rot', y*360/b1);
           model.geom('geom2').feature(num2str(-il_name)).set('base', 'center');
           model.geom('geom2').feature(num2str(-il_name)).set('pos', [pos(1),pos(2),pos(3)]);  
           model.geom('geom2').feature(num2str(-il_name)).set('type', 'solid');
           model.geom('geom2').run(num2str(-il_name));          
           
           il_bd_cell(con2)={num2str(-il_name)}; 
           con2=con2+1;
           il_name=il_name+1;
       end
   elseif MM(i)~=1&&RR(i)==1
%% Non BD PIT       
       if i<=a1*b1*T_hk
           
           pos = [cos(y*2*pi/b1)*(SBT/2+2+0.4*(4-z+0.5)),(x-50)*0.4,sin(y*2*pi/b1)*(SBT/2+2+0.4*(4-z+0.5))];
           model.geom('geom2').feature.create(num2str(hf_name), 'Block');
           model.geom('geom2').lengthUnit('nm');
           model.geom('geom2').feature(num2str(hf_name)).set('size', [0.4,0.4,0.4]);
           model.geom('geom2').feature(num2str(hf_name)).set('axistype', 'y');
           model.geom('geom2').feature(num2str(hf_name)).set('rot', y*360/b1);
           model.geom('geom2').feature(num2str(hf_name)).set('base', 'center');
           model.geom('geom2').feature(num2str(hf_name)).set('pos', [pos(1),pos(2),pos(3)]);  
           model.geom('geom2').feature(num2str(hf_name)).set('type', 'solid');
           model.geom('geom2').run(num2str(hf_name));          
           
           hf_pit_cell(con3)={num2str(hf_name)}; 
           con3=con3+1;
           hf_name=hf_name+1;
       else
           
           pos = [cos(y*2*pi/b1)*(SBT/2+0.4*(8-z+0.5)),(x-50)*0.4,sin(y*2*pi/b1)*(SBT/2+0.4*(8-z+0.5))];
           model.geom('geom2').feature.create(num2str(-il_name), 'Block');
           model.geom('geom2').lengthUnit('nm');
           model.geom('geom2').feature(num2str(-il_name)).set('size', [0.4,0.4,0.4]);
           model.geom('geom2').feature(num2str(-il_name)).set('axistype', 'y');
           model.geom('geom2').feature(num2str(-il_name)).set('rot', y*360/b1);
           model.geom('geom2').feature(num2str(-il_name)).set('base', 'center');
           model.geom('geom2').feature(num2str(-il_name)).set('pos', [pos(1),pos(2),pos(3)]);  
           model.geom('geom2').feature(num2str(-il_name)).set('type', 'solid');
           model.geom('geom2').run(num2str(-il_name));          
           
           il_pit_cell(con4)={num2str(-il_name)}; 
           con4=con4+1;
           il_name=il_name+1;
       end
   end
end


%% 

model.geom('geom2').feature.create('HF_BD', 'Union');
model.geom('geom2').feature('HF_BD').selection('input').set(hf_bd_cell);
model.geom('geom2').feature('HF_BD').set('createselection', 'on');
model.geom('geom2').run('HF_BD');

model.geom('geom2').feature.create('IL_BD', 'Union');
model.geom('geom2').feature('IL_BD').selection('input').set(il_bd_cell);
model.geom('geom2').feature('IL_BD').set('createselection', 'on');
model.geom('geom2').run('IL_BD');

model.geom('geom2').feature.create('HK_PIT', 'Union');
model.geom('geom2').feature('HK_PIT').selection('input').set(hf_pit_cell);
model.geom('geom2').feature('HK_PIT').set('createselection', 'on');
model.geom('geom2').run('HK_PIT');

model.geom('geom2').feature.create('IL_PIT', 'Union');
model.geom('geom2').feature('IL_PIT').selection('input').set(il_pit_cell);
model.geom('geom2').feature('IL_PIT').set('createselection', 'on');
model.geom('geom2').run('IL_PIT');

% model.material('mat7').selection.named('geom2_HF_BD_dom');
% model.material('mat8').selection.named('geom2_IL_BD_dom');
% model.material('mat5').selection.named('geom2_HK_PIT_dom');
% model.material('mat6').selection.named('geom2_IL_PIT_dom');
%%

model.geom('geom2').geomRep('comsol');
model.geom('geom2').runAll;




mphgeom(model, 'geom2');




out = model;

model = ModelUtil.model('GAA');

