function out = model
%
% TwoDNVM.m
%
% Model exported on Dec 11 2015, 15:06 by COMSOL 5.1.0.136.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\SUTD-ZJU\Desktop');

model.comments(['Untitled\n\n']);

model.modelNode.create('comp1');

model.geom.create('geom1', 2);

model.mesh.create('mesh1', 'geom1');


L_NC=7;
W_NC=1;
NC_d=2;
NC_i=2;
Edge=10;
OXH=10;
MGH=10;
OXOFFSET=2;


l_x=(L_NC+1)/2;
l_y=(W_NC+1)/2;


model.param.set('NC_d', '2');
model.param.set('OXH', '10');
model.param.set('NC_i', '2');
model.param.set('MGH', '10');
model.param.set('L_NC', '7');
model.param.set('W_NC', '1');
model.param.set('Edge', '5');
model.param.set('OXOFFSET', '2');

NC_number=1;
NC_cell=cell(1,L_NC);


%% NC

for i=1:(L_NC*1)

    NC_name=num2str(NC_number);
    x=i;
    
    model.geom('geom1').lengthUnit('nm');
    model.geom('geom1').create(NC_name, 'Circle');
    model.geom('geom1').feature(NC_name).set('r', 'NC_d/2');
    model.geom('geom1').feature(NC_name).set('pos', [(x-l_x)*(NC_d+NC_i),0]);
    model.geom('geom1').feature(NC_name).set('type', 'solid');
    model.geom('geom1').feature(NC_name).set('base', 'center');
%    model.geom('geom1').feature(NC_name).set('createselection', 'on');
    
%    model.geom('geom1').feature(NC_name).set('createselection', 'on');


    NC_cell(NC_number)={NC_name};
    NC_number=NC_number+1;
    
end

model.geom('geom1').feature.create('NC', 'Union');
model.geom('geom1').feature('NC').selection('input').set(NC_cell);
model.geom('geom1').feature('NC').set('createselection', 'on');
model.geom('geom1').run('NC');


%% Oxide

model.geom('geom1').create('r1', 'Rectangle');
model.geom('geom1').lengthUnit('nm');
model.geom('geom1').feature('r1').set('size', {'(NC_d+NC_i)*L_NC+Edge', 'OXH'});
model.geom('geom1').feature('r1').set('pos', {'0', 'OXOFFSET'});
model.geom('geom1').feature('r1').set('base', 'center');
model.geom('geom1').feature('r1').set('createselection', 'on');


%% Electrode

model.geom('geom1').create('r2', 'Rectangle');
model.geom('geom1').lengthUnit('nm');
model.geom('geom1').feature('r2').set('size', {'(NC_d+NC_i)*L_NC+Edge', 'MGH'});
model.geom('geom1').feature('r2').set('pos', {'0', 'OXOFFSET+MGH/2+OXH/2'});
model.geom('geom1').feature('r2').set('base', 'center');
model.geom('geom1').feature('r2').set('createselection', 'on');

model.geom('geom1').create('r3', 'Rectangle');
model.geom('geom1').lengthUnit('nm');
model.geom('geom1').feature('r3').set('size', {'(NC_d+NC_i)*L_NC+Edge', 'MGH'});
model.geom('geom1').feature('r3').set('pos', {'0', 'OXOFFSET-MGH/2-OXH/2'});
model.geom('geom1').feature('r3').set('base', 'center');
model.geom('geom1').feature('r3').set('createselection', 'on');



%% Defects


S=load('C:\Users\Mason\Documents\Project\Matlab Project\NC\Results\Test_C_20_40_60.mat');
M=S.NC_list{2};
[a,b,c]=size(M(:,:,1:end-2));
u=zeros(a,c);
for i=1:b/2
    u=u+reshape(M(:,i,1:c),[a,c]);
end

u=u';
hf_name=-1;
il_name=1;
hf_name_cell={};
il_name_cell={};
hf_fila_cell={};
il_fila_cell={};
hf_gb_cell={};
hf_gbn_cell={};
con1=1;
con2=1;
con3=1;
con4=1;
con5=1;
con6=1;
for i=1:a*c
   [x,y]=ind2sub([c,a],i);

   if u(x,y)>=1
       cunt=0;
       for ik= 1:L_NC
           distance=sqrt(((y-0.5-a/2)*0.4-(ik-l_x)*(NC_d+NC_i))^2+((-c/2+x-0.5)*0.4+OXOFFSET)^2);
           if distance <= NC_d/2+0.2
               cunt=cunt+1;
           end
       end
       if cunt == 0
           model.geom('geom1').feature.create(num2str(hf_name), 'Rectangle');
           model.geom('geom1').lengthUnit('nm');
           model.geom('geom1').feature(num2str(hf_name)).set('size', [0.4,0.4]);
           model.geom('geom1').feature(num2str(hf_name)).set('base', 'center');
           model.geom('geom1').feature(num2str(hf_name)).set('pos', [(y-0.5-a/2)*0.4,(-c/2+x-0.5)*0.4+OXOFFSET]);
           model.geom('geom1').feature(num2str(hf_name)).set('type', 'solid');
           model.geom('geom1').run(num2str(hf_name));

           hf_name_cell(con1)={num2str(hf_name)};
           con1=con1+1;               

           hf_name=hf_name-1;
       end

       

       
%    elseif R(x,y)==1&&M(x,y)~=1
%        
%        model.geom('geom1').feature.create(num2str(-il_name), 'Rectangle');
%        model.geom('geom1').lengthUnit('nm');
%        model.geom('geom1').feature(num2str(-il_name)).set('size', [0.8,0.8]);
%        model.geom('geom1').feature(num2str(-il_name)).set('base', 'center');
%        model.geom('geom1').feature(num2str(-il_name)).set('pos', [(EGL+2*TNOX+2*FGL/3-2*CGSNOL)/2+(y-0.5)*0.8,FGOX/2+(x-0.5-a1/2)*0.8]);
%        model.geom('geom1').feature(num2str(-il_name)).set('type', 'solid');
%        model.geom('geom1').run(num2str(-il_name));
%        
%        il_trap_cell(con6)={num2str(hf_name)};
%        con6=con6+1;
%        
%        il_name=il_name+1;
              
   end
end

model.geom('geom1').feature.create('hf', 'Union');
model.geom('geom1').feature('hf').selection('input').set(hf_name_cell);
model.geom('geom1').run('hf');

% model.geom('geom1').feature.create('ilfila', 'Union');
% model.geom('geom1').feature('ilfila').selection('input').set(il_fila_cell);
% model.geom('geom1').run('ilfila');

% model.geom('geom1').feature.create('iltrap', 'Union');
% model.geom('geom1').feature('iltrap').selection('input').set(il_trap_cell);
% model.geom('geom1').run('iltrap');

model.geom('geom1').feature('hf').set('createselection', 'on');
% model.geom('geom1').feature('ilfila').set('createselection', 'on');
% model.geom('geom1').feature('iltrap').set('createselection', 'on');




%% Materials

% NC
model.material.create('mat1', 'Common', 'comp1');
model.material('mat1').label('NC');
model.material('mat1').set('family', 'gold');
model.material('mat1').propertyGroup('def').set('dL', '(dL(T[1/K])-dL(Tempref[1/K]))/(1+dL(Tempref[1/K]))');
model.material('mat1').propertyGroup('def').set('CTE', 'CTE(T[1/K])[1/K]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', 'k_solid_1(T[1/K])[W/(m*K)]');
model.material('mat1').propertyGroup('def').set('resistivity', 'res_solid_1(T[1/K])[ohm*m]');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'1e4' '0' '0' '0' '1e4' '0' '0' '0' '1e4'});
model.material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '(alpha(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha(T[1/K])[1/K]-alpha(Tempref[1/K])[1/K])/(T-Tempref),d(alpha(T[1/K]),T)[1/K]))/(1+alpha(Tempref[1/K])[1/K]*(Tempref-293[K]))');
model.material('mat1').propertyGroup('def').set('heatcapacity', 'C_solid_1(T[1/K])[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('HC', 'HC_solid_1(T[1/K])[J/(mol*K)]');
model.material('mat1').propertyGroup('def').set('electricconductivity', 'sigma_solid_1(T[1/K])[S/m]');
model.material('mat1').propertyGroup('def').set('density', 'rho(T[1/K])[kg/m^3]');
model.material('mat1').propertyGroup('def').set('VP', 'VP_solid_1(T[1/K])[Pa]');
model.material('mat1').propertyGroup('def').func.create('dL', 'Piecewise');
model.material('mat1').propertyGroup('def').func('dL').set('funcname', 'dL');
model.material('mat1').propertyGroup('def').func('dL').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('dL').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('dL').set('pieces', {'0.0' '92.0' '-0.003245794-2.191545E-7*T^1+9.931873E-8*T^2-3.27165E-10*T^3+4.101526E-13*T^4'; '92.0' '1003.0' '-0.003767805+1.179943E-5*T^1+3.792275E-9*T^2-5.959619E-13*T^3'});
model.material('mat1').propertyGroup('def').func.create('CTE', 'Piecewise');
model.material('mat1').propertyGroup('def').func('CTE').set('funcname', 'CTE');
model.material('mat1').propertyGroup('def').func('CTE').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('CTE').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('CTE').set('pieces', {'86.0' '200.0' '4.040319E-6+1.106436E-7*T^1-4.389255E-10*T^2+5.990814E-13*T^3'; '200.0' '1003.0' '1.237291E-5+5.169284E-9*T^1-5.528527E-14*T^2'});
model.material('mat1').propertyGroup('def').func.create('k_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('k_solid_1').set('funcname', 'k_solid_1');
model.material('mat1').propertyGroup('def').func('k_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('k_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('k_solid_1').set('pieces', {'0.0' '14.0' '496.8118*T^1+41.00857*T^2-10.46095*T^3+0.5579006*T^4-0.009429161*T^5'; '14.0' '45.0' '7977.526-554.7566*T^1+14.56751*T^2-0.1468484*T^3+3.231278E-4*T^4'; '45.0' '85.0' '3568.586-199.4435*T^1+5.171042*T^2-0.06880706*T^3+4.612091E-4*T^4-1.23215E-6*T^5'; '85.0' '1338.0' '330.6431-0.02536626*T^1-8.191375E-5*T^2+6.792908E-8*T^3-2.15362E-11*T^4'});
model.material('mat1').propertyGroup('def').func.create('res_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('res_solid_1').set('funcname', 'res_solid_1');
model.material('mat1').propertyGroup('def').func('res_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('res_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('res_solid_1').set('pieces', {'1.0' '25.0' '2.194815E-10+1.027514E-12*T^1-3.409278E-13*T^2+3.010164E-14*T^3'; '25.0' '60.0' '1.390641E-9-1.228888E-10*T^1+4.182944E-12*T^2-2.773737E-14*T^3'; '60.0' '400.0' '-2.210068E-9+9.057611E-11*T^1-4.632985E-14*T^2+6.950205E-17*T^3'; '400.0' '1338.0' '-1.145028E-9+7.877041E-11*T^1-4.720065E-16*T^2+1.275961E-17*T^3'});
model.material('mat1').propertyGroup('def').func.create('alpha', 'Piecewise');
model.material('mat1').propertyGroup('def').func('alpha').set('funcname', 'alpha');
model.material('mat1').propertyGroup('def').func('alpha').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('alpha').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('alpha').set('pieces', {'0.0' '74.0' '1.109296E-5+4.023822E-8*T^1-2.476496E-10*T^2+7.218616E-13*T^3-7.164369E-16*T^4'; '74.0' '133.0' '1.221837E-5+1.15298E-8*T^1-1.573678E-11*T^2'; '133.0' '1003.0' '1.310696E-5+2.762787E-9*T^1'});
model.material('mat1').propertyGroup('def').func.create('C_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('C_solid_1').set('funcname', 'C_solid_1');
model.material('mat1').propertyGroup('def').func('C_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('C_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('C_solid_1').set('pieces', {'293.0' '1338.0' '399352.2*T^-2+114.8987+0.03228805*T^1'});
model.material('mat1').propertyGroup('def').func.create('HC_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('HC_solid_1').set('funcname', 'HC_solid_1');
model.material('mat1').propertyGroup('def').func('HC_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('HC_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('HC_solid_1').set('pieces', {'293.0' '1338.0' '78659.2*T^-2+22.63126+0.00635968*T^1'});
model.material('mat1').propertyGroup('def').func.create('sigma_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('sigma_solid_1').set('funcname', 'sigma_solid_1');
model.material('mat1').propertyGroup('def').func('sigma_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('sigma_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('sigma_solid_1').set('pieces', {'1.0' '25.0' '1/(3.010164E-14*T^3-3.409278E-13*T^2+1.027514E-12*T+2.194815E-10)'; '25.0' '60.0' '1/(-2.773737E-14*T^3+4.182944E-12*T^2-1.228888E-10*T+1.390641E-09)'; '60.0' '400.0' '1/(6.950205E-17*T^3-4.632985E-14*T^2+9.057611E-11*T-2.210068E-09)'; '400.0' '1338.0' '1/(1.275961E-17*T^3-4.720065E-16*T^2+7.877041E-11*T-1.145028E-09)'});
model.material('mat1').propertyGroup('def').func.create('rho', 'Piecewise');
model.material('mat1').propertyGroup('def').func('rho').set('funcname', 'rho');
model.material('mat1').propertyGroup('def').func('rho').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('rho').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('rho').set('pieces', {'0.0' '86.0' '19471.18+0.1076816*T^1-0.00817626*T^2+3.830366E-5*T^3-7.266799E-8*T^4'; '86.0' '1003.0' '19501.44-0.6933844*T^1-2.041944E-4*T^2+4.297982E-8*T^3'});
model.material('mat1').propertyGroup('def').func.create('VP_solid_1', 'Piecewise');
model.material('mat1').propertyGroup('def').func('VP_solid_1').set('funcname', 'VP_solid_1');
model.material('mat1').propertyGroup('def').func('VP_solid_1').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('VP_solid_1').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('VP_solid_1').set('pieces', {'293.0' '1337.0' '(exp((-1.934300e+004/T-7.479000e-001*log10(T)+1.203281e+001)*log(10.0)))*1.333200e+002'});
model.material('mat1').propertyGroup('def').addInput('temperature');
model.material('mat1').propertyGroup('def').addInput('strainreferencetemperature');
model.material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat1').propertyGroup('Enu').set('youngsmodulus', 'E(T[1/K])[Pa]');
model.material('mat1').propertyGroup('Enu').set('poissonsratio', 'nu(T[1/K])');
model.material('mat1').propertyGroup('Enu').func.create('E', 'Piecewise');
model.material('mat1').propertyGroup('Enu').func('E').set('funcname', 'E');
model.material('mat1').propertyGroup('Enu').func('E').set('arg', 'T');
model.material('mat1').propertyGroup('Enu').func('E').set('extrap', 'constant');
model.material('mat1').propertyGroup('Enu').func('E').set('pieces', {'0.0' '93.0' '8.190431E10-2.049635E7*T^1-8937.063*T^2'; '93.0' '293.0' '8.249999E10-3.153455E7*T^1+48127.92*T^2-72.49967*T^3+0.02725845*T^4'; '293.0' '1280.0' '7.736733E10+425187.2*T^1-20070.07*T^2'});
model.material('mat1').propertyGroup('Enu').func.create('nu', 'Piecewise');
model.material('mat1').propertyGroup('Enu').func('nu').set('funcname', 'nu');
model.material('mat1').propertyGroup('Enu').func('nu').set('arg', 'T');
model.material('mat1').propertyGroup('Enu').func('nu').set('extrap', 'constant');
model.material('mat1').propertyGroup('Enu').func('nu').set('pieces', {'293.0' '1280.0' '0.4403494+8.884688E-6*T^1-1.158876E-8*T^2+1.023495E-11*T^3'});
model.material('mat1').propertyGroup('Enu').addInput('temperature');
model.material('mat1').propertyGroup.create('KG', 'Bulk modulus and shear modulus');
model.material('mat1').propertyGroup('KG').set('G', 'mu(T[1/K])[Pa]');
model.material('mat1').propertyGroup('KG').set('K', 'kappa(T[1/K])[Pa]');
model.material('mat1').propertyGroup('KG').func.create('mu', 'Piecewise');
model.material('mat1').propertyGroup('KG').func('mu').set('funcname', 'mu');
model.material('mat1').propertyGroup('KG').func('mu').set('arg', 'T');
model.material('mat1').propertyGroup('KG').func('mu').set('extrap', 'constant');
model.material('mat1').propertyGroup('KG').func('mu').set('pieces', {'293.0' '1280.0' '2.682905E10+142041.4*T^1-7037.486*T^2'});
model.material('mat1').propertyGroup('KG').func.create('kappa', 'Piecewise');
model.material('mat1').propertyGroup('KG').func('kappa').set('funcname', 'kappa');
model.material('mat1').propertyGroup('KG').func('kappa').set('arg', 'T');
model.material('mat1').propertyGroup('KG').func('kappa').set('extrap', 'constant');
model.material('mat1').propertyGroup('KG').func('kappa').set('pieces', {'293.0' '1280.0' '2.224306E11-913189.7*T^1-37420.01*T^2'});
model.material('mat1').propertyGroup('KG').addInput('temperature');
model.material('mat1').set('family', 'gold');

% OXIDE
model.material.create('mat2');
model.material('mat2').name('HfO2');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'25' '0' '0' '0' '25' '0' '0' '0' '25'});
model.material('mat2').propertyGroup('def').set('electricconductivity', {'1e-15' '0' '0' '0' '1e-15' '0' '0' '0' '1e-15'});
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'0.8' '0' '0' '0' '0.8' '0' '0' '0' '0.8'});
model.material('mat2').propertyGroup('def').set('density', '9.68e3');
model.material('mat2').propertyGroup('def').set('heatcapacity', '120');
model.material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'6e-6' '0' '0' '0' '6e-6' '0' '0' '0' '6e-6'});
model.material('mat2').set('family', 'plastic');

% Substrate
model.material.create('mat3');
model.material('mat3').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat3').name('Highly Doped Silicon');
model.material('mat3').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').propertyGroup('def').set('electricconductivity', {'2.5e4[S/m]' '0' '0' '0' '2.5e4[S/m]' '0' '0' '0' '2.5e4[S/m]'});
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

% BD HK
model.material.create('mat4');
model.material('mat4').name('BD HfO2');
model.material('mat4').propertyGroup('def').set('electricconductivity', {'1.42857' '0' '0' '0' '1.42857' '0' '0' '0' '1.42857'});
model.material('mat4').propertyGroup('def').set('relpermittivity', {'25' '0' '0' '0' '25' '0' '0' '0' '25'});
model.material('mat4').propertyGroup('def').set('thermalconductivity', {'0.8' '0' '0' '0' '0.8' '0' '0' '0' '0.8'});
model.material('mat4').propertyGroup('def').set('density', '8.68e3');
model.material('mat4').propertyGroup('def').set('heatcapacity', '120');


%% 

% model.physics.create('es', 'Electrostatics', 'geom1');
% model.physics('es').feature.create('gnd1', 'Ground', 2);
% model.physics('es').feature.create('pot1', 'ElectricPotential', 2);
% model.physics('es').feature.create('sfcd1', 'SurfaceChargeDensity', 2);
% model.physics('es').feature('pot1').set('V0', '6');
% model.physics('es').feature('sfcd1').selection.named('geom1_NC_bnd');

model.name('combine with matlab.mph');

mphgeom(model);

model.geom('geom1').run;


out = model;

model1 = ModelUtil.model('Model');
