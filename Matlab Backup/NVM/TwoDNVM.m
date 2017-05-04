function out = model
%
% TwoDNVM.m
%
% Model exported on Dec 11 2015, 15:06 by COMSOL 5.1.0.136.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\Meisen\Desktop');

model.comments(['Untitled\n\n']);

model.modelNode.create('comp1');

model.geom.create('geom1', 2);

model.mesh.create('mesh1', 'geom1');


model.param.set('EGL', '178');
model.param.set('EGH', '60');
model.param.set('EGOXH', '45');
model.param.set('MOL', '740');
model.param.set('MOH', '360');
model.param.set('OXH', '200');
model.param.set('CAPH', '18');

model.param.set('TNOX', '11');
model.param.set('FGL', '108');
model.param.set('FGH', '33');
model.param.set('FGOX', '9.6');
model.param.set('CGSNOL', '10');
model.param.set('CGSNOH', '5');
model.param.set('CGH', '55');

model.param.set('WLOXL', '19');
model.param.set('WLOXH', '5');
model.param.set('WLL', '86');
model.param.set('WLH', '70');

model.param.set('CAPSNOH', '20');

model.param.set('DH', '45');
model.param.set('WLCH', '60');
model.param.set('SH', '30');




model.geom('geom1').create('r1', 'Rectangle');
model.geom('geom1').feature('r1').set('size', {'MOL' 'MOH'});
model.geom('geom1').feature('r1').set('pos', {'0' 'MOH/2-100'});
model.geom('geom1').feature('r1').set('base', 'center');
model.geom('geom1').feature('r1').setIndex('layer', '30', 0);
model.geom('geom1').feature('r1').setIndex('layer', '70', 1);
model.geom('geom1').feature('r1').set('createselection', 'on');


model.geom('geom1').create('e1', 'Ellipse');
model.geom('geom1').feature('e1').set('semiaxes', {'EGL/2' 'EGOXH/2'});
model.geom('geom1').feature('e1').set('pos', {'0' 'EGOXH/6'});
model.geom('geom1').feature('e1').set('createselection', 'on');


model.geom('geom1').create('r2', 'Rectangle');
model.geom('geom1').feature('r2').set('size', {'EGL' 'EGH'});
model.geom('geom1').feature('r2').set('pos', {'0' 'EGOXH/3+EGH/2'});
model.geom('geom1').feature('r2').set('base', 'center');
model.geom('geom1').feature('r2').set('createselection', 'on');


model.geom('geom1').create('dif1', 'Difference');
model.geom('geom1').feature('dif1').selection('input').set({'r2'});
model.geom('geom1').feature('dif1').selection('input2').set({'e1'});
model.geom('geom1').feature('dif1').set('createselection', 'on');


model.geom('geom1').feature.duplicate('e2', 'e1');


model.geom('geom1').create('r3', 'Rectangle');
model.geom('geom1').feature('r3').set('size', {'MOL' 'OXH'});
model.geom('geom1').feature('r3').set('pos', {'0' 'OXH/2'});
model.geom('geom1').feature('r3').set('base', 'center');
model.geom('geom1').feature('r3').set('createselection', 'on');


model.geom('geom1').create('uni1', 'Union');
model.geom('geom1').feature('uni1').selection('input').set({'e2' 'r3'});
model.geom('geom1').feature('uni1').set('intbnd', 'off');
model.geom('geom1').feature('uni1').set('createselection', 'on');


model.geom('geom1').create('r4', 'Rectangle');
model.geom('geom1').feature('r4').set('size', {'EGL' 'CAPH'});
model.geom('geom1').feature('r4').set('pos', {'0' 'EGOXH/3+EGH+CAPH/2'});
model.geom('geom1').feature('r4').set('base', 'center');
model.geom('geom1').feature('r4').set('createselection', 'on');


model.geom('geom1').create('r5', 'Rectangle');
model.geom('geom1').feature('r5').set('size', {'FGL' '2*FGH/3'});
model.geom('geom1').feature('r5').set('pos', {'EGL/2+TNOX+FGL/2' 'FGH/3+FGOX'});
model.geom('geom1').feature('r5').set('base', 'center');
model.geom('geom1').feature('r5').set('createselection', 'on');


model.geom('geom1').create('r6', 'Rectangle');
model.geom('geom1').feature('r6').set('size', {'2*FGL/3' 'FGH/3'});
model.geom('geom1').feature('r6').set('pos', {'EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2' '5*FGH/6+FGOX'});
model.geom('geom1').feature('r6').set('base', 'center');
model.geom('geom1').feature('r6').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r7', 'r5');
model.geom('geom1').feature('r7').set('pos', {'-(EGL/2+TNOX+FGL/2)' 'FGH/3+FGOX'});
model.geom('geom1').feature('r7').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r8', 'r6');
model.geom('geom1').feature('r8').set('pos', {'-(EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2)' '5*FGH/6+FGOX'});
model.geom('geom1').feature('r8').set('createselection', 'on');


model.geom('geom1').create('uni2', 'Union');
model.geom('geom1').feature('uni2').selection('input').set({'r5' 'r6'});
model.geom('geom1').feature('uni2').set('intbnd', 'off');
model.geom('geom1').feature('uni2').set('createselection', 'on');


model.geom('geom1').create('uni3', 'Union');
model.geom('geom1').feature('uni3').selection('input').set({'r7' 'r8'});
model.geom('geom1').feature('uni3').set('intbnd', 'off');
model.geom('geom1').feature('uni3').set('createselection', 'on');


model.geom('geom1').create('r9', 'Rectangle');
model.geom('geom1').feature('r9').set('size', {'2*FGL/3' 'CGSNOH'});
model.geom('geom1').feature('r9').set('pos', {'EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2' 'FGH+FGOX+3*CGSNOH/2'});
model.geom('geom1').feature('r9').set('base', 'center');
model.geom('geom1').feature('r9').set('createselection', 'on');


model.geom('geom1').create('r10', 'Rectangle');
model.geom('geom1').feature('r10').set('size', {'2*FGL/3' 'CGH'});
model.geom('geom1').feature('r10').set('pos', {'EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2' 'FGH+FGOX+3*CGSNOH+CGH/2'});
model.geom('geom1').feature('r10').set('base', 'center');
model.geom('geom1').feature('r10').set('createselection', 'on');


model.geom('geom1').create('r11', 'Rectangle');
model.geom('geom1').feature('r11').set('size', {'CGSNOL' 'CGH+3*CGSNOH'});
model.geom('geom1').feature('r11').set('pos', {'EGL/2+TNOX+FGL-CGSNOL/2' 'FGH+FGOX+(CGH+3*CGSNOH)/2'});
model.geom('geom1').feature('r11').set('base', 'center');
model.geom('geom1').feature('r11').set('createselection', 'on');


model.geom('geom1').create('r12', 'Rectangle');
model.geom('geom1').feature('r12').set('size', {'CGSNOL' 'CGH+3*CGSNOH'});
model.geom('geom1').feature('r12').set('pos', {'EGL/2+TNOX+FGL-2*FGL/3-2.5*CGSNOL' 'FGH+FGOX+(CGH+3*CGSNOH)/2'});
model.geom('geom1').feature('r12').set('base', 'center');
model.geom('geom1').feature('r12').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r13', 'r9');
model.geom('geom1').feature('r13').set('pos', {'-(EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2)' 'FGH+FGOX+3*CGSNOH/2'});
model.geom('geom1').feature('r13').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r14', 'r10');
model.geom('geom1').feature('r14').set('pos', {'-(EGL/2+TNOX+FGL-FGL/3-3*CGSNOL/2)' 'FGH+FGOX+3*CGSNOH+CGH/2'});
model.geom('geom1').feature('r14').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r15', 'r11');
model.geom('geom1').feature('r15').set('pos', {'-(EGL/2+TNOX+FGL-CGSNOL/2)' 'FGH+FGOX+(CGH+3*CGSNOH)/2'});
model.geom('geom1').feature('r15').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r16', 'r12');
model.geom('geom1').feature('r16').set('pos', {'-(EGL/2+TNOX+FGL-2*FGL/3-2.5*CGSNOL)' 'FGH+FGOX+(CGH+3*CGSNOH)/2'});
model.geom('geom1').feature('r16').set('createselection', 'on');


model.geom('geom1').create('r17', 'Rectangle');
model.geom('geom1').feature('r17').set('size', {'WLL' 'WLH'});
model.geom('geom1').feature('r17').set('pos', {'EGL/2+TNOX+FGL+WLOXL+WLL/2' 'WLOXH+WLH/2'});
model.geom('geom1').feature('r17').set('base', 'center');
model.geom('geom1').feature('r17').set('createselection', 'on');


model.geom('geom1').create('r18', 'Rectangle');
model.geom('geom1').feature('r18').set('size', {'WLL' 'CAPH'});
model.geom('geom1').feature('r18').set('pos', {'EGL/2+TNOX+FGL+WLOXL+WLL/2' 'WLOXH+WLH+CAPH/2'});
model.geom('geom1').feature('r18').set('base', 'center');
model.geom('geom1').feature('r18').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r19', 'r17');
model.geom('geom1').feature('r19').set('pos', {'-(EGL/2+TNOX+FGL+WLOXL+WLL/2)' 'WLOXH+WLH/2'});
model.geom('geom1').feature('r19').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r20', 'r18');
model.geom('geom1').feature('r20').set('pos', {'-(EGL/2+TNOX+FGL+WLOXL+WLL/2)' 'WLOXH+WLH+CAPH/2'});
model.geom('geom1').feature('r20').set('createselection', 'on');


model.geom('geom1').feature.create('unisel1', 'UnionSelection');
model.geom('geom1').feature('unisel1').name('Gate');
model.geom('geom1').feature('unisel1').set('input', {'dif1' 'uni1' 'uni2' 'r10' 'r14' 'r17' 'r19'});
model.geom('geom1').run('unisel1');


model.geom('geom1').create('r21', 'Rectangle');
model.geom('geom1').feature('r21').set('size', {'CGSNOL' 'FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH)'});
model.geom('geom1').feature('r21').set('pos', {'EGL/2-CGSNOL/2-3' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r21').set('base', 'center');
model.geom('geom1').feature('r21').set('createselection', 'on');


model.geom('geom1').create('r22', 'Rectangle');
model.geom('geom1').feature('r22').set('size', {'CGSNOL' '(FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r22').set('pos', {'EGL/2-3*CGSNOL/2-3' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH))/4'});
model.geom('geom1').feature('r22').set('base', 'center');
model.geom('geom1').feature('r22').set('createselection', 'on');


model.geom('geom1').create('uni4', 'Union');
model.geom('geom1').feature('uni4').selection('input').set({'r21' 'r22'});
model.geom('geom1').feature('uni4').set('intbnd', 'off');
model.geom('geom1').feature('uni4').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r23', 'r21');
model.geom('geom1').feature('r23').set('pos', {'-(EGL/2-CGSNOL/2-3)' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r23').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r24', 'r22');
model.geom('geom1').feature('r24').set('pos', {'-(EGL/2-3*CGSNOL/2-3)' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(EGOXH/3+EGH+CAPH+CGSNOH))/4'});
model.geom('geom1').feature('r24').set('createselection', 'on');


model.geom('geom1').create('uni5', 'Union');
model.geom('geom1').feature('uni5').selection('input').set({'r23' 'r24'});
model.geom('geom1').feature('uni5').set('intbnd', 'off');
model.geom('geom1').feature('uni5').set('createselection', 'on');


model.geom('geom1').create('r25', 'Rectangle');
model.geom('geom1').feature('r25').set('size', {'CGSNOL' 'FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH)'});
model.geom('geom1').feature('r25').set('pos', {'EGL/2+TNOX+FGL+WLOXL+CGSNOL/2' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r25').set('base', 'center');
model.geom('geom1').feature('r25').set('createselection', 'on');


model.geom('geom1').create('r26', 'Rectangle');
model.geom('geom1').feature('r26').set('size', {'CGSNOL' '(FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r26').set('pos', {'EGL/2+TNOX+FGL+WLOXL+3*CGSNOL/2' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH))/4'});
model.geom('geom1').feature('r26').set('base', 'center');
model.geom('geom1').feature('r26').set('createselection', 'on');


model.geom('geom1').create('uni6', 'Union');
model.geom('geom1').feature('uni6').selection('input').set({'r25' 'r26'});
model.geom('geom1').feature('uni6').set('intbnd', 'off');
model.geom('geom1').feature('uni6').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r27', 'r25');
model.geom('geom1').feature('r27').set('pos', {'-(EGL/2+TNOX+FGL+WLOXL+CGSNOL/2)' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r27').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r28', 'r26');
model.geom('geom1').feature('r28').set('pos', {'-(EGL/2+TNOX+FGL+WLOXL+3*CGSNOL/2)' 'EGOXH/3+EGH+CAPH+CGSNOH+(FGH+FGOX+3*CGSNOH+CGH-(WLOXH+WLH+CAPH+CGSNOH))/4'});
model.geom('geom1').feature('r28').set('createselection', 'on');


model.geom('geom1').create('uni7', 'Union');
model.geom('geom1').feature('uni7').selection('input').set({'r27' 'r28'});
model.geom('geom1').feature('uni7').set('intbnd', 'off');
model.geom('geom1').feature('uni7').set('createselection', 'on');


model.geom('geom1').create('r29', 'Rectangle');
model.geom('geom1').feature('r29').set('size', {'FGL+TNOX+WLOXL+6*CGSNOL' 'CAPSNOH'});
model.geom('geom1').feature('r29').set('pos', {'EGL/2+TNOX/2+FGL/2+WLOXL/2' 'FGH+FGOX+3*CGSNOH+CGH+CAPSNOH/2'});
model.geom('geom1').feature('r29').set('base', 'center');
model.geom('geom1').feature('r29').set('createselection', 'on');


model.geom('geom1').create('r30', 'Rectangle');
model.geom('geom1').feature('r30').set('size', {'CAPSNOH' 'FGH+FGOX+3*CGSNOH+CGH+CAPSNOH-(EGOXH/3+EGH+CAPH+CGSNOH)'});
model.geom('geom1').feature('r30').set('pos', {'EGL/2-3*CGSNOL-CAPSNOH/2' '(FGH+FGOX+3*CGSNOH+CGH+CAPSNOH+(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r30').set('base', 'center');
model.geom('geom1').feature('r30').set('createselection', 'on');


model.geom('geom1').create('r31', 'Rectangle');
model.geom('geom1').feature('r31').set('size', {'CAPSNOH' 'FGH+FGOX+3*CGSNOH+CGH+CAPSNOH-(EGOXH/3+EGH+CAPH+CGSNOH)'});
model.geom('geom1').feature('r31').set('pos', {'EGL/2+FGL+TNOX+WLOXL+3*CGSNOL+CAPSNOH/2' '(FGH+FGOX+3*CGSNOH+CGH+CAPSNOH+(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r31').set('base', 'center');
model.geom('geom1').feature('r31').set('createselection', 'on');


model.geom('geom1').create('r32', 'Rectangle');
model.geom('geom1').feature('r32').set('size', {'WLL-3*CGSNOL-CAPSNOH' 'CAPSNOH'});
model.geom('geom1').feature('r32').set('pos', {'EGL/2+FGL+TNOX+WLOXL+(WLL+3*CGSNOL+CAPSNOH)/2' 'WLH+WLOXH+CGSNOH+CAPH+CAPSNOH/2'});
model.geom('geom1').feature('r32').set('base', 'center');
model.geom('geom1').feature('r32').set('createselection', 'on');


model.geom('geom1').create('r33', 'Rectangle');
model.geom('geom1').feature('r33').set('size', {'CAPSNOH' 'WLH+WLOXH+CGSNOH+CAPH+CAPSNOH'});
model.geom('geom1').feature('r33').set('pos', {'EGL/2+FGL+TNOX+WLOXL+WLL+CAPSNOH/2' '(WLH+WLOXH+CGSNOH+CAPH+CAPSNOH)/2'});
model.geom('geom1').feature('r33').set('base', 'center');
model.geom('geom1').feature('r33').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r34', 'r29');
model.geom('geom1').feature('r34').set('pos', {'-(EGL/2+TNOX/2+FGL/2+WLOXL/2)' 'FGH+FGOX+3*CGSNOH+CGH+CAPSNOH/2'});
model.geom('geom1').feature('r34').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r35', 'r30');
model.geom('geom1').feature('r35').set('pos', {'-(EGL/2-3*CGSNOL-CAPSNOH/2)' '(FGH+FGOX+3*CGSNOH+CGH+CAPSNOH+(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r35').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r36', 'r31');
model.geom('geom1').feature('r36').set('pos', {'-(EGL/2+FGL+TNOX+WLOXL+3*CGSNOL+CAPSNOH/2)' '(FGH+FGOX+3*CGSNOH+CGH+CAPSNOH+(EGOXH/3+EGH+CAPH+CGSNOH))/2'});
model.geom('geom1').feature('r36').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r37', 'r32');
model.geom('geom1').feature('r37').set('pos', {'-(EGL/2+FGL+TNOX+WLOXL+(WLL+3*CGSNOL+CAPSNOH)/2)' 'WLH+WLOXH+CGSNOH+CAPH+CAPSNOH/2'});
model.geom('geom1').feature('r37').set('createselection', 'on');


model.geom('geom1').feature.duplicate('r38', 'r33');
model.geom('geom1').feature('r38').set('pos', {'-(EGL/2+FGL+TNOX+WLOXL+WLL+CAPSNOH/2)' '(WLH+WLOXH+CGSNOH+CAPH+CAPSNOH)/2'});
model.geom('geom1').feature('r38').set('createselection', 'on');


model.geom('geom1').create('r39', 'Rectangle');
model.geom('geom1').feature('r39').set('size', {'EGL-6*CGSNOL-2*CAPSNOH' 'CAPSNOH'});
model.geom('geom1').feature('r39').set('pos', {'0' 'EGOXH/3+EGH+CAPH+CGSNOH+CAPSNOH/2'});
model.geom('geom1').feature('r39').set('base', 'center');
model.geom('geom1').feature('r39').set('createselection', 'on');


model.geom('geom1').create('uni8', 'Union');
model.geom('geom1').feature('uni8').selection('input').set({'r29' 'r30' 'r31' 'r32' 'r33' 'r34' 'r35' 'r36' 'r37' 'r38' 'r39'});
model.geom('geom1').feature('uni8').set('intbnd', 'off');
model.geom('geom1').feature('uni8').set('createselection', 'on');


model.geom('geom1').feature.create('unisel2', 'UnionSelection');
model.geom('geom1').feature('unisel2').name('SNO');
model.geom('geom1').feature('unisel2').set('input', {'uni4' 'uni5' 'uni6' 'uni7' 'r9' 'r11' 'r12' 'r13' 'r15' 'r16'});
model.geom('geom1').run('unisel2');


model.geom('geom1').create('r40', 'Rectangle');
model.geom('geom1').feature('r40').set('size', {'MOL/2-(EGL/2+FGL+TNOX+WLOXL)' 'WLCH'});
model.geom('geom1').feature('r40').set('pos', {'(MOL/2+(EGL/2+FGL+TNOX+WLOXL))/2' '-WLCH/2'});
model.geom('geom1').feature('r40').set('base', 'center');
model.geom('geom1').feature('r40').set('createselection', 'on');
model.geom('geom1').run('r40');


model.geom('geom1').create('r41', 'Rectangle');
model.geom('geom1').feature('r41').set('size', {'MOL/2-(EGL/2+FGL+TNOX+WLOXL+WLL)' 'SH'});
model.geom('geom1').feature('r41').set('pos', {'(MOL/2+(EGL/2+FGL+TNOX+WLOXL+WLL))/2' '-SH/2'});
model.geom('geom1').feature('r41').set('base', 'center');
model.geom('geom1').feature('r41').set('createselection', 'on');
model.geom('geom1').run('r41');


model.geom('geom1').create('r42', 'Rectangle');
model.geom('geom1').feature('r42').set('size', {'EGL+2*TNOX+2*FGL/3-2*CGSNOL' 'DH'});
model.geom('geom1').feature('r42').set('pos', {'0' '-DH/2'});
model.geom('geom1').feature('r42').set('base', 'center');
model.geom('geom1').feature('r42').set('createselection', 'on');
model.geom('geom1').run('r42');


model.geom('geom1').feature.duplicate('r43', 'r40');
model.geom('geom1').feature('r43').set('pos', {'-(MOL/2+(EGL/2+FGL+TNOX+WLOXL))/2' '-WLCH/2'});
model.geom('geom1').feature('r43').set('createselection', 'on');
model.geom('geom1').run('r43');


model.geom('geom1').feature.duplicate('r44', 'r41');
model.geom('geom1').feature('r44').set('pos', {'-(MOL/2+(EGL/2+FGL+TNOX+WLOXL+WLL))/2' '-SH/2'});
model.geom('geom1').feature('r44').set('createselection', 'on');
model.geom('geom1').run('r44');


model.geom('geom1').create('fil1', 'Fillet');
model.geom('geom1').feature('fil1').set('radius', 'CAPH/3');
model.geom('geom1').feature('fil1').selection('point').set('uni8', [2 7 8 9 10 12 13 16 17 19 20 21 22 24 25 28 29 34]);
model.geom('geom1').feature('fil1').set('createselection', 'on');


model.geom('geom1').create('fil2', 'Fillet');
model.geom('geom1').feature('fil2').set('radius', '2*CGSNOH/3');
model.geom('geom1').feature('fil2').selection('point').set('uni3', [1 2 3 4 5 6 7 8]);
model.geom('geom1').feature('fil2').selection('point').set('uni2', [1 2 3 4 5 6 7 8]);
model.geom('geom1').feature('fil2').selection('point').set('r10', [1 2]);
model.geom('geom1').feature('fil2').selection('point').set('r14', [1 2]);
model.geom('geom1').feature('fil2').selection('point').set('r19', [2]);
model.geom('geom1').feature('fil2').selection('point').set('r17', [1]);
model.geom('geom1').feature('fil2').selection('point').set('dif1', [1 6]);
model.geom('geom1').feature('fil2').set('createselection', 'on');


model.geom('geom1').create('fil3', 'Fillet');
model.geom('geom1').feature('fil3').set('radius', 'SH');
model.geom('geom1').feature('fil3').selection('point').set('r44', [2]);
model.geom('geom1').feature('fil3').selection('point').set('r41', [1]);
model.geom('geom1').feature('fil3').set('createselection', 'on');


model.geom('geom1').create('fil4', 'Fillet');
model.geom('geom1').feature('fil4').set('radius', 'WLCH');
model.geom('geom1').feature('fil4').selection('point').set('r43', [2]);
model.geom('geom1').feature('fil4').selection('point').set('r40', [1]);
model.geom('geom1').feature('fil4').set('createselection', 'on');


model.geom('geom1').create('fil5', 'Fillet');
model.geom('geom1').feature('fil5').set('radius', 'DH');
model.geom('geom1').feature('fil5').selection('point').set('r42', [1 2]);
model.geom('geom1').feature('fil5').set('createselection', 'on');
model.geom('geom1').run('fil5');


model.geom('geom1').create('r45', 'Rectangle');
model.geom('geom1').feature('r45').set('size', {'MOL/2' 'MOH'});
model.geom('geom1').feature('r45').set('pos', {'-MOL/4' 'MOH/2-100'});
model.geom('geom1').feature('r45').set('base', 'center');
model.geom('geom1').feature('r45').set('createselection', 'on');


model.geom('geom1').create('b1', 'BezierPolygon');
model.geom('geom1').feature('b1').set('degree', [1]);
model.geom('geom1').feature('b1').set('p', {'0' '0'; '0' '0'});
model.geom('geom1').feature('b1').set('w', {'1' '1'});
model.geom('geom1').feature('b1').setIndex('p', 'MOH-100', 1, 0);
model.geom('geom1').feature('b1').setIndex('p', '-100', 1, 1);
model.geom('geom1').feature('b1').set('createselection', 'on');


model.geom('geom1').feature.create('del1', 'Delete');
model.geom('geom1').feature('del1').selection('input').init(1);
model.geom('geom1').feature('del1').selection('input').set('r1', [6]);
model.geom('geom1').feature('del1').selection('input').set('fil5', [1]);
model.geom('geom1').feature('del1').set('createselection', 'on');


model.geom('geom1').create('par1', 'Partition');
model.geom('geom1').feature('par1').selection('tool').set({'b1'});
model.geom('geom1').feature('par1').selection('input').set({'fil1' 'fil2(7)' 'del1(1)' 'del1(2)' 'r4' 'uni1'});
model.geom('geom1').feature('par1').set('createselection', 'on');


model.geom('geom1').create('dif2', 'Difference');
model.geom('geom1').feature('dif2').selection('input').set({'fil2(1)' 'fil2(4)' 'fil2(5)' 'fil3(1)' 'fil4(1)' 'par1' 'r13' 'r15' 'r16' 'r20'  ...
'uni5' 'uni7'});
model.geom('geom1').feature('dif2').selection('input2').set({'r45'});
model.geom('geom1').feature('dif2').set('createselection', 'on');

model.geom('geom1').create('r46', 'Rectangle');
model.geom('geom1').feature('r46').set('size', {'MOL/2-(EGL/2+TNOX+FGL+WLOXL+WLL+CAPSNOH)' 'OXH'});
model.geom('geom1').feature('r46').set('pos', {'(MOL/2+(EGL/2+TNOX+FGL+WLOXL+WLL+CAPSNOH))/2' 'OXH/2'});
model.geom('geom1').feature('r46').set('base', 'center');
model.geom('geom1').feature('r46').set('createselection', 'on');


%% Failure Setting
EGL=178;
TNOX=11;
FGL=108;
FGOX=9.6;
CGSNOL=10;

S=load('C:\D\Program Files\MATLAB\projeects\NVM\New folder\Data\10.mat');
M=S.Results_record{3,6};
s=S.Results_record{3,2};
R=(s(:,1,:)+s(:,2,:)+s(:,3,:))/3;
check_volum=S.Results_record{3,4};
t=length(check_volum(:,1));
[a1,b1]=size(M);
R=reshape(R,a1,b1);
R(R<0.5)=0;
R(R>0.5)=1;
R=R';
U=zeros(a1,b1);
M(M>0.5)=1;
M(M<=0.5)=0;
for i=1:t
    U(check_volum(i,1),check_volum(i,3))=1;
    M(check_volum(i,1),check_volum(i,3))=1;
end
M=M';
U=U';
u1=1:b1;
u2=a1:-1:1;
M=M(u1,u2);
U=U(u1,u2);
hf_name=1;
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
for i=1:a1*b1
   [x,y]=ind2sub([b1,a1],i);
   if 0
       
       model.geom('geom1').feature.create(num2str(hf_name), 'Rectangle');
       model.geom('geom1').lengthUnit('nm');
       model.geom('geom1').feature(num2str(hf_name)).set('size', [0.8,0.8]);
       model.geom('geom1').feature(num2str(hf_name)).set('base', 'center');
       model.geom('geom1').feature(num2str(hf_name)).set('pos', [FGOX+2*FGH/3+(y-0.5-a1/2)*0.8,EG/2+(0.5+b1/2-x)*0.8]);
       model.geom('geom1').feature(num2str(hf_name)).set('type', 'solid');
       model.geom('geom1').run(num2str(hf_name));
       
       if U(x,y,z)==1
          hf_fila_cell(con1)={num2str(hf_name)}; 
          con1=con1+1;
       elseif R(x,y,z)==1&&U(x,y,z)~=1
          hf_gb_cell(con2)={num2str(hf_name)}; 
          con2=con2+1;          
       else
          hf_name_cell(con3)={num2str(hf_name)};
          con3=con3+1;
       end
       
       hf_name=hf_name+1;
       
   elseif M(x,y)==1

       model.geom('geom1').feature.create(num2str(-il_name), 'Rectangle');
       model.geom('geom1').lengthUnit('nm');
       model.geom('geom1').feature(num2str(-il_name)).set('size', [0.8,0.8]);
       model.geom('geom1').feature(num2str(-il_name)).set('base', 'center');
       model.geom('geom1').feature(num2str(-il_name)).set('pos', [EG/2+(x-0.5)*0.8,FGOX+2*FGH/3+(y-0.5-a1/2)*0.8]);
       model.geom('geom1').feature(num2str(-il_name)).set('type', 'solid');
       model.geom('geom1').run(num2str(-il_name));
       
       if U(x,y)==1
           il_fila_cell(con4)={num2str(-il_name)};
           con4=con4+1;       
       else
           il_name_cell(con5)={num2str(-il_name)};
           con5=con5+1;
       end
       
       il_name=il_name+1;
       
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

model.geom('geom1').feature.create('il', 'Union');
model.geom('geom1').feature('il').selection('input').set(il_name_cell);
model.geom('geom1').run('il');

model.geom('geom1').feature.create('ilfila', 'Union');
model.geom('geom1').feature('ilfila').selection('input').set(il_fila_cell);
model.geom('geom1').run('ilfila');

% model.geom('geom1').feature.create('iltrap', 'Union');
% model.geom('geom1').feature('iltrap').selection('input').set(il_trap_cell);
% model.geom('geom1').run('iltrap');

model.geom('geom1').feature('il').set('createselection', 'on');
model.geom('geom1').feature('ilfila').set('createselection', 'on');
% model.geom('geom1').feature('iltrap').set('createselection', 'on');





%% MAT setting

model.material.create('mat1', 'Common', 'comp1');
model.material('mat1').label('Substrate');
model.material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'12.3e-6' '0' '0' '0' '12.3e-6' '0' '0' '0' '12.3e-6'});
model.material('mat1').propertyGroup('def').set('density', '2.648e3');
model.material('mat1').propertyGroup('def').set('electricconductivity', {'1e-14' '0' '0' '0' '1e-14' '0' '0' '0' '1e-14'});
model.material('mat1').propertyGroup('def').set('heatcapacity', '700');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'6' '0' '0' '0' '6' '0' '0' '0' '6'});
model.material('mat1').propertyGroup('def').set('resistivity', {'1e17' '0' '0' '0' '1e17' '0' '0' '0' '1e17'});
model.material('mat1').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat1').propertyGroup('def').set('poissonsratio', '0.17');
model.material('mat1').set('family', 'plastic');


model.material.create('mat2', 'Common', 'comp1');
model.material('mat2').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat2').propertyGroup.create('SemicondMaterial', 'Semiconductor material');
model.material('mat2').label('P Channel Well');
model.material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat2').propertyGroup('def').set('electricconductivity', {'4.52e-2[S/m]' '0' '0' '0' '4.52e-2[S/m]' '0' '0' '0' '4.52e-2[S/m]'});
model.material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]'});
model.material('mat2').propertyGroup('def').set('heatcapacity', '700[J/(kg*K)]');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.material('mat2').propertyGroup('def').set('density', '2329[kg/m^3]');
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]'});
model.material('mat2').propertyGroup('Enu').set('youngsmodulus', '170e9[Pa]');
model.material('mat2').propertyGroup('Enu').set('poissonsratio', '0.28');
model.material('mat2').propertyGroup('RefractiveIndex').set('n', '');
model.material('mat2').propertyGroup('RefractiveIndex').set('ki', '');
model.material('mat2').propertyGroup('RefractiveIndex').set('n', {'3.48' '0' '0' '0' '3.48' '0' '0' '0' '3.48'});
model.material('mat2').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.material('mat2').propertyGroup('SemicondMaterial').set('Eg0', '1.14');
model.material('mat2').propertyGroup('SemicondMaterial').set('chi0', '4.05');
model.material('mat2').propertyGroup('SemicondMaterial').set('Nv', '3.1e25*(T/300)^1.85');
model.material('mat2').propertyGroup('SemicondMaterial').set('Nc', '2.86e25*(T/300)^1.58');
model.material('mat2').propertyGroup('SemicondMaterial').set('mun', '0.1413*(T/300)^-2.6');
model.material('mat2').propertyGroup('SemicondMaterial').set('mup', '0.047*(T/300)^-2.3');
model.material('mat2').propertyGroup('SemicondMaterial').addInput('temperature');
model.material('mat2').set('family', 'plastic');


model.material.create('mat3', 'Common', 'comp1');
model.material('mat3').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat3').propertyGroup.create('SemicondMaterial', 'Semiconductor material');
model.material('mat3').label('WL D4 channel');
model.material('mat3').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').propertyGroup('def').set('electricconductivity', {'1.13e-1[S/m]' '0' '0' '0' '1.13e-1[S/m]' '0' '0' '0' '1.13e-1[S/m]'});
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
model.material('mat3').propertyGroup('SemicondMaterial').set('Eg0', '1.14');
model.material('mat3').propertyGroup('SemicondMaterial').set('chi0', '4.05');
model.material('mat3').propertyGroup('SemicondMaterial').set('Nv', '3.1e25*(T/300)^1.85');
model.material('mat3').propertyGroup('SemicondMaterial').set('Nc', '2.86e25*(T/300)^1.58');
model.material('mat3').propertyGroup('SemicondMaterial').set('mun', '0.1322*(T/300)^-2.6');
model.material('mat3').propertyGroup('SemicondMaterial').set('mup', '0.046*(T/300)^-2.3');
model.material('mat3').propertyGroup('SemicondMaterial').addInput('temperature');
model.material('mat3').set('family', 'plastic');


model.material.create('mat4', 'Common', 'comp1');
model.material('mat4').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat4').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat4').propertyGroup.create('SemicondMaterial', 'Semiconductor material');
model.material('mat4').label('Source');
model.material('mat4').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat4').propertyGroup('def').set('electricconductivity', {'6.4e1[S/m]' '0' '0' '0' '6.4e1[S/m]' '0' '0' '0' '6.4e1[S/m]'});
model.material('mat4').propertyGroup('def').set('thermalexpansioncoefficient', {'2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]'});
model.material('mat4').propertyGroup('def').set('heatcapacity', '700[J/(kg*K)]');
model.material('mat4').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.material('mat4').propertyGroup('def').set('density', '2329[kg/m^3]');
model.material('mat4').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]'});
model.material('mat4').propertyGroup('Enu').set('youngsmodulus', '170e9[Pa]');
model.material('mat4').propertyGroup('Enu').set('poissonsratio', '0.28');
model.material('mat4').propertyGroup('RefractiveIndex').set('n', '');
model.material('mat4').propertyGroup('RefractiveIndex').set('ki', '');
model.material('mat4').propertyGroup('RefractiveIndex').set('n', {'3.48' '0' '0' '0' '3.48' '0' '0' '0' '3.48'});
model.material('mat4').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.material('mat4').propertyGroup('SemicondMaterial').set('Eg0', '1.14');
model.material('mat4').propertyGroup('SemicondMaterial').set('chi0', '4.05');
model.material('mat4').propertyGroup('SemicondMaterial').set('Nv', '3.1e25*(T/300)^1.85');
model.material('mat4').propertyGroup('SemicondMaterial').set('Nc', '2.86e25*(T/300)^1.58');
model.material('mat4').propertyGroup('SemicondMaterial').set('mun', '0.1297*(T/300)^-2.6');
model.material('mat4').propertyGroup('SemicondMaterial').set('mup', '0.045*(T/300)^-2.3');
model.material('mat4').propertyGroup('SemicondMaterial').addInput('temperature');
model.material('mat4').set('family', 'plastic');


model.material.create('mat5', 'Common', 'comp1');
model.material('mat5').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat5').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat5').propertyGroup.create('SemicondMaterial', 'Semiconductor material');
model.material('mat5').label('Drain');
model.material('mat5').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat5').propertyGroup('def').set('electricconductivity', {'4.26e1[S/m]' '0' '0' '0' '4.26e1[S/m]' '0' '0' '0' '4.26e1[S/m]'});
model.material('mat5').propertyGroup('def').set('thermalexpansioncoefficient', {'2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]' '0' '0' '0' '2.6e-6[1/K]'});
model.material('mat5').propertyGroup('def').set('heatcapacity', '700[J/(kg*K)]');
model.material('mat5').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.material('mat5').propertyGroup('def').set('density', '2329[kg/m^3]');
model.material('mat5').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]' '0' '0' '0' '130[W/(m*K)]'});
model.material('mat5').propertyGroup('Enu').set('youngsmodulus', '170e9[Pa]');
model.material('mat5').propertyGroup('Enu').set('poissonsratio', '0.28');
model.material('mat5').propertyGroup('RefractiveIndex').set('n', '');
model.material('mat5').propertyGroup('RefractiveIndex').set('ki', '');
model.material('mat5').propertyGroup('RefractiveIndex').set('n', {'3.48' '0' '0' '0' '3.48' '0' '0' '0' '3.48'});
model.material('mat5').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.material('mat5').propertyGroup('SemicondMaterial').set('Eg0', '1.14');
model.material('mat5').propertyGroup('SemicondMaterial').set('chi0', '4.05');
model.material('mat5').propertyGroup('SemicondMaterial').set('Nv', '3.1e25*(T/300)^1.85');
model.material('mat5').propertyGroup('SemicondMaterial').set('Nc', '2.86e25*(T/300)^1.58');
model.material('mat5').propertyGroup('SemicondMaterial').set('mun', '0.1325*(T/300)^-2.6');
model.material('mat5').propertyGroup('SemicondMaterial').set('mup', '0.045*(T/300)^-2.3');
model.material('mat5').propertyGroup('SemicondMaterial').addInput('temperature');
model.material('mat5').set('family', 'plastic');


model.material.create('mat6', 'Common', 'comp1');
model.material('mat6').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat6').label('PS Gate');
model.material('mat6').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat6').propertyGroup('def').set('electricconductivity', {'1e5[S/m]' '0' '0' '0' '1e5[S/m]' '0' '0' '0' '1e5[S/m]'});
model.material('mat6').propertyGroup('def').set('thermalexpansioncoefficient', {'4.4e-6[1/K]' '0' '0' '0' '4.4e-6[1/K]' '0' '0' '0' '4.4e-6[1/K]'});
model.material('mat6').propertyGroup('def').set('heatcapacity', '678[J/(kg*K)]');
model.material('mat6').propertyGroup('def').set('relpermittivity', {'4.5' '0' '0' '0' '4.5' '0' '0' '0' '4.5'});
model.material('mat6').propertyGroup('def').set('density', '2648[kg/m^3]');
model.material('mat6').propertyGroup('def').set('thermalconductivity', {'22.4[W/(m*K)]' '0' '0' '0' '22.4[W/(m*K)]' '0' '0' '0' '22.4[W/(m*K)]'});
model.material('mat6').propertyGroup('Enu').set('youngsmodulus', '160e9[Pa]');
model.material('mat6').propertyGroup('Enu').set('poissonsratio', '0.22');
model.material('mat6').set('family', 'plastic');


model.material.create('mat7', 'Common', 'comp1');
model.material('mat7').propertyGroup.create('KG', 'Bulk modulus and shear modulus');
model.material('mat7').label('Gate CAP');
model.material('mat7').propertyGroup('def').set('electricconductivity', {'2.941e6' '0' '0' '0' '2.941e6' '0' '0' '0' '2.941e6'});
model.material('mat7').propertyGroup('def').set('density', '7200');
model.material('mat7').propertyGroup('def').set('resistivity', {'3.4e-7' '0' '0' '0' '3.4e-7' '0' '0' '0' '3.4e-7'});
model.material('mat7').propertyGroup('def').set('relpermittivity', {'1.685' '0' '0' '0' '1.685' '0' '0' '0' '1.685'});
model.material('mat7').propertyGroup('def').set('thermalexpansioncoefficient', {'1.6e-5' '0' '0' '0' '1.6e-5' '0' '0' '0' '1.6e-5'});
model.material('mat7').propertyGroup('def').set('thermalconductivity', {'33' '0' '0' '0' '33' '0' '0' '0' '33'});
model.material('mat7').propertyGroup('def').set('youngsmodulus', '1.5e11');
model.material('mat7').propertyGroup('def').set('heatcapacity', '191.739');
model.material('mat7').propertyGroup('KG').set('K', '');
model.material('mat7').propertyGroup('KG').set('G', '');
model.material('mat7').propertyGroup('KG').set('K', '1.406e11');
model.material('mat7').propertyGroup('KG').set('G', '9e10');
model.material('mat7').set('family', 'plastic');


model.material.create('mat8', 'Common', 'comp1');
model.material('mat8').label('SIN');
model.material('mat8').propertyGroup('def').set('thermalexpansioncoefficient', {'3.3e-6' '0' '0' '0' '3.3e-6' '0' '0' '0' '3.3e-6'});
model.material('mat8').propertyGroup('def').set('density', '3.27e3');
model.material('mat8').propertyGroup('def').set('electricconductivity', {'1e-10' '0' '0' '0' '1e-10' '0' '0' '0' '1e-10'});
model.material('mat8').propertyGroup('def').set('heatcapacity', '710.6');
model.material('mat8').propertyGroup('def').set('relpermittivity', {'9.5' '0' '0' '0' '9.5' '0' '0' '0' '9.5'});
model.material('mat8').propertyGroup('def').set('resistivity', {'1e10' '0' '0' '0' '1e10' '0' '0' '0' '1e10'});
model.material('mat8').propertyGroup('def').set('thermalconductivity', {'29' '0' '0' '0' '29' '0' '0' '0' '29'});
model.material('mat8').propertyGroup('def').set('poissonsratio', '0.24');
model.material('mat8').set('family', 'plastic');


model.material.create('mat9', 'Common', 'comp1');
model.material('mat9').label('OX');
model.material('mat9').propertyGroup('def').set('thermalexpansioncoefficient', {'12.3e-6' '0' '0' '0' '12.3e-6' '0' '0' '0' '12.3e-6'});
model.material('mat9').propertyGroup('def').set('density', '2.648e3');
model.material('mat9').propertyGroup('def').set('electricconductivity', {'1e-14' '0' '0' '0' '1e-14' '0' '0' '0' '1e-14'});
model.material('mat9').propertyGroup('def').set('heatcapacity', '700');
model.material('mat9').propertyGroup('def').set('relpermittivity', {'6' '0' '0' '0' '6' '0' '0' '0' '6'});
model.material('mat9').propertyGroup('def').set('resistivity', {'1e14' '0' '0' '0' '1e14' '0' '0' '0' '1e14'});
model.material('mat9').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat9').propertyGroup('def').set('poissonsratio', '0.17');
model.material('mat9').set('family', 'plastic');


model.material.create('mat10', 'Common', 'comp1');
model.material('mat10').label('CAPSIN');
model.material('mat10').propertyGroup('def').set('thermalexpansioncoefficient', {'3.3e-6' '0' '0' '0' '3.3e-6' '0' '0' '0' '3.3e-6'});
model.material('mat10').propertyGroup('def').set('density', '3.27e3');
model.material('mat10').propertyGroup('def').set('electricconductivity', {'1e-10' '0' '0' '0' '1e-10' '0' '0' '0' '1e-10'});
model.material('mat10').propertyGroup('def').set('heatcapacity', '710.6');
model.material('mat10').propertyGroup('def').set('relpermittivity', {'9.5' '0' '0' '0' '9.5' '0' '0' '0' '9.5'});
model.material('mat10').propertyGroup('def').set('resistivity', {'1e10' '0' '0' '0' '1e10' '0' '0' '0' '1e10'});
model.material('mat10').propertyGroup('def').set('thermalconductivity', {'29' '0' '0' '0' '29' '0' '0' '0' '29'});
model.material('mat10').propertyGroup('def').set('poissonsratio', '0.24');
model.material('mat10').set('family', 'plastic');


model.material.create('mat11');
model.material('mat11').label('Interconnects');
model.material('mat11').propertyGroup('def').set('resistivity', {'3.4e-7'});
model.material('mat11').propertyGroup('def').set('density', {'7200'});
model.material('mat11').propertyGroup('def').set('electricconductivity', {'2.941e6'});
model.material('mat11').propertyGroup('def').set('relpermittivity', {'1.685'});
model.material('mat11').propertyGroup.create('KG', 'Bulk modulus and shear modulus');
model.material('mat11').propertyGroup('KG').set('G', {});
model.material('mat11').propertyGroup('KG').set('K', {'1.406e11'});
model.material('mat11').propertyGroup('KG').set('G', {'9e10'});
model.material('mat11').propertyGroup('def').set('thermalconductivity', {'33'});
model.material('mat11').propertyGroup('def').set('youngsmodulus', {'1.5e11'});
model.material('mat11').propertyGroup('def').set('heatcapacity', {'191.739'});


model.material.create('mat12');
model.material('mat12').name('Fila SiO2');
model.material('mat12').propertyGroup('def').set('relpermittivity', {'8.4' '0' '0' '0' '8.4' '0' '0' '0' '8.4'});
model.material('mat12').propertyGroup('def').set('density', '2.648e3');
model.material('mat12').propertyGroup('def').set('heatcapacity', '700');
model.material('mat12').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat12').propertyGroup('def').set('electricconductivity', {'1.2e3' '0' '0' '0' '1.2e3' '0' '0' '0' '1.2e3'});
model.material('mat12').set('family', 'plastic');

model.material.create('mat13');
model.material('mat13').name('Post SiO2');
model.material('mat13').propertyGroup('def').set('relpermittivity', {'8.4' '0' '0' '0' '8.4' '0' '0' '0' '8.4'});
model.material('mat13').propertyGroup('def').set('density', '2.648e3');
model.material('mat13').propertyGroup('def').set('heatcapacity', '700');
model.material('mat13').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
model.material('mat13').propertyGroup('def').set('electricconductivity', {'1.2e-3' '0' '0' '0' '1.2e-3' '0' '0' '0' '1.2e-3'});
model.material('mat13').set('family', 'plastic');


% model.material.create('mat13');
% model.material('mat13').name('Trap SiO2');
% model.material('mat13').propertyGroup('def').set('relpermittivity', {'8.4' '0' '0' '0' '8.4' '0' '0' '0' '8.4'});
% model.material('mat13').propertyGroup('def').set('density', '2.648e3');
% model.material('mat13').propertyGroup('def').set('heatcapacity', '700');
% model.material('mat13').propertyGroup('def').set('thermalconductivity', {'1.3' '0' '0' '0' '1.3' '0' '0' '0' '1.3'});
% model.material('mat13').propertyGroup('def').set('electricconductivity', {'1.2e-5' '0' '0' '0' '1.2e-5' '0' '0' '0' '1.2e-5'});
% model.material('mat13').set('family', 'plastic');


model.geom('geom1').run;


model.material('mat1').selection.set([1]);
model.material('mat2').selection.set([2]);
model.material('mat3').selection.set([16]);
model.material('mat4').selection.set([3]);
model.material('mat5').selection.set([20]);
model.material('mat6').selection.set([5 11 14 17]);
model.material('mat7').selection.set([6 18]);
model.material('mat8').selection.set([10 12 13 15 19]);
model.material('mat9').selection.set([4 8]);
model.material('mat10').selection.set([7]);
model.material('mat11').selection.set([9 21]);
model.material('mat13').selection.named('geom1_il_dom');
model.material('mat12').selection.named('geom1_ilfila_dom');
% model.material('mat13').selection.named('geom1_iltrap_dom');

%% Testing setting (Program)


model.physics.create('semi', 'Semiconductor', 'geom1');
model.physics('semi').selection.set([2 3 4 7 16 20]);
model.physics('semi').feature.create('ccn1', 'ChargeConservation', 2);
model.physics('semi').feature('ccn1').selection.set([4 7]);
model.physics('semi').feature.create('adm1', 'AnalyticDopingModel', 2);
model.physics('semi').feature('adm1').label('P-Well doping');
model.physics('semi').feature('adm1').set('NAc', '6e12[1/cm^3]');
% model.physics('semi').feature('adm1').selection.set([2 3 16 20]);

model.physics('semi').feature.create('adm2', 'AnalyticDopingModel', 2);
model.physics('semi').feature('adm2').label('Source doping');
model.physics('semi').feature('adm2').set('impurityType', 'donor');
model.physics('semi').feature('adm2').set('NDc', '3.08e15[1/cm^3]');
% model.physics('semi').feature('adm2').selection.set([3]);

model.physics('semi').feature.create('adm3', 'AnalyticDopingModel', 2);
model.physics('semi').feature('adm3').label('Drain doping');
model.physics('semi').feature('adm3').set('impurityType', 'donor');
model.physics('semi').feature('adm3').set('NDc', '2e15[1/cm^3]');
% model.physics('semi').feature('adm3').selection.set([20]);

model.physics('semi').feature.create('adm4', 'AnalyticDopingModel', 2);
model.physics('semi').feature('adm4').label('WLCH doping');
model.physics('semi').feature('adm4').set('NAc', '9e12[1/cm^3]');
% model.physics('semi').feature('adm4').selection.set([16]);

%model.physics('semi').feature.create('gdm1', 'GeometricDopingModel', 2);
%model.physics('semi').feature('gdm1').label('Source doping');
%model.physics('semi').feature('gdm1').set('impurityType', 'donor');
%model.physics('semi').feature('gdm1').set('NDgen', '3.08e15[1/cm^3]');
%model.physics('semi').feature('gdm1').set('JunctionOrLength', 'decay_length');
%model.physics('semi').feature('gdm1').set('l_gen', '0.01[um]');
%model.physics('semi').feature('gdm1').selection.set([3]);
%model.physics('semi').feature('gdm1').feature('gdmbs1').selection.set([5 6 32 104 111]);

%model.physics('semi').feature.create('gdm2', 'GeometricDopingModel', 2);
%model.physics('semi').feature('gdm2').label('Drain doping');
%model.physics('semi').feature('gdm2').set('impurityType', 'donor');
%model.physics('semi').feature('gdm2').set('NDgen', '2e15[1/cm^3]');
%model.physics('semi').feature('gdm2').set('JunctionOrLength', 'decay_length');
%model.physics('semi').feature('gdm2').set('l_gen', '0.02[um]');
%model.physics('semi').feature('gdm2').selection.set([20]);
%model.physics('semi').feature('gdm2').feature('gdmbs1').selection.set([88 94 97 101 126]);

%model.physics('semi').feature.create('gdm3', 'GeometricDopingModel', 2);
%model.physics('semi').feature('gdm3').label('WLCH doping');
%model.physics('semi').feature('gdm3').set('impurityType', 'acceptor');
%model.physics('semi').feature('gdm3').set('NAgen', '9e12[1/cm^3]');
%model.physics('semi').feature('gdm3').set('JunctionOrLength', 'decay_length');
%model.physics('semi').feature('gdm3').set('l_gen', '0.005[um]');
%model.physics('semi').feature('gdm3').selection.set([16]);
%model.physics('semi').feature('gdm3').feature('gdmbs1').selection.set([66 86 97 100 122 126]);

model.physics('semi').feature.create('mc1', 'MetalContact', 1);
model.physics('semi').feature('mc1').label('Source voltage');
model.physics('semi').feature('mc1').set('V0', '4.5[V]');
% model.physics('semi').feature('mc1').selection.set([5]);

model.physics('semi').feature.create('mc2', 'MetalContact', 1);
model.physics('semi').feature('mc2').label('Drain setting');
model.physics('semi').feature('mc2').set('V0', '0.1[V]');
% model.physics('semi').feature('mc2').selection.set([101]);

model.physics('semi').feature.create('mc3', 'MetalContact', 1);
model.physics('semi').feature('mc3').label('Base');
model.physics('semi').feature('mc3').set('V0', '0.0[V]');
% model.physics('semi').feature('mc3').selection.set([4]);

model.physics('semi').feature.create('term1', 'Terminal', 1);
model.physics('semi').feature('term1').set('TerminalType', 'Voltage');
model.physics('semi').feature('term1').set('Phic', '4.05[V]');
model.physics('semi').feature('term1').label('CG');
model.physics('semi').feature('term1').set('V0', '10.5');
% model.physics('semi').feature('term1').selection.set([50 51 53 57 116 118]);

model.physics('semi').feature.create('term2', 'Terminal', 1);
model.physics('semi').feature('term2').set('TerminalType', 'Voltage');
model.physics('semi').feature('term2').set('Phic', '4.05[V]');
model.physics('semi').feature('term2').label('EG');
model.physics('semi').feature('term2').set('V0', '4.5');
% model.physics('semi').feature('term2').selection.set([12 33 36 37 105 110]);

model.physics('semi').feature.create('term3', 'Terminal', 1);
model.physics('semi').feature('term3').set('TerminalType', 'Voltage');
model.physics('semi').feature('term3').set('Phic', '4.05[V]');
model.physics('semi').feature('term3').label('WL');
model.physics('semi').feature('term3').set('V0', '0.7');
% model.physics('semi').feature('term3').selection.set([67 68 70 74 89 90 123]);

model.physics('semi').feature.create('ii2', 'InsulatorInterface', 1);
model.physics('semi').feature('ii2').label('Tunneling OX');
% model.physics('semi').feature('ii2').selection.set([54]);
model.physics('semi').feature('ii2').set('TunnelingType', 'FowlerNordheim');
model.physics('semi').feature('ii2').set('IncludeTraps', true);
model.physics('semi').feature('ii2').set('SpecifyNeutralEnergyLevel', 'FromConduction');
model.physics('semi').feature('ii2').set('E0in', '0.5[V]');

model.physics('semi').feature.create('fg1', 'FloatingGate', 1);
% model.physics('semi').feature('fg1').selection.set([38 39 40 46 52 55 59 63 112 113 114 115 117 119 120 121]);


model.physics.create('ec', 'ConductiveMedia', 'geom1');

model.physics('ec').feature.create('pot1', 'ElectricPotential', 1);
model.physics('ec').feature('pot1').label('WL');
% model.physics('ec').feature('pot1').selection.set([67 69 74 89 123]);
model.physics('ec').feature('pot1').set('V0', '0.7');

model.physics('ec').feature.create('pot2', 'ElectricPotential', 1);
model.physics('ec').feature('pot2').label('CG');
% model.physics('ec').feature('pot2').selection.set([50 51 53 57 116 118]);
model.physics('ec').feature('pot2').set('V0', '10.5');

model.physics('ec').feature.create('pot3', 'ElectricPotential', 1);
model.physics('ec').feature('pot3').label('Source');
% model.physics('ec').feature('pot3').selection.set([5]);
model.physics('ec').feature('pot3').set('V0', '4.5');

model.physics('ec').feature.create('pot4', 'ElectricPotential', 1);
model.physics('ec').feature('pot4').label('Ground');
% model.physics('ec').feature('pot4').selection.set([2]);
model.physics('ec').feature('pot4').set('V0', '0');

model.physics('ec').feature.create('pot5', 'ElectricPotential', 1);
model.physics('ec').feature('pot5').label('EG');
% model.physics('ec').feature('pot5').selection.set([8 10 33 36 105 110]);
model.physics('ec').feature('pot5').set('V0', '4.5');

model.physics('ec').feature.create('bcs1', 'BoundaryCurrentSource', 1);
model.physics('ec').feature('bcs1').set('Qjs', '-5e5');
% model.physics('ec').feature('bcs1').selection.set([97 126]);


model.physics.create('ht', 'HeatTransfer', 'geom1');

model.physics('ht').feature.create('open1', 'OpenBoundary', 1);
% model.physics('ht').feature('open1').selection.set([2 15 17 19 102 103]);
model.physics('ht').feature.create('hs1', 'HeatSource', 2);
model.physics('ht').feature('hs1').selection.all;
model.physics('ht').feature('hs1').set('Q0_src', 'root.comp1.semi.Q_tot');


model.multiphysics.create('emh1', 'ElectromagneticHeatSource', 'geom1', 2);

model.multiphysics('emh1').selection.all;





model.geom('geom1').run;

mphgeom(model);




model.label('Two-D-NVM.mph');

out = model;
