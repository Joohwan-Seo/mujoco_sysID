function Js_mat = Je(q1,q2,q3,q4,q5,q6)
%JE
%    JS_MAT = JE(Q1,Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    02-Mar-2023 12:49:24

t2 = conj(q2);
t3 = conj(q3);
t4 = conj(q4);
t5 = conj(q5);
t6 = conj(q6);
t7 = cos(q1);
t8 = cos(q2);
t9 = cos(q3);
t10 = cos(q4);
t11 = cos(q5);
t12 = cos(q6);
t13 = sin(q1);
t14 = sin(q2);
t15 = sin(q3);
t16 = sin(q4);
t17 = sin(q5);
t18 = sin(q6);
t29 = -q3;
t34 = pi./2.0;
t19 = cos(t2);
t20 = cos(t3);
t21 = cos(t4);
t22 = cos(t5);
t23 = cos(t6);
t24 = sin(t2);
t25 = sin(t3);
t26 = sin(t4);
t27 = sin(t5);
t28 = sin(t6);
t30 = t7.*t10;
t31 = t7.*t16;
t32 = t10.*t13;
t33 = t13.*t16;
t35 = q2+t29;
t36 = -t34;
t47 = t11.*(7.0./2.0e+2);
t48 = t11.*(7.0./4.0e+2);
t50 = (t10.*t12)./2.0e+1;
t51 = (t10.*t18)./2.0e+1;
t52 = (t16.*t17)./2.0e+1;
t63 = t9.*t11.*(1.1e+1./2.5e+1);
t64 = t10.*t12.*(7.0./4.0e+2);
t69 = t10.*t17.*(2.1e+1./5.0e+1);
t70 = t12.*t16.*(2.1e+1./5.0e+1);
t71 = t12.*t17.*(7.0./2.0e+2);
t76 = t16.*t18.*(2.1e+1./5.0e+1);
t77 = t17.*t18.*(7.0./2.0e+2);
t78 = t16.*t17.*(7.0./4.0e+2);
t82 = t10.*t12.*(2.3e+1./1.0e+2);
t85 = t10.*t18.*(2.3e+1./1.0e+2);
t86 = t12.*t17.*(2.3e+1./1.0e+2);
t95 = (t11.*t12.*t16)./2.0e+1;
t96 = (t11.*t16.*t18)./2.0e+1;
t99 = t14.*t16.*t17.*(1.1e+1./2.5e+1);
t100 = t15.*t16.*t18.*(1.1e+1./2.5e+1);
t112 = t10.*t11.*t12.*(2.1e+1./5.0e+1);
t113 = t10.*t12.*t14.*(1.1e+1./2.5e+1);
t114 = t9.*t12.*t17.*(1.1e+1./2.5e+1);
t115 = t10.*t11.*t18.*(2.1e+1./5.0e+1);
t119 = t10.*t14.*t18.*(1.1e+1./2.5e+1);
t120 = t10.*t15.*t17.*(1.1e+1./2.5e+1);
t121 = t12.*t15.*t16.*(1.1e+1./2.5e+1);
t122 = t9.*t17.*t18.*(1.1e+1./2.5e+1);
t127 = t11.*t12.*t16.*(2.3e+1./1.0e+2);
t128 = t9.*t11.*t14.*(-7.0./4.0e+2);
t130 = t11.*t16.*t18.*(2.3e+1./1.0e+2);
t134 = t8.*t9.*t10.*t12.*(2.1e+1./5.0e+1);
t135 = t10.*t11.*t12.*t15.*(1.1e+1./2.5e+1);
t136 = t8.*t9.*t10.*t18.*(2.1e+1./5.0e+1);
t137 = t8.*t10.*t12.*t15.*(7.0./2.0e+2);
t138 = t9.*t10.*t12.*t14.*(7.0./2.0e+2);
t139 = t8.*t9.*t10.*t17.*(7.0./4.0e+2);
t140 = t8.*t9.*t12.*t16.*(7.0./4.0e+2);
t141 = t11.*t12.*t14.*t16.*(1.1e+1./2.5e+1);
t142 = t10.*t11.*t15.*t18.*(1.1e+1./2.5e+1);
t143 = t8.*t9.*t16.*t17.*(2.1e+1./5.0e+1);
t144 = t10.*t12.*t14.*t15.*(2.1e+1./5.0e+1);
t145 = t8.*t10.*t15.*t18.*(7.0./2.0e+2);
t146 = t9.*t10.*t14.*t18.*(7.0./2.0e+2);
t147 = t11.*t14.*t16.*t18.*(1.1e+1./2.5e+1);
t148 = t10.*t14.*t15.*t18.*(2.1e+1./5.0e+1);
t149 = t8.*t15.*t16.*t17.*(7.0./2.0e+2);
t150 = t9.*t14.*t16.*t17.*(7.0./2.0e+2);
t151 = t10.*t14.*t15.*t17.*(7.0./4.0e+2);
t152 = t12.*t14.*t15.*t16.*(7.0./4.0e+2);
t153 = t8.*t15.*t17.*t18.*(7.0./4.0e+2);
t154 = t9.*t14.*t17.*t18.*(7.0./4.0e+2);
t156 = t14.*t15.*t16.*t17.*(2.1e+1./5.0e+1);
t160 = t8.*t9.*t12.*t16.*(2.3e+1./1.0e+2);
t162 = t8.*t9.*t16.*t18.*(2.3e+1./1.0e+2);
t167 = t12.*t14.*t15.*t16.*(2.3e+1./1.0e+2);
t168 = t8.*t15.*t17.*t18.*(2.3e+1./1.0e+2);
t169 = t9.*t14.*t17.*t18.*(2.3e+1./1.0e+2);
t172 = t14.*t15.*t16.*t18.*(2.3e+1./1.0e+2);
t184 = t8.*t12.*t15.*t17.*(-2.3e+1./1.0e+2);
t190 = t8.*t9.*t11.*t16.*t18.*(-2.1e+1./5.0e+1);
t192 = t8.*t11.*t12.*t15.*t16.*(-7.0./2.0e+2);
t194 = t9.*t11.*t14.*t16.*t18.*(-7.0./2.0e+2);
t195 = t11.*t14.*t15.*t16.*t18.*(-2.1e+1./5.0e+1);
t37 = cos(t35);
t38 = sin(t35);
t39 = t20.*t22;
t40 = t21.*t23;
t41 = t22.*t25;
t42 = t21.*t28;
t43 = t23.*t26;
t44 = t11.*t31;
t45 = t26.*t28;
t46 = t11.*t33;
t49 = q2+t36;
t53 = t8.*t15.*t31;
t54 = t9.*t14.*t31;
t55 = t8.*t15.*t32;
t56 = t9.*t14.*t32;
t57 = t8.*t15.*t33;
t58 = t9.*t14.*t33;
t65 = t20.*t21.*t27;
t66 = t20.*t23.*t27;
t72 = t21.*t25.*t27;
t73 = t20.*t27.*t28;
t74 = t23.*t25.*t27;
t79 = t25.*t27.*t28;
t83 = -t64;
t84 = -t69;
t87 = -t76;
t88 = t8.*t11.*t15.*t30;
t89 = t9.*t11.*t14.*t30;
t102 = -t86;
t116 = t8.*t15.*t48;
t117 = t9.*t14.*t48;
t118 = -t96;
t123 = t16.*t18.*t48;
t124 = -t100;
t129 = -t120;
t159 = -t130;
t161 = -t137;
t163 = t8.*t15.*t86;
t164 = t9.*t14.*t86;
t165 = -t145;
t166 = -t147;
t170 = -t149;
t171 = -t153;
t173 = t8.*t9.*t11.*t70;
t174 = t8.*t9.*t10.*t18.*t48;
t175 = t8.*t9.*t11.*t76;
t176 = t8.*t12.*t15.*t16.*t47;
t177 = t9.*t12.*t14.*t16.*t47;
t178 = t11.*t14.*t15.*t70;
t179 = t8.*t15.*t16.*t18.*t47;
t180 = t9.*t14.*t16.*t18.*t47;
t181 = t10.*t14.*t15.*t18.*t48;
t182 = t11.*t14.*t15.*t76;
t183 = -t162;
t185 = -t168;
t186 = -t172;
t188 = t8.*t9.*t11.*t82;
t189 = t8.*t9.*t11.*t85;
t191 = t11.*t14.*t15.*t82;
t193 = t11.*t14.*t15.*t85;
t256 = t70+t77+t78+t85+t115+t127;
t59 = cos(t49);
t60 = -t44;
t61 = sin(t49);
t62 = t22.*t40;
t67 = t22.*t42;
t68 = t22.*t43;
t75 = t22.*t45;
t80 = t11.*t38;
t81 = t17.*t38;
t90 = t11.*t55;
t91 = t11.*t56;
t92 = -t53;
t93 = -t55;
t94 = -t58;
t98 = t10.*t11.*t37;
t101 = -t66;
t103 = t7.*t17.*t37;
t104 = t11.*t13.*t37;
t105 = t10.*t17.*t37;
t106 = t12.*t16.*t37;
t107 = -t72;
t109 = t13.*t17.*t37;
t110 = t16.*t18.*t37;
t111 = -t79;
t125 = -t89;
t196 = t41+t65;
t202 = t48+t102;
t241 = t47+t83+t84+t123;
t252 = t71+t82+t87+t112+t159;
t270 = t121+t122+t142+t256;
t271 = t52+t99+t140+t143+t150+t152+t154+t156+t170+t171+t174+t181;
t272 = t51+t95+t119+t136+t141+t146+t148+t160+t165+t167+t169+t173+t177+t178+t185+t189+t192+t193;
t273 = t50+t113+t116+t118+t128+t134+t138+t139+t144+t151+t161+t164+t166+t179+t183+t184+t186+t188+t190+t191+t194+t195;
t97 = -t62;
t108 = -t75;
t126 = -t91;
t131 = t7.*t15.*t61;
t132 = -t105;
t133 = -t110;
t155 = t7.*t9.*t59;
t157 = t7.*t9.*t61;
t158 = t7.*t15.*t59;
t197 = t42+t68;
t198 = t43+t67;
t199 = t39+t107;
t203 = t19.*t196;
t206 = t81+t98;
t212 = t31+t56+t93;
t213 = t32+t54+t92;
t214 = t30+t57+t94;
t242 = t46+t88+t103+t125;
t246 = -t18.*(t44-t90+t91-t109);
t250 = -t12.*(t44-t90+t91-t109);
t251 = t12.*(t44-t90+t91-t109);
t255 = t63+t129+t241;
t266 = t114+t124+t135+t252;
t187 = -t158;
t200 = t40+t108;
t201 = t45+t97;
t204 = t20.*t198;
t205 = t25.*t198;
t208 = t24.*t199;
t210 = t80+t132;
t211 = t18.*t206;
t215 = t12.*t206;
t217 = t12.*t213;
t218 = t131+t155;
t219 = t17.*t212;
t220 = t12.*t214;
t221 = t18.*t213;
t222 = t18.*t214;
t243 = t60+t90+t109+t126;
t244 = t12.*t242;
t245 = t18.*t242;
t207 = t20.*t201;
t209 = t25.*t201;
t216 = -t208;
t223 = t157+t187;
t224 = t10.*t218;
t225 = t26.*t27.*t210;
t226 = t73+t205;
t230 = t111+t204;
t234 = -t19.*(t79-t204);
t235 = t106+t211;
t238 = t133+t215;
t239 = t104+t219;
t249 = -t245;
t253 = -t197.*(t110-t215);
t257 = t221+t244;
t261 = t220+t246;
t263 = t222+t251;
t227 = -t224;
t228 = t11.*t223;
t229 = t74+t207;
t231 = t24.*t226;
t232 = t101+t209;
t237 = -t24.*(t66-t209);
t247 = t203+t216;
t248 = t26.*t27.*t239;
t254 = t200.*t235;
t260 = t217+t249;
t265 = t197.*t257;
t268 = t197.*t263;
t269 = t200.*t261;
t233 = t19.*t229;
t236 = t33+t227;
t258 = t231+t234;
t267 = t200.*t260;
t240 = t17.*t236;
t262 = t233+t237;
t259 = t228+t240;
t264 = t26.*t27.*t259;
mt1 = [-t257.*t272-t259.*t271-t260.*t273,t239.*t271+t261.*t273+t263.*t272,-t210.*t271-t235.*t273-t272.*(t110-t215),-t247.*t259+t258.*t260+t257.*t262,t239.*t247-t258.*t261-t262.*t263,-t210.*t247+t235.*t258+t262.*(t110-t215),-t255.*t259+t257.*t266-t260.*t270,t239.*t255-t263.*t266+t261.*t270,-t210.*t255-t235.*t270+t266.*(t110-t215),-t264-t265-t267,t248+t268+t269,-t225+t253-t254,t241.*t259-t252.*t257+t256.*t260,-t239.*t241+t252.*t263-t256.*t261];
mt2 = [t210.*t241+t235.*t256-t252.*(t110-t215),t264+t265+t267,-t248-t268-t269,t225+t254+t197.*(t110-t215),t202.*t260-t17.*t18.*t257.*(2.3e+1./1.0e+2)-t17.*t18.*t259.*(7.0./4.0e+2),-t202.*t261+t17.*t18.*t239.*(7.0./4.0e+2)+t17.*t18.*t263.*(2.3e+1./1.0e+2),t202.*t235-t17.*t18.*t210.*(7.0./4.0e+2)-t17.*t18.*(t110-t215).*(2.3e+1./1.0e+2),t22.*t259-t23.*t27.*t257+t27.*t28.*t260,-t22.*t239+t23.*t27.*t263-t27.*t28.*t261,t22.*t210+t27.*t28.*t235-t23.*t27.*(t110-t215)];
mt3 = [t12.*t257.*(-2.3e+1./1.0e+2)-t12.*t259.*(7.0./4.0e+2)+t18.*t260.*(2.3e+1./1.0e+2),t12.*t239.*(7.0./4.0e+2)+t12.*t263.*(2.3e+1./1.0e+2)-t18.*t261.*(2.3e+1./1.0e+2),t12.*t210.*(-7.0./4.0e+2)+t18.*t235.*(2.3e+1./1.0e+2)-t12.*(t110-t215).*(2.3e+1./1.0e+2),t23.*t260+t28.*t257,-t23.*t261-t28.*t263,t23.*t235+t28.*(t110-t215),t217.*(7.0./4.0e+2)-t245.*(7.0./4.0e+2),t220.*(-7.0./4.0e+2)+t18.*(t44-t90+t91-t109).*(7.0./4.0e+2),t106.*(7.0./4.0e+2)+t211.*(7.0./4.0e+2),t259,-t239,t210];
Js_mat = reshape([mt1,mt2,mt3],6,6);
