#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

/* Generalised timing model for tempo2                 
 *
 * Based on the DD model, but includes
 *
 *  conversion to ELL1 model if EPS1 and EPS2 are set
 *  BT model                                          (set allTerms = 0)
 *  jumps from BTJ model
 *  use of SHAPMAX (i.e. DDS model)
 *  extra terms implemented in DDK model
 *  flag to convert to DDGR model
 *  multiple binary systems (e.g. planetary systems)
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"


double get_sum(int Num,double*array){
    double sum=0;
    int i;
    for(i=0;i<Num;i++){
        sum+=array[i];
    }
    return sum;
}

#define OUTER_CORR_NODES 16384

/* 外轨道修正表格 */
static struct {
    double time;          /* 相对于 T0_2 的时间 (秒) */
    double x2;            /* 修正后的投影半长轴 (lt-s) */
    double e2;            /* 修正后的偏心率 */
    double om2;           /* 修正后的近星点经度 (rad) */
    double dphase;        /* 相位偏移 δM2/(2π) */
} outer_corr_table[OUTER_CORR_NODES];
static int outer_nodes = 0;
static int outer_builded = 0;

#define INNER_CORR_NODES 16384

static struct {
    double time;      // relative to TASC (sec)
    double dx;        // Δx1 cumulative (lt-s)
    double dorbits;   // ΔΦ1 cumulative (orbits)
} inner_corr_table[INNER_CORR_NODES];
static int inner_nodes = 0;
static int inner_builded = 0;

static void build_outer_corr_table(pulsar *psr, int p, int com_inner, int com_outer);
static void interp1_outer(double t, double *x2, double *e2, double *om2, double *dphase);

static void build_inner_corr_table(pulsar *psr, int p, int com_inner, int com_outer);
static void interp1_inner(double t, double *dx, double *dorbits);

/* ---------- 外轨道方向投影计算 ---------- */
static void compute_unit_vectors(double i1, double i2, double dOmega,
                                 double u,   /* omega2 + theta2 */
                                 double *r_ex, double *r_ey, double *r_ez,
                                 double *t_ex, double *t_ey, double *t_ez)
{
    double cosOm = cos(dOmega), sinOm = sin(dOmega);
    double cosi1 = cos(i1), sini1 = sin(i1);
    double cosi2 = cos(i2), sini2 = sin(i2);

    /* 外轨道基矢 (天空系) */
    double l2_x = cosOm, l2_y = sinOm, l2_z = 0.0;          // 升交点方向
    double n2_x = sini2 * sinOm, n2_y = -sini2 * cosOm, n2_z = cosi2; // 法向
    double m2_x = -cosi2 * sinOm, m2_y = cosi2 * cosOm, m2_z = sini2; // n2 x l2

    /* 径向和横向单位矢 (天空系) */
    double r_x = cos(u) * l2_x + sin(u) * m2_x;
    double r_y = cos(u) * l2_y + sin(u) * m2_y;
    double r_z = cos(u) * l2_z + sin(u) * m2_z;

    double t_x = -sin(u) * l2_x + cos(u) * m2_x;
    double t_y = -sin(u) * l2_y + cos(u) * m2_y;
    double t_z = -sin(u) * l2_z + cos(u) * m2_z;

    /* 投影到内轨道坐标系 (Ω₁ = 0) */
    // e_x = (1,0,0)
    *r_ex = r_x;
    // e_y = (0, cosi1, sini1)
    *r_ey = cosi1 * r_y + sini1 * r_z;
    // e_z = (0, -sini1, cosi1)
    *r_ez = -sini1 * r_y + cosi1 * r_z;

    *t_ex = t_x;
    *t_ey = cosi1 * t_y + sini1 * t_z;
    *t_ez = -sini1 * t_y + cosi1 * t_z;
    double norm = (*r_ex)*(*r_ex) + (*r_ey)*(*r_ey) + (*r_ez)*(*r_ez);
    if (fabs(norm - 1.0) > 1e-6) {
        printf("WARNING: r_hat not normalized! norm = %.6f\n", norm);
        printf("  r_ex = %.6f, r_ey = %.6f, r_ez = %.6f\n", *r_ex, *r_ey, *r_ez);
    }
}

double get_theta2(pulsar *psr, int p, double tt0,int com1, int com2){
    double pb2=psr[p].param[param_pb].val[com2]*SECDAY;
    double ecc2=psr[p].param[param_ecc].val[com2];

    tt0+=(psr[p].param[param_tasc].val[com1]-psr[p].param[param_t0].val[com2])*SECDAY;
    double orbits = tt0/pb2;

    int norbits = (int)orbits;
    if (orbits < 0.0) norbits--;

    double phase = 2.0*M_PI * (orbits-norbits);

    /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
    double ep = phase + ecc2*sin(phase)*(1.0+ecc2*cos(phase));

    /* This line is wrong in the original tempo: should be inside the do loop */
    /*  denom = 1.0 - ecc*cos(ep);*/
    double dep;
    do {
        dep = (phase - (ep-ecc2*sin(ep)))/(1.0 - ecc2*cos(ep));
        ep += dep;
    } while (fabs(dep) > 1.0e-12);

    return 2.0*atan(sqrt((1+ecc2)/(1-ecc2))*tan(0.5*ep));
}

/* ---------- 外轨道完整状态（含径向/横向速度） ---------- */
static void compute_outer_state_and_vel(double i1_rad, double i2_rad, double dOmega_rad,
                                        double e2, double omega2_rad, double theta2_rad,
                                        double Porb2_day, double Mt,
                                        double *S2_lt, double *theta_dot,
                                        double *r_ex, double *r_ey, double *r_ez,
                                        double *t_ex, double *t_ey, double *t_ez,
                                        double *vr, double *vt)
{
    double P_sec = Porb2_day * SECDAY;
    double a2_lt = pow(GM_C3 * Mt * P_sec * P_sec / (4.0 * M_PI * M_PI), 1.0/3.0);
    double k = (1.0 - e2 * e2) / (1.0 + e2 * cos(theta2_rad));
    *S2_lt = a2_lt * k;
    *theta_dot = 2.0 * M_PI / P_sec * pow(1.0 + e2 * cos(theta2_rad), 2.0)
                 / pow(1.0 - e2 * e2, 1.5);

    double u = omega2_rad + theta2_rad;
    compute_unit_vectors(i1_rad, i2_rad, dOmega_rad, u, r_ex, r_ey, r_ez, t_ex, t_ey, t_ez);

    /* 径向与横向速度 */
    double vel_factor = sqrt(GM_C3 * Mt / (a2_lt * (1.0 - e2 * e2)));
    *vr = vel_factor * e2 * sin(theta2_rad);
    *vt = vel_factor * (1.0 + e2 * cos(theta2_rad));
}

/* ---------- 倾角变化率引起的 dot{x}_1 ---------- */
static double inclination_xdot_from_state(double i1_rad, double Mstar, double Mwd, double Mt,
                                          double Porbi_day, double a1_lt,
                                          double S2_lt, double r_ex, double r_ey, double r_ez)
{
    double Min = Mt - Mstar;  
    double mass_factor = Mstar * Min * Min / (Mwd * Mwd * Mwd);
    double PORB1 = Porbi_day * SECDAY;
    double di1_dt = (3.0 * M_PI / PORB1) * mass_factor 
                    * pow(a1_lt / S2_lt, 3.0) * (r_ex * r_ez);
    return a1_lt * cos(i1_rad) * di1_dt;
}

/* ---------- 四极矩平均势引起的 f_orb1（接收预计算状态） ---------- */
static double quadrupole_forb_from_state(double i1_rad, double Mstar, double Mwd, double Mt,
                                         double Porbi_day, double a1_lt,
                                         double S2_lt, double r_ex, double r_ey, double r_ez)
{
    double Min = Mt - Mstar;  /* 同上，Mt 可参数化 */
    double a_in_lt = a1_lt * Min / Mwd;
    double FORB1 = 1.0 / (Porbi_day * SECDAY);
    double ratio = a_in_lt / S2_lt;
    double Q = (3.0 * r_ez * r_ez - 1.0) - 1.5 / tan(i1_rad) * r_ez * r_ey;
    return FORB1 * (Mstar / Min) * pow(ratio, 3.0) * Q;
}

/* ---------- 横向多普勒效应（接收预计算状态） ---------- */
static void transverse_doppler_from_state(double i1_rad, double Mstar, double Mt,
                                          double a1_lt,
                                          double r_ex, double r_ey, double t_ex, double t_ey,
                                          double vr, double vt,
                                          double *delta_x, double *delta_phi)
{
    double v_rel_x = vr * r_ex + vt * t_ex;
    double v_rel_y = vr * r_ey + vt * t_ey;
    double q = Mstar / Mt;
    double v2x = v_rel_x * q;
    double v2y = v_rel_y * q;

    *delta_x   = a1_lt * v2y;
    *delta_phi = v2x / sin(i1_rad);
}

/* ---------- 主修正函数（消除重复计算） ---------- */
void get_correction(pulsar *psr, int p, double tt0, int com1, int com2,
                    double *orbits, double *x, double *dRomer, double *mc1, double *m_in)
{
    /* ----- 提取固定参数 ----- */
    double i1 = psr[p].param[param_a0].val[0] * M_PI / 180.0;
    double i2 = psr[p].param[param_b0].val[0] * M_PI / 180.0;
    double dOmega = psr[p].param[param_kom].val[0] * M_PI / 180.0;
    double e2     = psr[p].param[param_ecc].val[com2];
    double omega2 = psr[p].param[param_om].val[com2] * M_PI / 180.0;
    double Mt     = psr[p].param[param_mtot].val[0];
    double Porbo_day = psr[p].param[param_pb].val[com2];

    double x1_lt  = psr[p].param[param_a1].val[com1];
    double x2_lt  = psr[p].param[param_a1].val[com2];
    double Porbi_day = psr[p].param[param_pb].val[com1];

    /* 推算质量 */
    double Mstar = pow(pow(x2_lt, 3.0) * pow(Porbo_day * SECDAY, -2.0)
                       * 4.0 * M_PI * M_PI / GM_C3 * Mt * Mt, 1.0/3.0) / sin(i2);
    
    *m_in=Mt - Mstar;
    double Mwd  = pow(pow(x1_lt, 3.0) * pow(Porbi_day * SECDAY, -2.0) 
                       * 4.0 * M_PI * M_PI * *m_in* *m_in / GM_C3, 1.0/3.0) / sin(i1);
    *mc1 = Mwd;

    double a1_lt = x1_lt / sin(i1);   /* 内轨道半长轴（光秒） */

    /* ----- 1. 横向多普勒瞬时修正（只算 tt0 时刻） ----- */
    double theta2_tt0 = get_theta2(psr, p, tt0, com1, com2);
    double S2, theta_dot, r_ex, r_ey, r_ez, t_ex, t_ey, t_ez, vr, vt;
    compute_outer_state_and_vel(i1, i2, dOmega, e2, omega2, theta2_tt0,
                                Porbo_day, Mt,
                                &S2, &theta_dot,
                                &r_ex, &r_ey, &r_ez,
                                &t_ex, &t_ey, &t_ez,
                                &vr, &vt);

    double dx_T, dPhi_T;
    transverse_doppler_from_state(i1, Mstar, Mt, a1_lt,
                                  r_ex, r_ey, t_ex, t_ey, vr, vt,
                                  &dx_T, &dPhi_T);
    *x      += dx_T;
    *orbits += dPhi_T / (2.0 * M_PI);


    double dx_incl,dorb_quad;
    
    interp1_inner(tt0,&dx_incl,&dorb_quad);
    
    *x+=dx_incl;
    *orbits+=dorb_quad;

    /* ----- 2. 短周期振荡项（用于 Romer 延迟修正） ----- */
    /* 计算内双星相对半长轴 a_in (lt-s) */
    double a_in_lt = a1_lt * (*m_in) / Mwd;

    /* 无量纲比例因子 ratio = (a_in / S2)^3 * (Mstar / Min) */
    double ratio = pow(a_in_lt / S2, 3.0) * Mstar / *m_in;

    /* 内轨道平经度 phi：从 orbits 得到相位 */
    double phi = 2.0 * M_PI * fmod(*orbits, 1.0);
    if (phi < 0.0) phi += 2.0 * M_PI;

    double x2 = r_ex * r_ex;
    double y2 = r_ey * r_ey;
    double xy = r_ex * r_ey;

    /* 半长轴振荡 δa_in (lt-s)，对应投影 δx = δa_in * (Mwd/M_in) * sin(i1) */
    double delta_a_in = a_in_lt * ratio * (
        3.0 * xy * sin(2.0 * phi) + 1.5 * (x2 - y2) * cos(2.0 * phi)
    );
    double delta_x = delta_a_in * Mwd / (*m_in) * sin(i1);

    /* 倾角瞬时振荡 δi（弧度）*/
    double delta_i = 0.75 * ratio * r_ez * (
        r_ex * sin(2.0 * phi) - r_ey * cos(2.0 * phi)
    );
    delta_x += a1_lt * cos(i1) * delta_i;

    /* 相位振荡 δphi (来自频率扰动的二阶积分) */
    double cot_i1 = 1.0 / tan(i1);
    double delta_phi = ratio * (
        (21.0/4.0 * xy + 3.0/4.0 * r_ez * r_ex * cot_i1) * cos(2.0 * phi)
        + (-21.0/8.0 * (x2 - y2) + 3.0/4.0 * r_ez * r_ey * cot_i1) * sin(2.0 * phi)
    );

    /* 偏心率矢量振荡 δk1 = δ(e sinω), δk2 = δ(e cosω) */
    double delta_k1 = ratio * (
        (1.0 - 3.75 * x2 + 0.75 * y2) * sin(phi)
        + 4.5 * xy * cos(phi)
        + 0.25 * (x2 - y2) * sin(3.0 * phi)
        - 0.5 * xy * cos(3.0 * phi)
    );
    double delta_k2 = ratio * (
        (1.0 + 0.75 * x2 - 3.75 * y2) * cos(phi)
        + 4.5 * xy * sin(phi)
        + 0.25 * (x2 - y2) * cos(3.0 * phi)
        + 0.5 * xy * sin(3.0 * phi)
    );

    /* Romer 延迟修正： x cosφ δφ + δx sinφ - 0.5 x δk1 (cos2φ+3) + 0.5 x δk2 sin2φ */
    *dRomer = x1_lt * cos(phi) * delta_phi
              + delta_x * sin(phi)
              - 0.5 * x1_lt * delta_k1 * (cos(2.0 * phi) + 3.0)
              + 0.5 * x1_lt * delta_k2 * sin(2.0 * phi);
    
}

double btmodel(pulsar *psr,int p,int ipos,int param,int arr,double *torb,double *vorb,int com){

    double tt0,orbits;
    double pb,pbdot,pb2dot;     /* Orbital period (sec) */
    double ecc,edot;    /* Orbital eccentricity */
    double x,xdot;
    double omz,omdot;
    double gamma;
    int    norbits;
    double phase;
    double ep,dep,bige,tt,sin_omega,cos_omega;
    double alpha,beta,sbe,cbe,q,r,s;
    int ii;

    pbdot=0.0;
    pb2dot=0.0;
    xdot = 0.0;
    omdot = 0.0;
    edot = 0.0;

    tt0 = (psr[p].obsn[ipos].bbat - psr[p].param[param_t0].val[com])*SECDAY;
    pb     = psr[p].param[param_pb].val[com] * SECDAY;
    ecc    = psr[p].param[param_ecc].val[com] + edot*tt0;

    tt0+=torb[0]+torb[1];

    if (ecc < 0.0 || ecc > 1.0)
    {
        ld_printf("BTmodel: problem with eccentricity = %Lg\n",psr[p].param[param_ecc].val[com]);
        exit(1);
    }

    x  = psr[p].param[param_a1].val[com];

    omz  = (psr[p].param[param_om].val[com])/180.0*M_PI;

    if (psr[p].param[param_gamma].paramSet[com]==1) gamma = psr[p].param[param_gamma].val[com];
    else gamma  = 0.0;

    double Mt=psr[p].param[param_mtot].val[0];
    double i2=psr[p].param[param_b0].val[0]*M_PI/180;
    double Mstar= pow(pow(x,3.0)*pow(pb,-2.0)*4*M_PI*M_PI/GM_C3*Mt*Mt,1.0/3)/sin(i2);


    gamma = ecc*pow(pb/(2*M_PI),1.0/3.0)*pow(GM_C3,2.0/3.0)*(Mstar+Mt)*Mstar/pow(Mt,4.0/3);

    /* Should ct be the barycentric arrival time? -- check bnrybt.f */
    orbits = tt0/pb - 0.5*(pbdot)*pow(tt0/pb,2.0) - 1./6.*pb2dot*pow(tt0/pb,3.0);

    double x_tab, e_tab, om_tab, ph_tab;
    interp1_outer(tt0, &x_tab, &e_tab,& om_tab, &ph_tab);

    x     += x_tab;
    ecc   += e_tab;
    omz   += om_tab;
    orbits += ph_tab;

    norbits = (int)orbits;
    if (orbits < 0.0) norbits--;

    phase = 2.0*M_PI * (orbits-norbits);

    /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
    ep = phase + ecc*sin(phase)*(1.0+ecc*cos(phase));

    /* This line is wrong in the original tempo: should be inside the do loop */
    /*  denom = 1.0 - ecc*cos(ep);*/

    do {
        dep = (phase - (ep-ecc*sin(ep)))/(1.0 - ecc*cos(ep));
        ep += dep;
    } while (fabs(dep) > 1.0e-12);
    bige = ep;

    tt = 1.0-ecc*ecc;
    sin_omega = sin(omz);
    cos_omega = cos(omz);

    alpha = x*sin_omega;
    beta = x*cos_omega*sqrt(tt);
    sbe = sin(bige);
    cbe = cos(bige);
    q = alpha * (cbe-ecc) + (beta+gamma)*sbe;
    r = -alpha*sbe + beta*cbe;
    s = 1.0/(1.0-ecc*cbe);

    torb[com] = -q;
    vorb[com] = (2*M_PI/pb)*(r*s+gamma*cbe*s);

    if (param!=-1 && com==arr){
        if (param==param_pb)
            return -2.0*M_PI*r*s/pb*tt0/pb;  /* fctn(12+j) */
        else if (param==param_a1)
            return (sin_omega*(cbe-ecc) + cos_omega*sbe*sqrt(tt));                /* fctn(9+j) */
        else if (param==param_ecc)
            return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt; /* fctn(10+j) */
        else if (param==param_om)
            return x*(cos_omega*(cbe-ecc) - sin_omega*sqrt(tt)*sbe);          /* fctn(13+j) */
        else if (param==param_t0)
            return -2.0*M_PI/pb*r*s*SECDAY;                           /* fctn(11+j) */
        else if (param==param_gamma) 
            return sbe;
    }
    
    return 0.0;
            
}

double ell1modl(pulsar *psr,int p,int ipos,int param,int arr,double *torb,double*vorb,int com){
    double an,x0,m2,tt0,orbits,phase,e1,e2,dre,drep,drepp,brace,dlogbr,ds,da,pb;
    double eps1,eps2,eps1dot,eps2dot,si,a0,b0;
    double d2bar,Csigma,Cx,Ceps1,Ceps2,Cm2,Csi,ct,t0asc,pbdot,pb2dot,xpbdot,x,xdot,x2dot,am2;
    int norbits;
    int ii;
    double SUNMASS = 4.925490947e-6;
    const char *CVS_verNum = "$Id: 853536c7e64efde12efd9a9b936741f9eb9bec4a $";

    if (displayCVSversion == 1) CVSdisplayVersion("ELL1model.C","ELL1model()",CVS_verNum);

    a0 = 0.0; /* WHAT SHOULD THESE BE? */
    b0 = 0.0;

    if (psr[p].param[param_fb].paramSet[0]==1) 
        pb    = 1.0/psr[p].param[param_fb].val[0];
    else
        pb    = psr[p].param[param_pb].val[com]*SECDAY;

    if (psr[p].param[param_pbdot].paramSet[com] == 1) pbdot = psr[p].param[param_pbdot].val[com];
    else pbdot=0.0;

    if (psr[p].param[param_pb2dot].paramSet[com] == 1) pb2dot  = psr[p].param[param_pb2dot].val[com];
    else pb2dot=0.0;

    an    = 2.0*M_PI/pb;

    si=sin(psr[p].param[param_a0].val[com]*M_PI/180);

    x0    = psr[p].param[param_a1].val[com];
    if (psr[p].param[param_a1dot].paramSet[com] == 1) 
        xdot  = psr[p].param[param_a1dot].val[com];
    else 
        xdot = 0.0;

    if (psr[p].param[param_a2dot].paramSet[0] == 1) 
        x2dot  = psr[p].param[param_a2dot].val[0];
    else 
        x2dot = 0.0;

    t0asc = psr[p].param[param_tasc].val[com];

    if (psr[p].param[param_m2].paramSet[com]==1) am2 = psr[p].param[param_m2].val[com];
    else am2 = 0.0;
    // m2 = am2*SUNMASS;


    eps1  = psr[p].param[param_eps1].val[com];
    eps2  = psr[p].param[param_eps2].val[com];
    if (psr[p].param[param_eps1dot].paramSet[com]==1) eps1dot = psr[p].param[param_eps1dot].val[com];
    else eps1dot=0;
    if (psr[p].param[param_eps2dot].paramSet[com]==1) eps2dot = psr[p].param[param_eps2dot].val[com];
    else eps2dot=0;

    ct = psr[p].obsn[ipos].bbat;

    tt0 = (ct-t0asc)*SECDAY;


    tt0+=torb[0]+torb[1];


    // --- Changes to handle higher orbital-frequency derivatives (FB1, FB2, ...) ---                                               
    // 04/2015, H. J. Pletsch                                                                                                       
    orbits = tt0/pb;

    orbits -= (0.5*pbdot*pow(tt0/pb,2) + 1.0/6*pb2dot*pow(tt0/pb,3));

    x = x0+xdot*tt0+0.5*x2dot*tt0*tt0;


    double dRomer,Mc1,M_in;
    get_correction(psr,p,tt0,0,1,&orbits,&x,&dRomer,&Mc1,&M_in);
    double omdot=3*pow(2*M_PI/pb,5.0/3)*pow(GM_C3*M_in,2.0/3);
    eps1dot=eps2*omdot;
    eps2dot=-eps1*omdot;


    m2=Mc1*SUNMASS;

    norbits = (int)orbits;
    if (orbits<0.0) norbits = norbits-1;
    phase = 2.0*M_PI*(orbits-norbits);

    e1 = eps1+eps1dot*tt0;
    e2 = eps2+eps2dot*tt0;
  
    long double ecc,omega,u,du,su,cu,onemecu,cae,sae,ae,sw,cw,alpha,beta,bg,anhat,sqr1me2,cume;
    longdouble Comega,Ce;
    
    ecc=sqrt(e1*e1+e2*e2);
    if(ecc!=0){
        omega=atan2(e1,e2);
        if(omega<0) omega+=2*M_PI;
        phase=phase-omega; 
        if(phase<0) phase+=2*M_PI;
    }
    else{
        omega=0;
    }

    x=x*(1+vorb[1]);

    /*  Compute eccentric anomaly u by iterating Kepler's equation. */
    u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
    do {
        du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
        u=u+du;
    } while (fabs(du)>1.0e-12);

    /*  DD equations 17b, 17c, 29, and 46 through 52 */
    
    su=sin(u);
    cu=cos(u);
    onemecu=1.0-ecc*cu;
    sw=sin(omega);
    cw=cos(omega);
    alpha=x*sw;
    beta=x*sqrt(1-pow(ecc,2))*cw;
    bg=beta;
    dre=alpha*(cu-ecc) + bg*su;
    drep=-alpha*su + bg*cu;
    drepp=-alpha*cu - bg*su;
    anhat=an/onemecu;

    /* DD equations 26, 27, 57: */
    sqr1me2=sqrt(1-pow(ecc,2));
    cume=cu-ecc;
    brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
    // printf("GEORGE: si = %g, brace = %g\n",(double)si,(double)brace);
    dlogbr=log(brace);
    ds=-2*m2*dlogbr;

    d2bar=dre + ds + dRomer;
    torb[com]=-d2bar;
    vorb[com]=anhat*drep;
    phase=phase+omega;
    
    Comega   = x*(cw*cume-sqr1me2*sw*su);
    Csigma   = x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    Ce       = su*Csigma-x*sw-ecc*x*cw*su/sqr1me2;
    Cx       = sw*(cu-ecc)+sqr1me2*cw*su;
    if(ecc==0){
        Ceps1    = -0.5*x*cos(2*phase);
        Ceps2    =  0.5*x*sin(2*phase);
    }
    else{
        Ceps1    = Ce*sin(omega)+(Comega-Csigma)*cos(omega)/ecc;
        Ceps2    = Ce*cos(omega)-(Comega-Csigma)*sin(omega)/ecc;
    }
    Cm2=-2*dlogbr;
    Csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace;

    

    Cx*=(1+vorb[1]);
    /* Now we need the partial derivatives. */
    if (param!=-1 && com==arr){
        if (param==param_pb)
            return -Csigma*an*tt0/pb; /* Pb    */
        else if (param==param_a1)
            return Cx;
        else if (param==param_eps1)
            return Ceps1;
        else if (param==param_tasc)
            return -Csigma*an*SECDAY;
        else if (param==param_eps2)
            return Ceps2;
        else if (param==param_eps1dot)
            return Ceps1*tt0;
        else if (param==param_eps2dot)
            return Ceps2*tt0;
        else if (param==param_pbdot)
            return 0.5*tt0*(-Csigma*an*tt0/pb);
        else if (param==param_a1dot)
            return Cx*tt0;      
        else if( param == param_a2dot)
            return 0.5*Cx*tt0*tt0;  
        else if (param==param_sini)
            return Csi;
        else if (param==param_m2)
            return Cm2*SUNMASS;
        else if (param==param_pb2dot)
            return  1./6.*tt0*(-Csigma*an*tt0/pb)*tt0/pb;
    }

    return 0.0;
}

double T2model(pulsar *psr,int p,int ipos,int param,int arr)
{
   

    const char *CVS_verNum = "$Id: 34a0d83c9c1e7a2002ac8e42218b95b3f729401e $";

    if (displayCVSversion == 1) 
        CVSdisplayVersion("T2model.C","T2model()",CVS_verNum);

    int com,com1,com2;

    if (param==-1) 
    {
        com1 = 0;
        com2 = psr[p].nCompanion;
    }
    else
    {
        com1 = arr;
        com2 = arr+1;
    }
    double torb[psr[p].nCompanion];
    double vorb[psr[p].nCompanion];
    for(com=0;com<psr[p].nCompanion;com++){
        torb[com]=0;
        vorb[com]=0;
    }

    if(ipos==0||inner_builded==0){
        build_inner_corr_table(psr,p,0,1);
        inner_builded=1;
    }

    if(ipos==0||outer_builded==0){
        build_outer_corr_table(psr,p,0,1);
        outer_builded=1;
    }
    
    int ii;
    double cparam=0;

    for (com = 1; com >= 0;com--){
        if(psr[p].param[param_t0].paramSet[com]==1){
            cparam=btmodel(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else if(psr[p].param[param_tasc].paramSet[com]==1){
            cparam=ell1modl(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else{
            printf("TO or TASC is required!!\n");
        }
    }

    for (com = 1; com >= 0;com--){
        if(psr[p].param[param_t0].paramSet[com]==1){
            cparam=btmodel(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else if(psr[p].param[param_tasc].paramSet[com]==1){
            cparam=ell1modl(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else{
            printf("TO or TASC is required!!\n");
        }
    }

    for (com = com2-1; com >= com1;com--){
        if(psr[p].param[param_t0].paramSet[com]==1){
            cparam=btmodel(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else if(psr[p].param[param_tasc].paramSet[com]==1){
            cparam=ell1modl(psr,p,ipos,param,arr,torb,vorb,com);
        }
        else{
            printf("TO or TASC is required!!\n");
        }
        if (param==-1 && com == com1){    
            return get_sum(psr[p].nCompanion,torb);
        }
        else if (param!=-1 && com==arr){
            return cparam*(1-vorb[0]-vorb[1]);
        }

    }
    
    return 0.0;
}

void updateT2(pulsar *psr,double val,double err,int pos,int arr)
{
    if (pos==param_pb)
    {
        psr->param[param_pb].val[arr] += val/SECDAY;
        psr->param[param_pb].err[arr]  = err/SECDAY;
    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_gamma || pos==param_edot 
                                             || pos==param_eps1 || pos==param_eps2 || pos==param_tasc )
    {
        psr->param[pos].val[arr] += val;
        psr->param[pos].err[arr]  = err;
    }
    else if (pos==param_om)
    {
        psr->param[pos].val[arr] += val*180.0/M_PI;
        psr->param[pos].err[arr]  = err*180.0/M_PI;
    }
    else if (pos==param_pbdot || pos==param_pb2dot || pos==param_eps1dot || pos==param_eps2dot)
    {
        psr->param[pos].val[arr] += val;
        psr->param[pos].err[arr]  = err;
    }
    else if (pos==param_omdot)
    {
        psr->param[pos].val[arr] += val*(SECDAY*365.25)*180.0/M_PI;
        psr->param[pos].err[arr]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
    else if (pos==param_a1dot||pos==param_a2dot)
    {
        psr->param[pos].val[arr] += val;
        psr->param[pos].err[arr]  = err;
    }
}

/* 线性插值 */
static void interp1_outer(double t, double *x2, double *e2, double *om2, double *dphase)
{
    if (outer_nodes == 0) {
        *x2 = *e2 = *om2 = *dphase = 0.0;
        return;
    }
    if (outer_nodes == 1 || t <= outer_corr_table[0].time) {
        *x2     = outer_corr_table[0].x2;
        *e2     = outer_corr_table[0].e2;
        *om2    = outer_corr_table[0].om2;
        *dphase = outer_corr_table[0].dphase;
        return;
    }
    int i = 0, j = outer_nodes - 1;
    while (i < j) {
        int mid = (i + j) / 2;
        if (outer_corr_table[mid].time < t) i = mid + 1;
        else j = mid;
    }
    if (i == 0 || i >= outer_nodes) {
        *x2     = outer_corr_table[i].x2;
        *e2     = outer_corr_table[i].e2;
        *om2    = outer_corr_table[i].om2;
        *dphase = outer_corr_table[i].dphase;
        return;
    }
    double frac = (t - outer_corr_table[i-1].time)
                  / (outer_corr_table[i].time - outer_corr_table[i-1].time);
    *x2     = outer_corr_table[i-1].x2     + frac * (outer_corr_table[i].x2     - outer_corr_table[i-1].x2);
    *e2     = outer_corr_table[i-1].e2     + frac * (outer_corr_table[i].e2     - outer_corr_table[i-1].e2);
    *om2    = outer_corr_table[i-1].om2    + frac * (outer_corr_table[i].om2    - outer_corr_table[i-1].om2);
    *dphase = outer_corr_table[i-1].dphase + frac * (outer_corr_table[i].dphase - outer_corr_table[i-1].dphase);
}

/* ---------- 由平近点角计算真近点角 ---------- */
static double true_anomaly(double M, double e) {
    double E = M + e * sin(M);
    for (int i = 0; i < 10; ++i) {
        double dE = (M - E + e * sin(E)) / (1.0 - e * cos(E));
        E += dE;
        if (fabs(dE) < 1e-12) break;
    }
    double cos_th = (cos(E) - e) / (1.0 - e * cos(E));
    double sin_th = sqrt(1.0 - e*e) * sin(E) / (1.0 - e * cos(E));
    return atan2(sin_th, cos_th);
}

/* ---------- 计算四极矩摄动函数的几何因子 ---------- */
static double eval_outer_R(double i1_rad, double dOmega_rad,
                           double M2, double a2_lt, double e2,
                           double i2_rad, double omega2_rad)
{
    double theta2 = true_anomaly(M2, e2);
    double u = omega2_rad + theta2;
    double S2 = a2_lt * (1.0 - e2*e2) / (1.0 + e2 * cos(theta2));
    double r_ex, r_ey, r_ez, t_ex, t_ey, t_ez;
    compute_unit_vectors(i1_rad, i2_rad, dOmega_rad, u,
                         &r_ex, &r_ey, &r_ez, &t_ex, &t_ey, &t_ez);
    double Q = 3.0 * r_ez * r_ez - 1.0;
    return -Q / (4.0 * S2 * S2 * S2);
}

/* ---------- 外轨道根数摄动变化率 ---------- */
static void compute_outer_perturbation_rates(
    double i1_rad, double dOmega_rad,
    double M2, double a2_lt, double e2,
    double i2_rad, double omega2_rad, double n2_0,
    double Mt_Msun, double Mstar_Msun, double Mwd_Msun,
    double a1_rel_lt,   // 内双星相对半长轴 (光秒)
    double *a2_dot, double *e2_dot, double *i2_dot,
    double *Omega2_dot, double *omega2_dot, double *deltaM2_dot)
{
    double Min = Mt_Msun - Mstar_Msun;
    double Mpsr = Min - Mwd_Msun;
    double n2 = sqrt(GM_C3 * Mt_Msun / (a2_lt * a2_lt * a2_lt));

    // 差分步长
    double da = a2_lt * 1e-4;
    double de = 1e-4;
    double di = 1e-4;
    double dOm = 1e-4;
    double dom = 1e-4;
    double dM = 1e-4;

    // 四极矩耦合系数
    double ratio = GM_C3 * (Mpsr * Mwd_Msun * Mt_Msun) / (Min * Min) * a1_rel_lt * a1_rel_lt;

    // 中心差分计算偏导数
    double dR_da = (eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt+da, e2, i2_rad, omega2_rad) -
                    eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt-da, e2, i2_rad, omega2_rad)) / (2.0 * da);
    double dR_de = (eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2+de, i2_rad, omega2_rad) -
                    eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2-de, i2_rad, omega2_rad)) / (2.0 * de);
    double dR_di = (eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2, i2_rad+di, omega2_rad) -
                    eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2, i2_rad-di, omega2_rad)) / (2.0 * di);
    double dR_dOm= (eval_outer_R(i1_rad, dOmega_rad+dOm, M2, a2_lt, e2, i2_rad, omega2_rad) -
                    eval_outer_R(i1_rad, dOmega_rad-dOm, M2, a2_lt, e2, i2_rad, omega2_rad)) / (2.0 * dOm);
    double dR_dom= (eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2, i2_rad, omega2_rad+dom) -
                    eval_outer_R(i1_rad, dOmega_rad, M2, a2_lt, e2, i2_rad, omega2_rad-dom)) / (2.0 * dom);
    double dR_dM = (eval_outer_R(i1_rad, dOmega_rad, M2+dM, a2_lt, e2, i2_rad, omega2_rad) -
                    eval_outer_R(i1_rad, dOmega_rad, M2-dM, a2_lt, e2, i2_rad, omega2_rad)) / (2.0 * dM);

    double sin_i2 = sin(i2_rad);
    double cos_i2 = cos(i2_rad);
    double sqrt_1me2 = sqrt(1.0 - e2*e2);
    double n2a2 = n2 * a2_lt * a2_lt;

    // 拉格朗日行星方程 (Murray & Dermott 1999)
    *a2_dot     = ratio * 2.0 / (n2 * a2_lt) * dR_dM;
    *e2_dot     = ratio * ((1.0 - e2*e2) / (n2a2 * e2) * dR_dM
                   - sqrt_1me2 / (n2a2 * e2) * dR_dom);
    *i2_dot     = ratio / (n2a2 * sqrt_1me2 * sin_i2) * (cos_i2 * dR_dom - dR_dOm);
    *Omega2_dot = ratio / (n2a2 * sqrt_1me2 * sin_i2) * dR_di;
    *omega2_dot = ratio * (sqrt_1me2 / (n2a2 * e2) * dR_de
                   - cos_i2 / (n2a2 * sqrt_1me2 * sin_i2) * dR_di);
    *deltaM2_dot     = n2-n2_0 + ratio * ( - (1.0 - e2*e2) / (n2a2 * e2) * dR_de
                   - 2.0 / (n2 * a2_lt) * dR_da );
}

static void build_outer_corr_table(pulsar *psr, int p, int com_inner, int com_outer) {
    double i1  = psr[p].param[param_a0].val[com_inner] * M_PI / 180.0;
    double i2  = psr[p].param[param_b0].val[0] * M_PI / 180.0;
    double dOmega = psr[p].param[param_kom].val[0] * M_PI / 180.0;
    double Porbi_day = psr[p].param[param_pb].val[com_inner];
    double Porbo_day = psr[p].param[param_pb].val[com_outer];
    double T0_2 = psr[p].param[param_t0].val[com_outer];
    double e2_0 = psr[p].param[param_ecc].val[com_outer];
    double omega2_0 = psr[p].param[param_om].val[com_outer] * M_PI / 180.0;
    double x2_0 = psr[p].param[param_a1].val[com_outer];
    double i2_0 = i2;   // 注意倾角使用 B0
    double Mt = psr[p].param[param_mtot].val[0];

    double Mstar = pow(pow(x2_0,3.0) / pow(Porbo_day*SECDAY,2.0) * 4.0*M_PI*M_PI/GM_C3 * Mt*Mt, 1.0/3.0) / sin(i2_0);
    double Min = Mt - Mstar;
    double Mwd = pow(pow(Min,2.0) * pow(psr[p].param[param_a1].val[com_inner],3.0) / pow(Porbi_day*SECDAY,2.0) * 4.0*M_PI*M_PI/GM_C3, 1.0/3.0) / sin(i1);
    double a1_rel_lt = (psr[p].param[param_a1].val[com_inner] / sin(i1)) * Min / Mwd;
    double a2_0 = pow(GM_C3 * Mt * pow(Porbo_day*SECDAY,2.0) / (4.0*M_PI*M_PI), 1.0/3.0);
    double n2 = 2.0*M_PI / (Porbo_day * SECDAY);

    double t_start = (psr[p].param[param_start].val[0] - T0_2 - 500) * SECDAY;
    double t_finish = (psr[p].param[param_finish].val[0] - T0_2 + 500) * SECDAY;

    double step = Porbo_day * SECDAY / 512;
    int N_before = (int)((0.0 - t_start) / step) + 1;
    int N_after = (int)((t_finish - 0.0) / step) + 1;
    if(N_before<0) N_before=0;
    if(N_after<0) N_after=0;
    if (N_before+N_after > OUTER_CORR_NODES){
        step=(t_finish-t_start)/(OUTER_CORR_NODES-4);
        N_before = (int)((0.0 - t_start) / step);
        N_after = (int)((t_finish - 0.0) / step) + 1;
    }
    if(N_before<0) N_before=0;
    if(N_after<0) N_after=0;
    outer_nodes = N_before+N_after;


    double t = 0;
    double om2 = omega2_0, e2 = e2_0, a2 = a2_0, i2c = i2_0, M2 = n2*0.5*step, dphase = 0.0;
    double a2_dot, e2_dot, i2_dot, Om2_dot, om2_dot, deltaM2_dot;
    int i;
    for (i = N_before; i < outer_nodes; i++) {
        outer_corr_table[i].time   = t;
        outer_corr_table[i].x2     = a2*Mstar/Mt*sin(i2c);   // 投影半长轴（考虑质量比）
        outer_corr_table[i].e2     = e2;
        outer_corr_table[i].om2    = om2;
        outer_corr_table[i].dphase = dphase;
        compute_outer_perturbation_rates(i1, dOmega, M2, a2, e2, i2c, om2, n2,
                                         Mt, Mstar, Mwd, a1_rel_lt,
                                         &a2_dot, &e2_dot, &i2_dot, &Om2_dot, &om2_dot, &deltaM2_dot);
        a2  += step * a2_dot;
        e2  += step * e2_dot;
        i2c += step * i2_dot;
        om2 += step * om2_dot;
        dphase += step * deltaM2_dot/ (2.0 * M_PI);   // 摄动部分相位
        M2 += step * (n2 +deltaM2_dot);    // 开普勒部分
        t  += step;
    }

    step=-step;
    // double om2 = omega2_0, e2 = e2_0, a2 = a2_0, i2c = i2_0, M2 = n2*(t_start+0.5*step), dphase = 0.0;

    om2 = omega2_0; 
    e2 = e2_0;
    a2 = a2_0;
    i2c = i2_0;
    M2 = n2*0.5*step;
    dphase = 0.0;
    for (i = N_before-1; i >= 0; i--) {
        compute_outer_perturbation_rates(i1, dOmega, M2, a2, e2, i2c, om2, n2,
                                         Mt, Mstar, Mwd, a1_rel_lt,
                                         &a2_dot, &e2_dot, &i2_dot, &Om2_dot, &om2_dot, &deltaM2_dot);
        a2  += step * a2_dot;
        e2  += step * e2_dot;
        i2c += step * i2_dot;
        om2 += step * om2_dot;
        dphase += step * deltaM2_dot/ (2.0 * M_PI);   // 摄动部分相位
        M2 += step * (n2 +deltaM2_dot);   // 开普勒部分
        t  += step;

        outer_corr_table[i].time   = t;
        outer_corr_table[i].x2     = a2*Mstar/Mt*sin(i2c);   // 投影半长轴（考虑质量比）
        outer_corr_table[i].e2     = e2;
        outer_corr_table[i].om2    = om2;
        outer_corr_table[i].dphase = dphase;
    }

    // 插值获取 T0 (t=0) 时刻
    double x2_T0,e2_T0,om2_T0,dphase_T0;  
    interp1_outer(0.0,&x2_T0,&e2_T0,&om2_T0,&dphase_T0);
    

    // 将所有节点的值减去 T0 时刻的值，变成相对于 T0 的偏移量
    for (int i = 0; i < outer_nodes; i++) {
        outer_corr_table[i].x2     -= x2_T0;
        outer_corr_table[i].e2     -= e2_T0;
        outer_corr_table[i].om2    -= om2_T0;
        outer_corr_table[i].dphase -= dphase_T0;
    }

}


/* 线性插值 */
static void interp1_inner(double t, double *dx, double *dorbits) {
    if (inner_nodes == 0) { *dx = *dorbits = 0.0; return; }
    if (inner_nodes == 1 || t <= inner_corr_table[0].time) {
        *dx = inner_corr_table[0].dx; *dorbits = inner_corr_table[0].dorbits; return;
    }
    int i = 0, j = inner_nodes - 1;
    while (i < j) {
        int mid = (i + j) / 2;
        if (inner_corr_table[mid].time < t) i = mid + 1; else j = mid;
    }
    if (i == 0 || i >= inner_nodes) { *dx = inner_corr_table[i].dx; *dorbits = inner_corr_table[i].dorbits; return; }
    double frac = (t - inner_corr_table[i-1].time) / (inner_corr_table[i].time - inner_corr_table[i-1].time);
    *dx      = inner_corr_table[i-1].dx      + frac * (inner_corr_table[i].dx      - inner_corr_table[i-1].dx);
    *dorbits = inner_corr_table[i-1].dorbits + frac * (inner_corr_table[i].dorbits - inner_corr_table[i-1].dorbits);
}

static void build_inner_corr_table(pulsar *psr, int p, int com_inner, int com_outer) {
    double i1 = psr[p].param[param_a0].val[com_inner] * M_PI / 180.0;
    double i2 = psr[p].param[param_b0].val[0] * M_PI / 180.0;
    double dOmega = psr[p].param[param_kom].val[0] * M_PI / 180.0;
    double e2 = psr[p].param[param_ecc].val[com_outer];
    double omega2 = psr[p].param[param_om].val[com_outer] * M_PI / 180.0;
    double Mt = psr[p].param[param_mtot].val[0];
    double Pb_day = psr[p].param[param_pb].val[com_outer];
    double x1_lt = psr[p].param[param_a1].val[com_inner];
    double x2_lt = psr[p].param[param_a1].val[com_outer];
    double Porbi_day = psr[p].param[param_pb].val[com_inner];

    double Mstar = pow(pow(x2_lt,3.0) / pow(Pb_day*SECDAY,2.0) * 4.0*M_PI*M_PI/GM_C3 * Mt*Mt, 1.0/3.0) / sin(i2);
    double Min = Mt - Mstar;
    double Mwd = pow(pow(Min,2.0) * pow(x1_lt,3.0) / pow(Porbi_day*SECDAY,2.0) * 4.0*M_PI*M_PI/GM_C3, 1.0/3.0) / sin(i1);
    double a1_lt = x1_lt / sin(i1);

    double tasc = psr[p].param[param_tasc].val[com_inner];
    double t_start = (psr[p].param[param_start].val[0] - tasc - 500) * SECDAY;
    double t_finish = (psr[p].param[param_finish].val[0] - tasc + 500) * SECDAY;

    double step = Porbi_day * SECDAY / 2048;
    int N = (int)((t_finish - t_start) / step) + 1;
    if (N > INNER_CORR_NODES) N = INNER_CORR_NODES;
    step = (t_finish - t_start) / N;
    inner_nodes = N;

    double accum_dx = 0.0, accum_dorb = 0.0;
    double t = t_start;

    for (int i = 0; i < N; i++) {
        inner_corr_table[i].time = t;
        inner_corr_table[i].dx = accum_dx;
        inner_corr_table[i].dorbits = accum_dorb;

        double tt = t + 0.5 * step;
        double theta2 = get_theta2(psr, p, tt, com_inner, com_outer);
        double S2, theta_dot, r_ex, r_ey, r_ez, t_ex, t_ey, t_ez, vr, vt;
        compute_outer_state_and_vel(i1, i2, dOmega, e2, omega2, theta2,
                                    Pb_day, Mt,
                                    &S2, &theta_dot,
                                    &r_ex, &r_ey, &r_ez,
                                    &t_ex, &t_ey, &t_ez,
                                    &vr, &vt);

        double dx_dt = inclination_xdot_from_state(i1, Mstar, Mwd, Mt, Porbi_day,
                                                   a1_lt, S2, r_ex, r_ey, r_ez);
        double forb = quadrupole_forb_from_state(i1, Mstar, Mwd, Mt, Porbi_day,
                                                 a1_lt, S2, r_ex, r_ey, r_ez);

        accum_dx += dx_dt * step;
        accum_dorb += forb * step;
        t += step;
    }

    // 使用插值获取 TASC (t=0) 时刻的累积值作为零点偏移
    double dx0, dorb0;
    interp1_inner(0.0, &dx0, &dorb0);

    // 转换为相对于 TASC 的偏移量
    for (int i = 0; i < N; i++) {
        inner_corr_table[i].dx      -= dx0;
        inner_corr_table[i].dorbits -= dorb0;
    }
}