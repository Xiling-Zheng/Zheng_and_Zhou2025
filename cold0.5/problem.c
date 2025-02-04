/**
 * Restricted three body problem.
 *
 * This example simulates a disk of test particles around 
 * a central object, being perturbed by a planet.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#define min(i, j) (((i) < (j)) ? (i) : (j))


double data[51];
/*
0: number
1: initial a
2-10: a
11-19: e
20-28: i
29-37: category
38: qmin
39: qpPmin
40: af
41: ef
42: incf
43: xf
44: yf
45: zf
46: vxf
47: vyf
48: vzf
49: Pxf
50: Pyf
*/
double distancem;
double distancempP;
double distance;
double distancepP;
int index = 2;
int collision_index(struct reb_simulation* const r, struct reb_collision c){
    index = min(c.p1, c.p2);
    if (c.p1<c.p2) {return 2;}
    else {return 1;}
}
void heartbeat(struct reb_simulation* r);
int main(int argc, char *argv[]){
    double a0_min = atof(argv[1]);
    double start, stop;
    start = clock();
    struct reb_simulation* r = reb_create_simulation();
    double mu = 0.0019;
    double e0 = 0.02;
    double inc0 = 0.01;
    char file_name[15];
    sprintf(file_name, "data%s.txt", argv[1]);
    int random_seed = atoi(argv[2]);
    srand(random_seed); // Set random seed
    FILE *f;

    // Setup constants
    // r->dt           = 1.0e-2*2.*M_PI;
    // r->integrator   = REB_INTEGRATOR_MERCURIUS;
    // r->ri_mercurius.hillfac = 14.;
    r->integrator   = REB_INTEGRATOR_IAS15;
    // r->integrator   = REB_INTEGRATOR_WHFAST;
    // r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    // r->boundary     = REB_BOUNDARY_OPEN;
    // r->ri_ias15.min_dt           = 1.0e-4*2.*M_PI;
    // r->ri_ias15.epsilon           = 1.0e-9;
    r->N_active     = 2;    // Only the star and the planet have non-zero mass
    r->heartbeat    = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = collision_index;

    int i = 1;
    while(i<=atoi(argv[3])){
        // Initial conditions for star
        reb_remove_all(r);
        r->t = 0;
        struct reb_particle star = {0};
        star.m = 1 - mu;
        star.r = 0.000895;
        reb_add(r, star);

        // Initial conditions for planet
        double planet_e = 0.2;
        struct reb_particle planet = {0};
        planet.x  = 1.-planet_e;
        planet.vy = sqrt(2./(1.-planet_e)-1.);
        planet.m  = mu;
        planet.r  = 0.00018;
        reb_add(r, planet);
        reb_move_to_com(r);

        struct reb_particle planetesimal = {0};
        double a0 = a0_min + ((double)rand() / (double)RAND_MAX) * 0.05;
        double Omega = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        double omega = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        double M = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI;
        // Initial conditions for planetesimal
        reb_add_fmt(r, "a e inc Omega omega M", a0, e0, inc0, Omega, omega, M);
        data[0] = i;
        data[1] = a0;
        struct reb_orbit o;
        o = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
        struct reb_particle op = r->particles[2];
        struct reb_particle oP = r->particles[1];
        struct reb_particle oS = r->particles[0];
        distancem = sqrt((op.x - oS.x) * (op.x - oS.x) + (op.y - oS.y) * (op.y - oS.y) + (op.z - oS.z) * (op.z - oS.z));
        distancempP = sqrt((op.x - oP.x) * (op.x - oP.x) + (op.y - oP.y) * (op.y - oP.y) + (op.z - oP.z) * (op.z - oP.z));
        data[40] = o.a;
        data[41] = o.e;
        data[42] = o.inc;
        data[43] = op.x;
        data[44] = op.y;
        data[45] = op.z;
        data[46] = op.vx;
        data[47] = op.vy;
        data[48] = op.vz;
        data[49] = oP.x;
        data[50] = oP.y;
        reb_integrate(r, 2.001e6*M_PI);// integrated period

        f = fopen(file_name,"a+");
        int j = 0;
        while(j<=50){
            fprintf(f,"%e\t",data[j]);
            j++;
        }
        fprintf(f,"\n");
        fclose(f);
        i++;
    }

    stop = clock();
    f = fopen("time_used.txt","a+");
    fprintf(f,"%e\t%e\n",a0_min,(stop-start)/CLOCKS_PER_SEC);
    fclose(f);
    return 0;
}

void heartbeat(struct reb_simulation* r){
    struct reb_orbit o;
    struct reb_orbit opP;
    double q;
    double qpP;
    // FILE *ft;
    if (reb_output_check(r, r->dt)){
        if (r->N == 3) {
            o = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
            opP = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[1]);
            struct reb_particle op = r->particles[2];
            struct reb_particle oP = r->particles[1];
            struct reb_particle oS = r->particles[0];
            distance = sqrt((op.x - oS.x) * (op.x - oS.x) + (op.y - oS.y) * (op.y - oS.y) + (op.z - oS.z) * (op.z - oS.z));
            distancepP = sqrt((op.x - oP.x) * (op.x - oP.x) + (op.y - oP.y) * (op.y - oP.y) + (op.z - oP.z) * (op.z - oP.z));
            q = o.a * (1. - o.e);
            qpP = opP.a * (1. - opP.e);
            if (distance > 0.05) {
                if (distance < distancem) {
                    distancem = distance;}
            }
            else {
                if (q < distancem) {
                    distancem = q;}
            }

            if (distancepP > 0.005) {
                if (distancepP < distancempP) {
                    distancempP = distancepP;}
            }
            else {
                if (qpP < distancempP) {
                    distancempP = qpP;}
            }
            data[38] = distancem;
            data[39] = distancempP;
            if (distance > 1.0e2 && o.e > 1) {
                data[40] = o.a;
                data[41] = o.e;
                data[42] = o.inc;
                data[43] = op.x;
                data[44] = op.y;
                data[45] = op.z;
                data[46] = op.vx;
                data[47] = op.vy;
                data[48] = op.vz;
                data[49] = oP.x;
                data[50] = oP.y;   
                reb_remove(r, 2, 1);}
        }
    }
    if (reb_output_check(r, 1.0e2*2.*M_PI)){
        reb_move_to_com(r);
        o = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
        struct reb_particle op = r->particles[2];
        struct reb_particle oP = r->particles[1];
        data[40] = o.a;
        data[41] = o.e;
        data[42] = o.inc;
        data[43] = op.x;
        data[44] = op.y;
        data[45] = op.z;
        data[46] = op.vx;
        data[47] = op.vy;
        data[48] = op.vz;
        data[49] = oP.x;
        data[50] = oP.y;
        int crit1 = 0;
        int w1 = round(r->t / (1.0e2*2.*M_PI));
        int x1;
        for (int wx1 = 0; wx1 <= 4; wx1++) {
            if (w1 == pow(10,wx1)) {crit1 = 1; x1 = (wx1 + 2) * 2;}
            else if (w1 > pow(10,wx1)) {x1 = (wx1 + 2) * 2;}
        }
        if (r->N == 3) {
            if (crit1 == 1) {
            o = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
            data[x1-2] = o.a;
            data[x1+7] = o.e;
            data[x1+16] = o.inc;
            if (o.e < 1.0) { data[x1+25] = (int)1; }
            else { data[x1+25] = (int)2; }
            }
        }
        else if (r->N == 2) {
            if (index == 0) {
                for (int k = x1 + 1; k < 13;) {
                    data[k-2] = 0.;
                    data[k+7] = 0.;
                    data[k+16] = 0.;
                    data[k+25] = (int)3;
                    k = k + 1;
                }
                index = 2;
                reb_remove_all(r);
            }
            else if (index == 1) {
                for (int k = x1 + 1; k < 13;) {
                    data[k-2] = 0.;
                    data[k+7] = 0.;
                    data[k+16] = 0.;
                    data[k+25] = (int)4;
                    k = k + 1;
                    index = 2;
                }
                reb_remove_all(r);
            }
            else {
                for (int k = x1 + 1; k < 13;) {
                    data[k-2] = data[40];
                    data[k+7] = data[41];
                    data[k+16] = data[42];
                    data[k+25] = (int)2;                        
                    k = k + 1;
                }
                reb_remove_all(r);
            }
        }
    }
    if (reb_output_check(r, pow(10,2.5)*2.*M_PI)){
        int crit2 = 0;
        int w2 = round(r->t / (pow(10,2.5)*2.*M_PI));
        int x2;
        for (int wx2 = 0; wx2 <= 3; wx2++) {
            if (w2 == pow(10,wx2)) {crit2 = 1; x2 = (wx2 + 2) * 2 + 1;}
            if (w2 > pow(10,wx2)) {x2 = (wx2 + 2) * 2 + 1;}
        }
        if (r->N == 3) {
            if (crit2 == 1) {
            o = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
            if (distance > 1.0e2 && o.e > 1) {reb_remove(r, 2, 1);}
            data[x2-2] = o.a;
            data[x2+7] = o.e;
            data[x2+16] = o.inc;
            if (o.e < 1.0) { data[x2+25] = (int)1; }
            else { data[x2+25] = (int)2; }
            }
        }
        else if (r->N == 2)  {
            if (index == 0) {
                for (int k = x2 + 1; k < 13;) {
                    data[k-2] = 0.;
                    data[k+7] = 0.;
                    data[k+16] = 0.;
                    data[k+25] = (int)3;
                    k = k + 1;
                }
                index = 2;
                reb_remove_all(r);
            }
            else if (index == 1) {
                for (int k = x2 + 1; k < 13;) {
                    data[k-2] = 0.;
                    data[k+7] = 0.;
                    data[k+16] = 0.;
                    data[k+25] = (int)4;
                    k = k + 1;
                    index = 2;
                }
                reb_remove_all(r);
            }
            else {
                for (int k = x2 + 1; k < 13;) {
                    data[k-2] = data[40];
                    data[k+7] = data[41];
                    data[k+16] = data[42];
                    data[k+25] = (int)2;
                    k = k + 1;
                }
                reb_remove_all(r);
            }
        }
    }
}