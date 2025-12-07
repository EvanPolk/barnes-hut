#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <limits>
#include <ctime>

#ifndef __RESTRICT
#  define __RESTRICT
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define NDIM 2
#define G 1.0
#define TINY (std::numeric_limits<double>::epsilon())
#define acc_array(i,j,n,acc) acc[(i) + (j)*n]
#define pos_array(i,j,n,pos) pos[(i) + (j)*n]
#define vel_array(i,j,n,vel) vel[(i) + (j)*n]

typedef double ValueType;

// Simple QuadNode for demo BH tree
struct QuadNode {
    bool is_leaf;
    int body_idx;
    ValueType com_x, com_y;
    ValueType mass;
    QuadNode* children[4];

    QuadNode() : is_leaf(true), body_idx(-1), com_x(0), com_y(0), mass(0) {
        for(int i=0;i<4;i++) children[i]=nullptr;
    }

    void compute_acc(int idx, ValueType* pos, ValueType* acc, int n) {
        if(is_leaf) {
            if(body_idx == -1 || body_idx == idx) return;
            ValueType dx = pos_array(body_idx,0,n,pos) - pos_array(idx,0,n,pos);
            ValueType dy = pos_array(body_idx,1,n,pos) - pos_array(idx,1,n,pos);
            ValueType dist2 = dx*dx + dy*dy + TINY*TINY;
            ValueType invDist = 1.0 / sqrt(dist2);
            ValueType invDist3 = invDist*invDist*invDist;
            acc_array(idx,0,n,acc) += G*dx*invDist3;
            acc_array(idx,1,n,acc) += G*dy*invDist3;
        } else {
            ValueType dx = com_x - pos_array(idx,0,n,pos);
            ValueType dy = com_y - pos_array(idx,1,n,pos);
            ValueType dist = sqrt(dx*dx + dy*dy + TINY*TINY);
            ValueType s = 1.0; // unit size
            if(s / dist < 0.5) {
                ValueType invDist3 = 1.0 / (dist*dist*dist);
                acc_array(idx,0,n,acc) += G*mass*dx*invDist3;
                acc_array(idx,1,n,acc) += G*mass*dy*invDist3;
            } else {
                for(int i=0;i<4;i++)
                    if(children[i]) children[i]->compute_acc(idx,pos,acc,n);
            }
        }
    }
};

// -------------------- Direct N^2 serial vectorized loop --------------------
void compute_acc_direct(int n, ValueType* pos, ValueType* acc) {
    for(int i=0;i<n;i++){
        ValueType ax=0, ay=0;
        #ifdef SERIAL_VEC
        #pragma omp simd reduction(+:ax,ay)
        #endif
        for(int j=0;j<n;j++){
            if(i==j) continue;
            ValueType dx = pos_array(j,0,n,pos) - pos_array(i,0,n,pos);
            ValueType dy = pos_array(j,1,n,pos) - pos_array(i,1,n,pos);
            ValueType dist2 = dx*dx + dy*dy + TINY*TINY;
            ValueType invDist3 = 1.0 / (dist2*sqrt(dist2));
            ax += dx*invDist3;
            ay += dy*invDist3;
        }
        acc_array(i,0,n,acc) = G*ax;
        acc_array(i,1,n,acc) = G*ay;
    }
}

// -------------------- Velocity statistics --------------------
void search(ValueType pos[], ValueType vel[], const int n) {
    ValueType minv = 1e10, maxv = 0, ave = 0;
    for(int i=0;i<n;i++){
        ValueType vmag = 0;
        for(int k=0;k<NDIM;k++)
            vmag += vel_array(i,k,n,vel)*vel_array(i,k,n,vel);
        vmag = sqrt(vmag);
        if(vmag>maxv) maxv = vmag;
        if(vmag<minv) minv = vmag;
        ave += vmag;
    }
    ave /= n;
    printf("min/max/ave velocity = %e, %e, %e\n", minv, maxv, ave);
}

// -------------------- Main simulation loop --------------------
int main(int argc, char* argv[]) {
    int n = 100;
    int num_steps = 100;
    double dt = 0.01;

    ValueType* pos = new ValueType[n*NDIM];
    ValueType* vel = new ValueType[n*NDIM];
    ValueType* acc = new ValueType[n*NDIM];

    // Random initial positions and zero velocity
    srand(n);
    for(int i=0;i<n;i++){
        for(int k=0;k<NDIM;k++){
            pos_array(i,k,n,pos) = 2*((ValueType)rand()/RAND_MAX-0.5);
            vel_array(i,k,n,vel) = 0.0;
            acc_array(i,k,n,acc) = 0.0;
        }
    }

    QuadNode* root = new QuadNode();
    root->body_idx = 0;
    root->mass = 1.0;
    root->com_x = pos_array(0,0,n,pos);
    root->com_y = pos_array(0,1,n,pos);

    clock_t t_start = clock();
    double t_acc = 0, t_update = 0;

    for(int step=0;step<num_steps;step++){
        clock_t t0 = clock();

        #ifdef OPENMP
        #pragma omp parallel for
        for(int i=0;i<n;i++)
            root->compute_acc(i,pos,acc,n);
        #elif defined(SERIAL_VEC)
        compute_acc_direct(n,pos,acc);
        #else
        for(int i=0;i<n;i++)
            root->compute_acc(i,pos,acc,n);
        #endif

        clock_t t1 = clock();

        // Simple Euler update
        for(int i=0;i<n;i++)
            for(int k=0;k<NDIM;k++){
                pos_array(i,k,n,pos) += vel_array(i,k,n,vel)*dt + 0.5*acc_array(i,k,n,acc)*dt*dt;
                vel_array(i,k,n,vel) += acc_array(i,k,n,acc)*dt;
            }

        clock_t t2 = clock();
        t_acc += (double)(t1-t0)/CLOCKS_PER_SEC;
        t_update += (double)(t2-t1)/CLOCKS_PER_SEC;

        if(step%10==0)
            search(pos,vel,n);
    }

    double t_total = (double)(clock()-t_start)/CLOCKS_PER_SEC;
    printf("Average time per step = %f ms\n", t_total*1000/num_steps);

    delete[] pos;
    delete[] vel;
    delete[] acc;
    delete root;
    return 0;
}
