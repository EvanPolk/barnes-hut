#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits>
#include <ctime>
#include <cstring>
#include <cctype>
#include <chrono>

#ifndef __RESTRICT
#  define __RESTRICT
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define NDIM 2
#define G 1.0
#define TINY (std::numeric_limits<double>::epsilon())

typedef double ValueType;

// Array indexing macros
#define acc_array(i,j,n,acc) acc[(i) + (j)*n]
#define pos_array(i,j,n,pos) pos[(i) + (j)*n]
#define vel_array(i,j,n,vel) vel[(i) + (j)*n]

template<typename T>
struct QuadNode {
    bool is_leaf;
    int body_idx;
    T com_x, com_y;
    T mass;
    T x_min, x_max, y_min, y_max;
    QuadNode<T>* children[4];

    QuadNode(T xmin, T xmax, T ymin, T ymax)
        : is_leaf(true), body_idx(-1), com_x(0), com_y(0), mass(0),
          x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax) {
        for(int i=0;i<4;i++) children[i]=nullptr;
    }

    bool contains(int idx, T* pos, int n) {
        T x = pos_array(idx,0,n,pos);
        T y = pos_array(idx,1,n,pos);
        return x>=x_min && x<x_max && y>=y_min && y<y_max;
    }

    void insert(int idx, T* pos, int n) {
        if(is_leaf) {
            if(body_idx==-1) { body_idx=idx; return; }
            is_leaf=false;
            T xm = 0.5*(x_min+x_max);
            T ym = 0.5*(y_min+y_max);
            children[0]=new QuadNode<T>(x_min,xm,y_min,ym);
            children[1]=new QuadNode<T>(xm,x_max,y_min,ym);
            children[2]=new QuadNode<T>(x_min,xm,ym,y_max);
            children[3]=new QuadNode<T>(xm,x_max,ym,y_max);

            for(int i=0;i<4;i++) if(children[i]->contains(body_idx,pos,n)) { children[i]->insert(body_idx,pos,n); break; }
            for(int i=0;i<4;i++) if(children[i]->contains(idx,pos,n)) { children[i]->insert(idx,pos,n); break; }
            body_idx=-1;
        } else {
            for(int i=0;i<4;i++) if(children[i]->contains(idx,pos,n)) { children[i]->insert(idx,pos,n); break; }
        }
    }

    void compute_com(T* pos, T* mass_array, int n) {
        if(is_leaf) {
            if(body_idx!=-1) {
                com_x = pos_array(body_idx,0,n,pos);
                com_y = pos_array(body_idx,1,n,pos);
                mass = mass_array[body_idx];
            }
        } else {
            mass = 0; com_x = 0; com_y = 0;
            for(int i=0;i<4;i++) {
                if(children[i]) {
                    children[i]->compute_com(pos,mass_array,n);
                    com_x += children[i]->com_x * children[i]->mass;
                    com_y += children[i]->com_y * children[i]->mass;
                    mass += children[i]->mass;
                }
            }
            if(mass>0) { com_x/=mass; com_y/=mass; }
        }
    }

    void compute_acc(int idx, T* pos, T* acc, int n, T theta=0.5) {
        if(is_leaf) {
            if(body_idx==-1 || body_idx==idx) return;
            T dx = pos_array(body_idx,0,n,pos) - pos_array(idx,0,n,pos);
            T dy = pos_array(body_idx,1,n,pos) - pos_array(idx,1,n,pos);
            T dist2 = dx*dx + dy*dy + TINY*TINY;
            T invDist3 = 1.0/(dist2*sqrt(dist2));
            acc_array(idx,0,n,acc)+=G*dx*invDist3;
            acc_array(idx,1,n,acc)+=G*dy*invDist3;
        } else {
            T dx = com_x - pos_array(idx,0,n,pos);
            T dy = com_y - pos_array(idx,1,n,pos);
            T dist = sqrt(dx*dx+dy*dy+TINY*TINY);
            T s = x_max - x_min;
            if(s/dist<theta){
                T invDist3 = 1.0/(dist*dist*dist);
                acc_array(idx,0,n,acc)+=G*mass*dx*invDist3;
                acc_array(idx,1,n,acc)+=G*mass*dy*invDist3;
            } else {
                for(int i=0;i<4;i++) if(children[i]) children[i]->compute_acc(idx,pos,acc,n,theta);
            }
        }
    }
};

template<typename T>
void compute_acc_vec(int n, T* pos, T* acc, QuadNode<T>* root){
    #pragma omp simd
    for(int i=0;i<n;i++){
        root->compute_acc(i,pos,acc,n);
    }
}

template<typename T>
void search(T pos[], T vel[], const int n){
    T minv = 1e10, maxv = 0, ave = 0;
    for(int i=0;i<n;i++){
        T vmag = 0;
        for(int k=0;k<NDIM;k++)
            vmag += vel_array(i,k,n,vel) * vel_array(i,k,n,vel);
        vmag = sqrt(vmag);
        if(vmag > maxv) maxv = vmag;
        if(vmag < minv) minv = vmag;
        ave += vmag;
    }
    ave /= n;
    printf("min/max/ave velocity = %e, %e, %e\n", (double)minv, (double)maxv, (double)ave);
}

void help(const char* prog){
    printf("Usage: %s [options]\n",prog);
    printf("  -n | --nparticles <int>    Number of particles\n");
    printf("  -s | --nsteps <int>        Number of steps\n");
    printf("  -t | --stepsize <float>    Time step size\n");
    printf("  -d | --double               Use double precision\n");
    printf("  -f | --float                Use single precision\n");
    printf("  -h | --help                 Show this help\n");
}

template<typename T>
int run_tests(int n, int num_steps, double dt) {
    T* pos = new T[n*NDIM];
    T* vel = new T[n*NDIM];
    T* acc = new T[n*NDIM];
    T* mass = new T[n];

    srand(n);
    for(int i=0;i<n;i++){
        for(int k=0;k<NDIM;k++){
            pos_array(i,k,n,pos)=2*((T)rand()/RAND_MAX-0.5);
            vel_array(i,k,n,vel)=0.0;
            acc_array(i,k,n,acc)=0.0;
        }
        mass[i]=1.0;
    }

    double t_start, t_end;

#ifdef _OPENMP
    t_start = omp_get_wtime();
#else
    auto t_start_chrono = std::chrono::high_resolution_clock::now();
#endif

    for(int step=0; step<num_steps; step++) {
        QuadNode<T>* root = new QuadNode<T>(-1,1,-1,1);

        // Build tree
        for(int i=0;i<n;i++) root->insert(i,pos,n);
        root->compute_com(pos,mass,n);

        // Compute acceleration
#ifdef OPENMP
#pragma omp parallel for
        for(int i=0;i<n;i++) root->compute_acc(i,pos,acc,n);
#elif defined(SERIAL_VEC)
        compute_acc_vec(n,pos,acc,root);
#else
        for(int i=0;i<n;i++) root->compute_acc(i,pos,acc,n);
#endif

        // Update positions and velocities
        for(int i=0;i<n;i++)
            for(int k=0;k<NDIM;k++){
                pos_array(i,k,n,pos) += vel_array(i,k,n,vel)*dt + 0.5*acc_array(i,k,n,acc)*dt*dt;
                vel_array(i,k,n,vel) += acc_array(i,k,n,acc)*dt;
            }

        // Print stats every 10 steps
        if(step%10==0) search(pos,vel,n);

        delete root;
    }

#ifdef _OPENMP
    t_end = omp_get_wtime();
    double t_total = t_end - t_start;
#else
    auto t_end_chrono = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end_chrono - t_start_chrono;
    double t_total = elapsed.count();
#endif

    printf("Average time per step = %f ms\n", t_total*1000/num_steps);

    delete[] pos; delete[] vel; delete[] acc; delete[] mass;
    return 0;
}

int main(int argc,char* argv[]){
    int n=100, num_steps=100;
    double dt=0.01;
    bool useDouble = true;

    for(int i=1;i<argc;i++){
        #define check_index(i,str) if((i)>=argc){fprintf(stderr,"Missing 2nd argument for %s\n",str); help(argv[0]); return 1;}

        if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0){
            help(argv[0]);
            return 0;
        }
        else if(strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--nparticles")==0){
            check_index(i+1,"--nparticles|-n"); i++;
            if(not isdigit(*argv[i])) { fprintf(stderr,"Invalid value for option \"%s\"\n",argv[i]); help(argv[0]); return 1; }
            n = atoi(argv[i]);
        }
        else if(strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--nsteps")==0){
            check_index(i+1,"--nsteps|-s"); i++;
            if(not isdigit(*argv[i])) { fprintf(stderr,"Invalid value for option \"%s\"\n",argv[i]); help(argv[0]); return 1; }
            num_steps = atoi(argv[i]);
        }
        else if(strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--stepsize")==0){
            check_index(i+1,"--stepsize|-t"); i++;
            if(not isdigit(*argv[i])) { fprintf(stderr,"Invalid value for option \"%s\"\n",argv[i]); help(argv[0]); return 1; }
            dt = atof(argv[i]);
        }
        else if(strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--double")==0){
            useDouble = true;
        }
        else if(strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--float")==0){
            useDouble = false;
        }
        else{
            fprintf(stderr,"Unknown option %s\n",argv[i]);
            help(argv[0]);
            return 1;
        }
    }

    if(useDouble)
        return run_tests<double>(n,num_steps,dt);
    else
        return run_tests<float>(n,num_steps,dt);
}
