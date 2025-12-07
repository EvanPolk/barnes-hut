#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <chrono>
#include <cstdlib>

constexpr int NDIM = 2;      // 2D simulation
using ValueType = double;
constexpr ValueType G = 1.0;
constexpr ValueType TINY = std::numeric_limits<ValueType>::epsilon();
constexpr ValueType THETA = 0.5;  // Barnes-Hut opening angle

// -------------------- Fixed Macros --------------------
#define _index(i,j,n) ((i) + (j)*(n))
#define pos_array(i,j,n,pos) pos[_index(i,j,n)]
#define vel_array(i,j,n,vel) vel[_index(i,j,n)]
#define acc_array(i,j,n,acc) acc[_index(i,j,n)]

// -------------------- Quadtree for Barnes-Hut --------------------
struct QuadNode {
    ValueType x_min, x_max, y_min, y_max;
    ValueType mass;
    ValueType com_x, com_y;
    int body_idx;
    bool is_leaf;
    QuadNode* children[4];

    QuadNode(ValueType xmin, ValueType xmax, ValueType ymin, ValueType ymax)
        : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax),
          mass(0.0), com_x(0.0), com_y(0.0), body_idx(-1), is_leaf(true)
    {
        for(int i=0;i<4;i++) children[i]=nullptr;
    }

    ~QuadNode() { for(int i=0;i<4;i++) delete children[i]; }

    int get_quadrant(ValueType x, ValueType y) {
        ValueType xm = 0.5*(x_min + x_max);
        ValueType ym = 0.5*(y_min + y_max);
        if(x < xm) return (y < ym ? 0 : 1); // SW : NW
        else return (y < ym ? 2 : 3);       // SE : NE
    }

    void subdivide() {
        ValueType xm = 0.5*(x_min + x_max);
        ValueType ym = 0.5*(y_min + y_max);
        children[0] = new QuadNode(x_min, xm, y_min, ym);
        children[1] = new QuadNode(x_min, xm, ym, y_max);
        children[2] = new QuadNode(xm, x_max, y_min, ym);
        children[3] = new QuadNode(xm, x_max, ym, y_max);
        is_leaf = false;
    }

    void insert(int idx, ValueType* pos, int n) {
        ValueType x = pos_array(idx,0,n,pos);
        ValueType y = pos_array(idx,1,n,pos);

        if(is_leaf) {
            if(body_idx == -1) body_idx = idx;
            else {
                int old_idx = body_idx;
                subdivide();
                body_idx = -1;
                children[get_quadrant(pos_array(old_idx,0,n,pos), pos_array(old_idx,1,n,pos))]->insert(old_idx,pos,n);
                children[get_quadrant(x,y)]->insert(idx,pos,n);
            }
        } else {
            children[get_quadrant(x,y)]->insert(idx,pos,n);
        }

        ValueType total_mass = mass + 1.0;
        com_x = (com_x*mass + x) / total_mass;
        com_y = (com_y*mass + y) / total_mass;
        mass = total_mass;
    }

    void compute_acc(int idx, ValueType* pos, ValueType* acc, int n) {
        if(is_leaf) {
            if(body_idx == -1 || body_idx == idx) return;
            ValueType dx = pos_array(body_idx,0,n,pos) - pos_array(idx,0,n,pos);
            ValueType dy = pos_array(body_idx,1,n,pos) - pos_array(idx,1,n,pos);
            ValueType dist2 = dx*dx + dy*dy + TINY*TINY;
            ValueType invDist = 1.0 / std::sqrt(dist2);
            ValueType invDist3 = invDist*invDist*invDist;
            acc_array(idx,0,n,acc) += G*dx*invDist3;
            acc_array(idx,1,n,acc) += G*dy*invDist3;
        } else {
            ValueType dx = com_x - pos_array(idx,0,n,pos);
            ValueType dy = com_y - pos_array(idx,1,n,pos);
            ValueType dist = std::sqrt(dx*dx + dy*dy + TINY*TINY);
            ValueType s = x_max - x_min;
            if(s / dist < THETA) {
                ValueType invDist3 = 1.0 / (dist*dist*dist);
                acc_array(idx,0,n,acc) += G*mass*dx*invDist3;
                acc_array(idx,1,n,acc) += G*mass*dy*invDist3;
            } else {
                for(int i=0;i<4;i++) if(children[i]) children[i]->compute_acc(idx,pos,acc,n);
            }
        }
    }
};

// -------------------- Compute accelerations --------------------
void compute_acc_tree(int n, ValueType* pos, ValueType* acc) {
    for(int i=0;i<n*NDIM;i++) acc[i]=0.0;
    ValueType xmin=pos_array(0,0,n,pos), xmax=pos_array(0,0,n,pos);
    ValueType ymin=pos_array(0,1,n,pos), ymax=pos_array(0,1,n,pos);
    for(int i=1;i<n;i++){
        xmin = std::min(xmin,pos_array(i,0,n,pos));
        xmax = std::max(xmax,pos_array(i,0,n,pos));
        ymin = std::min(ymin,pos_array(i,1,n,pos));
        ymax = std::max(ymax,pos_array(i,1,n,pos));
    }
    QuadNode root(xmin,xmax,ymin,ymax);
    for(int i=0;i<n;i++) root.insert(i,pos,n);
    for(int i=0;i<n;i++) root.compute_acc(i,pos,acc,n);
}

// -------------------- Integration --------------------
void integrate(int n, ValueType dt, ValueType* pos, ValueType* vel, ValueType* acc) {
    for(int i=0;i<n;i++) 
        for(int j=0;j<NDIM;j++)
            pos_array(i,j,n,pos) += vel_array(i,j,n,vel)*dt + 0.5*acc_array(i,j,n,acc)*dt*dt;

    compute_acc_tree(n,pos,acc);

    for(int i=0;i<n;i++)
        for(int j=0;j<NDIM;j++)
            vel_array(i,j,n,vel) += 0.5*acc_array(i,j,n,acc)*dt;
}

// -------------------- Velocity statistics --------------------
void velocity_stats(int n, ValueType* vel, ValueType& vmin, ValueType& vmax, ValueType& vave){
    vmin = std::numeric_limits<ValueType>::max();
    vmax = std::numeric_limits<ValueType>::lowest();
    vave = 0.0;
    for(int i=0;i<n;i++){
        ValueType v = std::sqrt(vel_array(i,0,n,vel)*vel_array(i,0,n,vel) +
                                vel_array(i,1,n,vel)*vel_array(i,1,n,vel));
        vmin = std::min(vmin,v);
        vmax = std::max(vmax,v);
        vave += v;
    }
    vave /= n;
}

// -------------------- Main --------------------
int main() {
    int n = 100;
    int nsteps = 100;
    ValueType dt = 0.01;

    std::cout << "Number Objects = " << n << "\n";
    std::cout << "Number Steps   = " << nsteps << "\n";
    std::cout << "Timestep size  = " << dt << "\n";
    std::cout << "Alignment      = " << sizeof(ValueType) << " bytes\n";
    std::cout << "ValueType      = double\n";
    std::cout << "Format         = StructureOfArrays\n";

    std::vector<ValueType> pos(n*NDIM);
    std::vector<ValueType> vel(n*NDIM);
    std::vector<ValueType> acc(n*NDIM);

    // Initialize positions and velocities
    for(int i=0;i<n;i++){
        pos_array(i,0,n,pos) = 2.0*(drand48()-0.5);
        pos_array(i,1,n,pos) = 2.0*(drand48()-0.5);
        vel_array(i,0,n,vel) = drand48()-0.5;
        vel_array(i,1,n,vel) = drand48()-0.5;
    }

    compute_acc_tree(n,pos.data(),acc.data());

    ValueType total_time_ms = 0.0;
    for(int step=0; step<nsteps; step++){
        auto t0 = std::chrono::high_resolution_clock::now();
        integrate(n,dt,pos.data(),vel.data(),acc.data());
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<ValueType,std::milli> step_dur = t1-t0;
        total_time_ms += step_dur.count();

        if(step % 10 == 0) { // Print every 10 steps
            ValueType vmin,vmax,vave;
            velocity_stats(n,vel.data(),vmin,vmax,vave);
            std::cout << "min/max/ave velocity = " << vmin << ", " << vmax << ", " << vave << "\n";
        }
    }

    ValueType avg_time = total_time_ms / nsteps;
    std::cout << "Average time = " << avg_time << " (ms) per step with " << n << " elements\n";

    return 0;
}