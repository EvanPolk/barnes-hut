// barnes_hut_nbody.cpp
// Compile: g++ -std=c++17 -O2 barnes_hut_nbody.cpp -o barnes_hut_nbody
// Run: ./barnes_hut_nbody

#include <bits/stdc++.h>
using namespace std;

struct Simulation {
    // N bodies, 2D
    size_t N;
    vector<double> x;   // x positions (size N)
    vector<double> y;   // y positions (size N)
    vector<double> vx;  // x velocities (size N)
    vector<double> vy;  // y velocities (size N)
    vector<double> mass;// masses (size N)

    // simulation parameters
    double G = 1.0;     // gravitational constant
    double eps = 1e-3;  // softening
    double theta = 0.5; // Barnes-Hut opening threshold

    Simulation(size_t n = 0) { resize(n); }

    void resize(size_t n) {
        N = n;
        x.assign(N, 0.0);
        y.assign(N, 0.0);
        vx.assign(N, 0.0);
        vy.assign(N, 0.0);
        mass.assign(N, 1.0);
    }
};

// Quadtree stored in arrays using index-based children: 0..3
struct QuadTree {
    struct Node {
        // bounding box:
        double cx, cy; // center of square
        double half;   // half-size (square)
        // mass properties:
        double mass;
        double com_x, com_y; // center of mass
        // a single body index if this node is a leaf with one body; -1 otherwise
        int body_index;
        // children indices in node array, -1 if empty
        int child[4];
        bool is_leaf;
        Node() {
            cx = cy = 0.0; half = 0.0;
            mass = 0.0; com_x = com_y = 0.0;
            body_index = -1;
            for (int i=0;i<4;i++) child[i] = -1;
            is_leaf = true;
        }
    };

    vector<Node> nodes;
    int root;

    QuadTree() { clear(); }

    void clear() {
        nodes.clear();
        nodes.emplace_back();
        root = 0;
    }

    // Create root bounding box that contains all bodies
    void init_root(double cx, double cy, double half) {
        clear();
        nodes[root].cx = cx;
        nodes[root].cy = cy;
        nodes[root].half = half;
        nodes[root].mass = 0.0;
        nodes[root].com_x = 0.0;
        nodes[root].com_y = 0.0;
        nodes[root].body_index = -1;
        nodes[root].is_leaf = true;
        for (int i=0;i<4;i++) nodes[root].child[i] = -1;
    }

    // Which quadrant (0..3) a point belongs to relative to node center
    inline int quadrant(const Node &n, double px, double py) const {
        int q = 0;
        if (px >= n.cx) q |= 1; // East
        if (py >= n.cy) q |= 2; // North
        return q;
    }

    // Create a new node and return its index
    int new_node(double cx, double cy, double half) {
        Node nd;
        nd.cx = cx; nd.cy = cy; nd.half = half;
        nd.mass = 0.0; nd.com_x = 0.0; nd.com_y = 0.0;
        nd.body_index = -1; nd.is_leaf = true;
        for (int i=0;i<4;i++) nd.child[i] = -1;
        nodes.push_back(nd);
        return (int)nodes.size() - 1;
    }

    // insert body index i at position (px,py) using simulation data
    void insert_body(int node_idx, int body_idx, double px, double py, const Simulation &sim) {
        Node &node = nodes[node_idx];

        // If node is empty leaf (no body), place it here
        if (node.body_index == -1 && node.is_leaf && node.mass == 0.0) {
            node.body_index = body_idx;
            node.mass = sim.mass[body_idx];
            node.com_x = px;
            node.com_y = py;
            return;
        }

        // If node is leaf and contains one body, we must subdivide and reinsert existing body & new body
        if (node.is_leaf) {
            // existing body
            int existing = node.body_index;
            double ex = sim.x[existing];
            double ey = sim.y[existing];

            // subdivide node
            node.is_leaf = false;
            node.body_index = -1;
            double h = node.half / 2.0;
            // create 4 children: q=0..3 mapping:
            for (int q=0;q<4;q++) {
                double ncx = node.cx + ((q & 1) ? h : -h);
                double ncy = node.cy + ((q & 2) ? h : -h);
                node.child[q] = new_node(ncx, ncy, h);
            }
            // insert existing into appropriate child
            int qex = quadrant(node, ex, ey);
            insert_body(node.child[qex], existing, ex, ey, sim);
            // then fall through to insert new body
        }

        // if internal node, insert into the appropriate child
        int q = quadrant(node, px, py);
        if (node.child[q] == -1) {
            double h = node.half / 2.0;
            double ncx = node.cx + ((q & 1) ? h : -h);
            double ncy = node.cy + ((q & 2) ? h : -h);
            node.child[q] = new_node(ncx, ncy, h);
        }
        insert_body(node.child[q], body_idx, px, py, sim);
    }

    // After inserting all bodies, compute mass and center of mass for nodes (bottom-up)
    void compute_mass_COM(int idx, const Simulation &sim) {
        Node &n = nodes[idx];
        if (n.is_leaf) {
            if (n.body_index != -1) {
                int b = n.body_index;
                n.mass = sim.mass[b];
                n.com_x = sim.x[b];
                n.com_y = sim.y[b];
            } else {
                n.mass = 0.0;
                n.com_x = n.com_y = 0.0;
            }
            return;
        }
        double mx = 0.0, my = 0.0, msum = 0.0;
        for (int q=0;q<4;q++) {
            int c = n.child[q];
            if (c != -1) {
                compute_mass_COM(c, sim);
                msum += nodes[c].mass;
                mx += nodes[c].mass * nodes[c].com_x;
                my += nodes[c].mass * nodes[c].com_y;
            }
        }
        n.mass = msum;
        if (msum > 0.0) {
            n.com_x = mx / msum;
            n.com_y = my / msum;
        } else {
            n.com_x = n.com_y = 0.0;
        }
    }

    // Build tree for given simulation (assumes bounding box already set)
    void build(const Simulation &sim) {
        // nodes[ root ] bounding box must be set externally before build()
        // clear children/content (but keep bounding box)
        int r = root;
        nodes[r].mass = 0.0;
        nodes[r].com_x = nodes[r].com_y = 0.0;
        nodes[r].body_index = -1;
        nodes[r].is_leaf = true;
        for (int i=0;i<4;i++) nodes[r].child[i] = -1;
        // insert all bodies
        for (size_t i = 0; i < sim.N; ++i) {
            insert_body(r, (int)i, sim.x[i], sim.y[i], sim);
        }
        // compute masses and center of mass
        compute_mass_COM(r, sim);
    }

    // compute force on body b index; returns pair (fx, fy)
    pair<double,double> compute_force_on_body(int idx, int b, const Simulation &sim) const {
        double fx = 0.0, fy = 0.0;
        compute_force_recursive(idx, b, sim, fx, fy);
        return {fx, fy};
    }

    void compute_force_recursive(int node_idx, int b, const Simulation &sim, double &fx, double &fy) const {
        const Node &n = nodes[node_idx];
        if (n.mass == 0.0) return;
        double px = sim.x[b], py = sim.y[b];
        // If this node is a leaf and contains the same body, skip
        if (n.is_leaf && n.body_index == b) return;

        double dx = n.com_x - px;
        double dy = n.com_y - py;
        double dist = sqrt(dx*dx + dy*dy + sim.eps*sim.eps);

        // opening criterion: s / d < theta, where s = node size (2*half)
        double s = n.half * 2.0;
        if (n.is_leaf || (s / dist) < sim.theta) {
            // treat as single body
            double invr3 = 1.0 / (dist*dist*dist);
            double f = sim.G * n.mass * invr3;
            fx += f * dx;
            fy += f * dy;
        } else {
            // recurse to children
            for (int q=0;q<4;q++) {
                int c = n.child[q];
                if (c != -1) compute_force_recursive(c, b, sim, fx, fy);
            }
        }
    }
};

// Utility: find bounding square (center, half) that contains all points (with margin)
void find_bounding_square(const Simulation &sim, double &cx, double &cy, double &half) {
    if (sim.N == 0) { cx = cy = 0; half = 1.0; return; }
    double minx = sim.x[0], maxx = sim.x[0], miny = sim.y[0], maxy = sim.y[0];
    for (size_t i=1;i<sim.N;i++) {
        minx = min(minx, sim.x[i]);
        maxx = max(maxx, sim.x[i]);
        miny = min(miny, sim.y[i]);
        maxy = max(maxy, sim.y[i]);
    }
    double dx = maxx - minx;
    double dy = maxy - miny;
    half = 0.5 * max(dx, dy);
    // if degenerate, ensure non-zero half
    if (half < 1e-6) half = 1.0;
    cx = 0.5 * (maxx + minx);
    cy = 0.5 * (maxy + miny);
    // add small margin
    half *= 1.1;
}

// advance simulation by dt (single leapfrog-like step: compute forces -> update velocities -> update positions)
void step_simulation(Simulation &sim, QuadTree &tree, double dt) {
    // build tree
    double cx, cy, half;
    find_bounding_square(sim, cx, cy, half);
    tree.init_root(cx, cy, half);
    tree.build(sim);

    // compute all forces
    vector<double> fx(sim.N, 0.0), fy(sim.N, 0.0);
    for (size_t i=0;i<sim.N;i++) {
        auto p = tree.compute_force_on_body(tree.root, (int)i, sim);
        fx[i] = p.first;
        fy[i] = p.second;
    }

    // update velocities and positions (simple explicit Euler or semi-implicit Euler)
    for (size_t i=0;i<sim.N;i++) {
        // acceleration
        double ax = fx[i] / sim.mass[i];
        double ay = fy[i] / sim.mass[i];
        // velocity update (semi-implicit)
        sim.vx[i] += ax * dt;
        sim.vy[i] += ay * dt;
        // position update
        sim.x[i] += sim.vx[i] * dt;
        sim.y[i] += sim.vy[i] * dt;
    }
}

// compute total kinetic + potential energy (approx potential via pairwise for small N or tree for large N)
// Simple O(N^2) potential for demonstration (works for moderate N)
double total_energy(const Simulation &sim) {
    double K = 0.0;
    for (size_t i=0;i<sim.N;i++) {
        K += 0.5 * sim.mass[i] * (sim.vx[i]*sim.vx[i] + sim.vy[i]*sim.vy[i]);
    }
    double U = 0.0;
    for (size_t i=0;i<sim.N;i++) {
        for (size_t j=i+1;j<sim.N;j++) {
            double dx = sim.x[i] - sim.x[j];
            double dy = sim.y[i] - sim.y[j];
            double r = sqrt(dx*dx + dy*dy + sim.eps*sim.eps);
            U -= sim.G * sim.mass[i] * sim.mass[j] / r;
        }
    }
    return K + U;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // simple example: random Plummer-like distribution or just uniform random
    size_t N = 1000;
    double dt = 0.01;
    int steps = 500;

    // create sim
    Simulation sim(N);
    sim.G = 1.0;
    sim.eps = 0.05;
    sim.theta = 0.6;

    // random initialization
    std::mt19937_64 rng(1234567);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);
    std::uniform_real_distribution<double> massd(0.5, 1.5);

    for (size_t i=0;i<N;i++) {
        // place in a disk
        double r = sqrt(ud(rng)*ud(rng) + ud(rng)*ud(rng));
        double angle = ud(rng) * 3.141592653589793 * 2.0;
        double radius = 0.2 + 0.8 * fabs(ud(rng));
        sim.x[i] = radius * cos(angle);
        sim.y[i] = radius * sin(angle);
        // small random velocities
        sim.vx[i] = 0.1 * ud(rng);
        sim.vy[i] = 0.1 * ud(rng);
        sim.mass[i] = massd(rng);
    }

    QuadTree tree;
    // initial energy
    cout << "# Step, Energy\n";
    cout << 0 << ", " << total_energy(sim) << "\n";

    for (int t=1; t<=steps; ++t) {
        step_simulation(sim, tree, dt);
        if (t % 10 == 0) {
            cout << t << ", " << total_energy(sim) << "\n";
        }
    }

    // print final positions of first 10 bodies
    cout << "# Final positions (first 10 bodies):\n";
    for (size_t i=0; i<min<size_t>(10, sim.N); ++i) {
        cout << i << ": " << sim.x[i] << " " << sim.y[i] << "  vx=" << sim.vx[i] << " vy=" << sim.vy[i] << " m=" << sim.mass[i] << "\n";
    }

    return 0;
}