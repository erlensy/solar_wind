#include <stdio.h>
#include <math.h>

struct Vec {
    double x, y, z;
};

double length(struct Vec v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double dot(struct Vec u, struct Vec v) {
    return v.x * u.x + v.y * u.y + v.z * u.z;
}

struct Vec cross(struct Vec v, struct Vec u) {
   struct Vec w;
   w.x = v.y * u.z - (v.z * u.y);
   w.y = - (v.x * u.z - (v.z * u.x));
   w.z = v.x * u.y - (v.y * u.x);
   return w;
}

struct Vec mag_field(struct Vec r) {
    double r_length = length(r);
    double r_third = r_length ** 3;
    r.x /= r_length; r.y /= r_length; r.z /= r_length;
    double 3_m_dot_r_hat = 3.0 * dot(m, r);
    
    struct Vec B;
    B.x = (r.x * 3_m_dot_r_hat - m.x) / r_third;
    B.y = (r.y * 3_m_dot_r_hat - m.y) / r_third;
    B.z = (r.z * 3_m_dot_r_hat - m.z) / r_third;
    return B;
}

void RK4(int steps, double h, struct Vec r0, struct Vec v0) {
    // solution array for r and v
    struct Vec r[steps + 1];
    struct Vec v[steps + 1];
    
    // initial conditions
    r[0] = r0;
    v[0] = v0;

    for (unsigned int i = 0; i < steps; i++) {




}


int main() {
    struct Vec v;
}
