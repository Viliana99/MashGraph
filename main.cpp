#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <vector>

using namespace std;

class vec {

  private:
  public:
    float x;
    float y;
    float z;
 

    vec(float a=0,float b=0, float c=0): x(a),y(b),z(c){};
    
    float norm() const {
        return sqrt(x*x + y*y + z*z); 
    }
    
    vec& normalize() {
        float norma = (*this).norm();
        this -> x = x / norma;
        this -> y = y / norma; 
        this -> z = z / norma;
        return *this;
    }
    
    float operator[](size_t i) {
        assert(i < 3);
        return i == 0 ? x : (i == 1 ? y : z);
    }
    
    friend vec operator+(const vec& a, const vec& b);

    friend vec operator-(const vec& a, const vec& b);
    
    friend float operator*(const vec& a, const vec& b);
    
    friend vec operator*(float a, const vec& b);
    
    friend vec operator*(const vec& a, float b);
    
    friend ostream& operator<<(ostream& out, const vec& v);
};





vec operator+(const vec& a, const vec& b) {
        return vec(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec operator-(const vec& a, const vec& b) {
        return vec(a.x - b.x, a.y - b.y, a.z - b.z);
}

float operator*(const vec& a, const vec& b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

vec operator*(float a, const vec& b) {
    return vec(a * b.x, a * b.y, a * b.z);
}

vec operator*(const vec& a, float b) {
    return vec(a.x * b, a.y * b, a.z * b);
}

ostream& operator<<(ostream& out, const vec& v){
    out << '(' << v.x << ", " << v.y  << ", " << v.z << ')';
    return out;
}

struct Material {
    vec color;
    Material(const vec &color1) : color(color1) {}
    Material() : color() {}
};


struct Sphere{
    vec center;
    float radius;
    Material material;

    Sphere(const vec &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}
    /*
    bool ray_intersection(const vec &orig, const vec &dir, float &t0) const {
        vec L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) zreturn false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    } */
    bool ray_intersection(const vec &orig, const vec &dir) const {
        vec vco = center - orig;
        float norm_dir = dir.norm();
        vec ndir = vec(dir.x / norm_dir, dir.y / norm_dir, dir.z / norm_dir);
        vec point_norm = orig + (ndir * vco) * ndir;
        if ((center - point_norm).norm() > radius) {
            return false;
        }
        return true;
    }
    friend vec cast_ray(const vec &orig, const vec &dir, const Sphere &sphere);
};


bool scene_intersection(const vec &orig, const vec &dir, const vector<Sphere> &spheres, vec &hit, vec &N, Material &material) {
    float spheres_dist = numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersection(orig, dir) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    return spheres_dist<1000;
}
    vec cast_ray(const vec &orig, const vec &dir, const vector<Sphere> &spheres) {
        vec point, N;
        Material material;

        if (scene_intersection(orig, dir, spheres, point, N, material)) {
            return material.color;
        }
        return vec(0.9, 0.9, 0.6); // background color
    }


void render(const vector<Sphere> &spheres) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    vector<vec> framebuffer(width*height);

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width
              - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            vec dir = vec(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(vec(0,0,0), dir, spheres);
        }
    }

    ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}


int main() {
    Material      ivory(vec(0.4, 0.4, 0.3));
    Material red_rubber(vec(0.3, 0.1, 0.1));

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(vec(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(vec( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(vec( 7, 5, -18), 4,      ivory));

    render(spheres);

    return 0;
}