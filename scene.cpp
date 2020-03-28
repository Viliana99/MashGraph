#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <vector>




#define WIDTH 1024
#define HEIGHT 768

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
	
	friend float operator*(const vec& a, const vec& b);
	
	friend vec operator*(float a, const vec& b);
	
	friend vec operator+(const vec& a, const vec& b);

	friend vec operator-(const vec& a, const vec& b);

	friend vec operator*(const vec& a, float b);
	

    vec& normalize() {
	    float norma = (*this).norm();
	    this -> x = x / norma;
	    this -> y = y / norma; 
	    this -> z = z / norma;
	    return *this;
    }

    friend ostream& operator<<(ostream& out, const vec& v);
	
	float operator[](size_t i) {
	    assert(i < 3);
	    return i == 0 ? x : (i == 1 ? y : z);
    }
	
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
	out << v.x << ' ' << v.y  << " " << v.z << endl;
	return out;
}
/*
struct Object{
	vec color;
};*/
struct Sphere{
	vec color;
	vec center;
	float radius;
	Sphere(vec c, float r, vec cc): center(c), radius(r), color(cc) {} 
};

struct Scene {
	std::vector<Sphere> spheres;
	vec background;

	Scene(vec c): background(c){}

	vec intersection_ray(const vec& camera, const vec& dir,  vec& normalb){
		vec color_pixel = background;
		float min_dist = numeric_limits<float>::max();
		for (int i = 0; i < spheres.size(); i++) {
			bool intersection;
			vec vco = spheres[i].center - camera;
	        vec point_norm = camera + (dir * vco) * dir;
	        if (!((spheres[i].center - point_norm).norm() > spheres[i].radius)) {
	        	if (!(((dir * vco) < 0) && ((vco * vco - spheres[i].radius * spheres[i].radius) > 0))) {
	        		if ((vco * vco - spheres[i].radius * spheres[i].radius) < 0) {
	        			if (sqrt(vco * vco)< min_dist){
	        				color_pixel = spheres[i].color;
	        				min_dist = sqrt(vco * vco);
	        			}
	        		} else {
	        			if (sqrt(vco * vco - spheres[i].radius * spheres[i].radius) < min_dist){
	        				color_pixel = spheres[i].color;
	        				min_dist = sqrt(vco * vco - spheres[i].radius * spheres[i].radius);
	        			}
	        		}
	        	} 
	        }
		}
		return color_pixel;
	}

	void print_scene(vec camera){
		float fov = M_PI / 6;
		std::vector<vec> frame(WIDTH * HEIGHT);
		for (int i = 0; i < WIDTH; i++){
			for (int j = 0; j < HEIGHT; j++){
				float new_height = 2 * tan(fov / 2.);
				float new_width = new_height * (float)WIDTH / (float)HEIGHT;
				float step = new_height / (float)HEIGHT;
				float x = j * step - new_height / 2;
				float y = -(i * step) + new_width / 2;
            	vec dir = vec(x, y, -1).normalize();
            	vec N;
            	frame[i+j*WIDTH] = intersection_ray(camera, dir, N);
			}
		}
		std::ofstream ofs; 
	    ofs.open("./out.ppm");
	    ofs << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";
	    for (int i = 0; i < HEIGHT*WIDTH; i++) {
	            ofs << (char)(frame[i][0]);
	            ofs << (char)(frame[i][1]);
	            ofs << (char)(frame[i][2]);
	    }
	    ofs.close();
	}

};


int main(){
	Scene my_scene(vec(0.,0.,0.));
	my_scene.spheres.push_back(Sphere(vec(0, 0, 25), 5, vec(255,255, 255)));
	my_scene.spheres.push_back(Sphere(vec(10, 10, 50), 10, vec(255,255,0.)));
	my_scene.spheres.push_back(Sphere(vec(-10, -5, 5), 2, vec(0.,0.,255)));
	//my_scene.spheres.push_light(Light(Vect3D(0, 0, 0), 0.5f));
	//my_scene.spheres.push_light(Light(Vect3D(5, 30, 30), 1.f));
	my_scene.print_scene(vec(0.,0.,80.));
	return 0;
}