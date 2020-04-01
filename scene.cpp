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

struct Material{
	vec color;
	float shine;
	float min_light;
	float mat;
	Material(const vec& c, const float s, const float m, float mm=0.f): color(c), shine(s), min_light(m) {
		mat = mm;
	};
};

struct Sphere{
	//vec color;
	Material material;
	vec center;
	float radius;

	Sphere(const vec& c, float r, const Material& cc): center(c), radius(r), material(cc) {} 
};

struct Light {
	float brightness;
	vec location;
	Light(vec a, float b) : location(a), brightness(b) {}
};

struct Scene {
	std::vector<Sphere> spheres;
	std::vector<Light> lamps;
	vec background;

	Scene(vec c): background(c){}

	bool condition_intersection_sphere(const vec& camera, const vec& dir, const Sphere & s, vec& hit, float& dist) {
		vec vco = s.center - camera;
	    vec point_norm = camera + (dir * vco) * dir;
        if (!((s.center - point_norm).norm() > s.radius)) {
        	if (!(((dir * vco) < 0) && ((vco * vco - s.radius * s.radius) > 0))) {
				if (sqrt(vco * vco - s.radius * s.radius) < dist){
    				dist = sqrt(vco * vco - s.radius * s.radius);
    				float dist_hit = sqrt(s.radius * s.radius - (s.center - point_norm)*(s.center - point_norm));
    				hit = point_norm - dir*(dist_hit);
    				return true;
	        	}
        	} 
        }
        return false;
	}

	bool intersection_ray(const vec& camera, const vec& dir, vec& hit,  vec& normalb, Material& material_pixel) {
		bool flag_sphere = false;
		float min_dist = numeric_limits<float>::max();
		for (int i = 0; i < spheres.size(); i++) {
			if (condition_intersection_sphere(camera, dir, spheres[i], hit, min_dist)) {
				material_pixel = spheres[i].material;
				flag_sphere = true;
				normalb = (hit - spheres[i].center).normalize();
	        }
		}
		return flag_sphere;
	}

	bool skip_light(int ind, vec& hit) {
		float old_dist = (hit - lamps[ind].location).norm();
		vec new_dir = (hit - lamps[ind].location).normalize();
		vec hit_new;
		for (int i = 0; i < spheres.size(); i++) {
			float dist = numeric_limits<float>::max();
			if (condition_intersection_sphere(lamps[ind].location, new_dir, spheres[i], hit_new, dist)) {
				if (old_dist > dist) {
					return true;
				}
			}
		}
		return false;
	}

	vec lighting_scene(const vec& camera, const vec& dir) {
		float sum_brightness = 0.f, sum_shine = 0.f;		
		vec normalb;
		Material material_pixel = Material(background, 0.f, 0.0f);
		vec hit;
		if (intersection_ray(camera, dir, hit, normalb, material_pixel)){
			for (int i = 0; i < lamps.size(); i++){
				if (skip_light(i, hit)) {
					continue;
				}
				vec a =  (lamps[i].location - hit).normalize();
				sum_brightness += lamps[i].brightness * (a * normalb);
				vec med_vec = ((lamps[i].location - hit).normalize() + (camera - hit).normalize()).normalize();
				sum_shine += pow((med_vec * normalb), material_pixel.shine)* lamps[i].brightness;
			}

			vec res = (material_pixel.color * max(material_pixel.min_light, min(1.f , sum_brightness)) + vec(255., 255., 255.)*sum_shine*material_pixel.mat);
			if ( (res[0] >= 256) || (res[1] >= 256) || (res[2] >= 256)){
				res = vec(min(255.f, res[0]), min(255.f, res[1]), min(255.f, res[2]));
			}
			return res; 
		}
		return material_pixel.color;
	}

	void print_scene(vec camera){
		float fov = M_PI / 6;
		std::vector<vec> frame(WIDTH * HEIGHT);
		for (int i = 0; i < WIDTH; i++){
			for (int j = 0; j < HEIGHT; j++){
				float new_height = 2 * tan(fov / 2.);
				float new_width = new_height * (float)WIDTH / (float)HEIGHT;
				float step = new_height / (float)HEIGHT;
				float x = (j + 0.5) * step - new_height / 2;
				float y = -((i + 0.5) * step) + new_width / 2;
            	vec dir = vec(x, y, -1).normalize();
            	//vec N;
            	frame[i+j*WIDTH] = lighting_scene(camera, dir);
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
	vec camera = vec(0., 0., 10.);
	Scene my_scene(vec(0., 0., 0.));
	Material m1 = Material(vec(255, 255, 255), 50.f, 0.05f, 1.f);
	Material m2 = Material(vec(255, 255 ,0.), 100.0f, 0.05f, 1.f);
	Material m3 = Material(vec(0., 0., 255), 60.f, 0.05f, 1.f);
	Material m4 = Material(vec(0., 255., 255), 11.f, 0.05f, 1.f);
	my_scene.spheres.push_back(Sphere(vec(-3,    0,   -16), 2, m1));
	my_scene.spheres.push_back(Sphere(vec(-1.0, -1.5, -12), 2, m2));
	my_scene.spheres.push_back(Sphere(vec(1.5, -0.5, -18), 3, m3));
	my_scene.spheres.push_back(Sphere(vec(7., 5., -18), 3, m4));
	my_scene.lamps.push_back(Light(vec(-30, 20,  20), 1.));//camera, 1.));//
	my_scene.lamps.push_back(Light(vec(10, -20,  40), 0.4));
	my_scene.lamps.push_back(Light(camera, 0.02));
	my_scene.print_scene(camera);
	return 0;
}