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
	Material(const vec& c, const float s): color(c), shine(s){};
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

	bool intersection_ray(const vec& camera, const vec& dir, vec& hit,  vec& normalb, Material& material_pixel){
		float min_dist = numeric_limits<float>::max();
		bool flag_sphere = false;
		float dist_hit;
		for (int i = 0; i < spheres.size(); i++) {
			vec vco = spheres[i].center - camera;
	        vec point_norm = camera + (dir * vco) * dir;

	        if (!((spheres[i].center - point_norm).norm() > spheres[i].radius)) {
	        	if (!(((dir * vco) < 0) && ((vco * vco - spheres[i].radius * spheres[i].radius) > 0))) {
	        		if ((vco * vco - spheres[i].radius * spheres[i].radius) < 0) {
	        			if (sqrt(vco * vco)< min_dist){
	        				material_pixel = spheres[i].material;
	        				flag_sphere = true;
	        				min_dist = sqrt(vco * vco);
	        				dist_hit = sqrt(spheres[i].radius * spheres[i].radius - (spheres[i].center - point_norm)*(spheres[i].center - point_norm));
	        				hit = camera + dir*(dist_hit);
            				normalb = (hit - spheres[i].center).normalize();
	        			}
	        		} else {
	        			if (sqrt(vco * vco - spheres[i].radius * spheres[i].radius) < min_dist){
	        				material_pixel = spheres[i].material;
	        				flag_sphere = true;
	        				min_dist = sqrt(vco * vco - spheres[i].radius * spheres[i].radius);
	        				dist_hit = sqrt(spheres[i].radius * spheres[i].radius - (spheres[i].center - point_norm)*(spheres[i].center - point_norm));
	        				hit = point_norm - dir*(dist_hit);
            				normalb = (hit - spheres[i].center).normalize();
	        			}
	        		}
	        	} 
	        }
		}
		return flag_sphere;
	}

	vec lighting_scene(const vec& camera, const vec& dir){
		float sum_brightness = 0., sum_shine = 0;		
		vec normalb;
		Material material_pixel = Material(background, 0.f);
		vec hit;
		if (intersection_ray(camera, dir, hit, normalb, material_pixel)){
			for (int i = 0; i < lamps.size(); i++){
				vec a =  (lamps[i].location - hit).normalize();
				sum_brightness += lamps[i].brightness * (a * normalb);
				vec med_vec = ((lamps[i].location - hit).normalize() + (camera - hit).normalize()).normalize();
				sum_shine += pow((med_vec * normalb), material_pixel.shine)* lamps[i].brightness;
			}
			vec res = (material_pixel.color * max(0.05f, min(1.f , sum_brightness)) + vec(255., 255., 255.)*sum_shine);
			if ( (res[0] >= 256) || (res[1] >= 256) || (res[2] >= 256)){
				//cout<<res<<endl;
				res = vec(min(255.f, res[0]), min(255.f, res[1]), min(255.f, res[2]));
				
				//res = vec(255, 255, 255);
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
				float x = j * step - new_height / 2;
				float y = -(i * step) + new_width / 2;
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
	vec camera = vec(0.,0.,10.);
	Scene my_scene(vec(0.,0.,0.));
	Material m1 = Material(vec(255,255, 255), 50.f);
	Material m2 = Material(vec(255,255,0.), 100.0f);
	Material m3 = Material(vec(0.,0.,255), 60.f);
	Material m4 = Material(vec(0.,255.,255), 11.f);
	my_scene.spheres.push_back(Sphere(vec(-3,    0,   -16), 2, m1));
	my_scene.spheres.push_back(Sphere(vec(-1.0, -1.5, -12), 2, m2));
	my_scene.spheres.push_back(Sphere(vec(1.5, -0.5, -18), 3, m3));
	my_scene.spheres.push_back(Sphere(vec(7., 5., -18), 3, m4));
	my_scene.lamps.push_back(Light(vec(-20, 20,  20), 1.));//camera, 1.));//
	my_scene.lamps.push_back(Light(vec(10, -20,  40), 0.3));
	my_scene.lamps.push_back(Light(camera, 0.02));
	my_scene.print_scene(camera);
	return 0;
}