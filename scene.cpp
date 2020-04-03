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
	    /*if (!norma) {
	    	cout << "noooooooooooorm" << endl;
	    }*/
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
	float matte;
	float mirror;
	float color_saturation;
	float refraction;
	float refraction_saturation;

	Material(const vec& color1, vec setting, const float shine1, const float min_light1, float refraction1=0.f, float refraction_saturation1=0.f): color(color1), shine(shine1), min_light(min_light1) {
		color_saturation = setting[0];
		matte = setting[1];
		mirror = setting[2];
		
		refraction = refraction1;
		refraction_saturation = refraction_saturation1;
	};
};
struct Plane {
	vec position;
	Material material_table;
	Plane(const vec& pos, const Material& mater) : position(pos), material_table(mater) {}
};

struct Sphere{
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
	int depth_scene;
	float eps;
	float eps2;
	bool WhereAmI; //in air (true) or in object (false)
	Plane plane;
	Scene(vec c, int d, const Plane& p, float e=0.01, float ee=0.001, bool where=true): plane(p), background(c), WhereAmI(where), depth_scene(d) {
		eps = e;
		eps2 = ee;
	}

	bool condition_intersection_sphere(const vec& camera, const vec& dir, const Sphere & s, vec& hit, float& dist) {
		vec vco = s.center - camera;
	    vec point_norm = camera + (dir * vco) * dir;
        if (!((s.center - point_norm).norm() > s.radius)) {
        	if (!(((dir * vco) < 0) && ((vco * vco - s.radius * s.radius) > 0))) {
				if ((vco * vco - s.radius * s.radius) > eps ) {
					if (sqrt(vco * vco - s.radius * s.radius) < dist) {
	    				dist = sqrt(vco * vco - s.radius * s.radius);
	    				float dist_hit = sqrt(s.radius * s.radius - (s.center - point_norm)*(s.center - point_norm) + eps2);
	    				hit = point_norm - dir*(dist_hit);
	    				WhereAmI = true;
	    				return true;
		        	} 
		        } else {
		        	//if (sqrt(vco * vco)< dist) {
        				dist = sqrt(vco * vco);
        				float dist_hit = 2.f * sqrt(s.radius * s.radius - (s.center - point_norm)*(s.center - point_norm));
        				hit = camera + dir*(dist_hit);
        				WhereAmI = false;
        				return true;
	        		//}/
		        }
        	} 
        }
        return false;
	}

	bool intersection_ray(const vec& camera, const vec& dir, vec& hit,  vec& normalb, Material& material_pixel ) {
		bool flag_sphere = false;
		float min_dist = numeric_limits<float>::max();
		for (int i = 0; i < spheres.size(); i++) {
			if (condition_intersection_sphere(camera, dir, spheres[i], hit, min_dist)) {
				material_pixel = spheres[i].material;
				flag_sphere = true;
				normalb = (hit - spheres[i].center).normalize();
	        }
		}
/*
		float plane_dist = numeric_limits<float>::max();
		vec a = vec(camera.x - plane.position.x, 0., 0.);
		float cosFi  = a * dir;
		if (cosFi > 0) {
			return flag_sphere;
		}
		float len_hypotenuse = a.norm() * (1 / cosFi);
		vec table_hit = dir * len_hypotenuse + camera;
		//cout << len_hypotenuse << min_dist << endl;
		if ((len_hypotenuse < min_dist) & (100<len_hypotenuse < 1000)) {
			hit = table_hit;
			normalb = vec(0., -1., 0.);
			material_pixel = plane.material_table;
			flag_sphere = true;
		}
*/

		/*float plane_dist = numeric_limits<float>::max();
		vec a = camera - plane.position;
		float cosFi  = abs(a.normalize() * dir);
		//cout << cosFi << endl;
		float len_hypotenuse = a.norm() * (1 / cosFi);
		vec table_hit = dir * len_hypotenuse + camera;
		//cout << len_hypotenuse << min_dist << endl;
		if (len_hypotenuse < min_dist) {
			hit = table_hit;
			normalb = vec(0., 0., 1.);
			material_pixel = plane.material_table;
			flag_sphere = true;
		}*/
		float dist_plane = std::numeric_limits<float>::max();
  		if (fabs(dir.x)>0.00001)  {
  			float d = (camera.x+plane.position.x)/(dir.x); // the checkerboard plane has equation y = -4
		    if (d > 0) {
			    vec hit_plane = camera + dir*d;
			    //cout << hit_plane.x << endl;
			    if (fabs(hit_plane.y)<30 && hit_plane.z<camera.z && hit_plane.z>-600 && d<min_dist) {
			        dist_plane = d;
			        hit = hit_plane;
			        normalb = vec(-1,0,0);
			        material_pixel = plane.material_table;
			       // cout << int(.5*hit.y+1) + int(.5*hit.z) << endl;
			        material_pixel.color = ((int(0.25 * hit.y - 30) + int(0.25 * hit.z)) & 1) ? vec(255,255,255) : vec(0, 0, 0);
			        //material_pixel.color = material_pixel.color*.3;
			    }
			}
		}
		return std::min(min_dist, dist_plane) < 100000;

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

	vec refract(const Material& material_pixel, const vec& dir, const vec& normalb) {
		float n_1 = WhereAmI > 0.f ? 1.f : material_pixel.refraction_saturation;
		float n_2 = WhereAmI > 0.f ?  material_pixel.refraction_saturation : 1.f;
		vec refraction_dir;
		float cosFi1 = (normalb * dir) < 0 ? -(normalb * dir) : (normalb * dir);
		if (!(1 - (1 - pow(cosFi1, 2)) * pow(n_1 / n_2, 2) < 0)) {
			float cosFi2 = sqrtf(1 - (1 - pow(cosFi1, 2)) * pow(n_1 / n_2, 2));
			vec parallel_transfer = (dir + (cosFi1) * normalb).normalize();
			refraction_dir = (dir *(n_1 / n_2) + normalb * ((n_1 / n_2) * (cosFi1) - cosFi2)).normalize(); 
		} else {
			 refraction_dir = vec(0., 0., 0.);
		}
		return refraction_dir;
	}
	
	vec lighting_scene(const vec& camera, const vec& dir, int depth=0) {
		float sum_brightness = 0.f, sum_shine = 0.f;		
		vec normalb;
		Material material_pixel = Material(background, vec(1., 0., 0.), 0.f, 0.0f);
		vec hit;
		if ((depth < depth_scene) && intersection_ray(camera, dir, hit, normalb, material_pixel)) {
			vec hit_shift = hit + eps2 * normalb;
			vec reflection_color;
			vec refraction_color;
			if (material_pixel.mirror) {
				vec reflection_dir = dir - 2.f * (normalb * dir) * normalb;
				reflection_color =  lighting_scene(hit_shift, reflection_dir, depth + 1);
		    }
			if (material_pixel.refraction) {
				refraction_color =  lighting_scene(hit_shift, refract(material_pixel, dir, normalb), depth + 1);
			}
			for (int i = 0; i < lamps.size(); i++){
				if (skip_light(i, hit)) {
					continue;
				}
				vec a =  (lamps[i].location - hit).normalize();
				sum_brightness += lamps[i].brightness * (a * normalb);
				vec med_vec;
				if ((camera - hit).norm() < 0.1) {
					vec med_vec = (lamps[i].location - hit).normalize();
				} else med_vec = ((lamps[i].location - hit).normalize() + (camera - hit).normalize()).normalize();
				sum_shine += pow((med_vec * normalb), material_pixel.shine)* lamps[i].brightness;
			}
			vec res = material_pixel.color * max(material_pixel.min_light, min(1.f , sum_brightness))*material_pixel.color_saturation;
			res = res + refraction_color*material_pixel.refraction;
			res = res + vec(255., 255., 255.)*sum_shine*material_pixel.matte;
	        res =  res + reflection_color*material_pixel.mirror;
			
			//cout << refraction_color << endl;
			if ( (res[0] >= 256) || (res[1] >= 256) || (res[2] >= 256)) {
				res = vec(min(255.f, max(0.f, res[0])), min(255.f, max(0.f,res[1])), min(255.f, max(0.f,res[2])));
			}
			return res; 
		}
		return material_pixel.color;
	}

	void print_scene(vec camera){
		float fov = M_PI / 6;
		std::vector<vec> frame(WIDTH * HEIGHT);
		for (int i = 0; i < WIDTH; i++) {
			for (int j = 0; j < HEIGHT; j++) {
				WhereAmI = true;
				float new_height = 2 * tan(fov / 2.);
				float new_width = new_height * (float)WIDTH / (float)HEIGHT;
				float step = new_height / (float)HEIGHT;
				float x = (j + 0.5) * step - new_height / 2;
				float y = -((i + 0.5) * step) + new_width / 2;
            	vec dir = vec(x, y, -1).normalize();
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
	//vec color; const float shine; const float min_light; float matte1=0.f; float mirror1=0.f;
	//float color_saturation1=1.f, float refraction ; float refraction_saturation)
	vec setting1 =  vec(1.,1., 0.3); //c_s, no_matte (0), mirror
	//fdsfsdfdsfsdfdkffdkjngkjfbfjhvbhfbvhkdbfhvbfdkjbvjds

	vec setting2 =  vec(0. ,0.5, 0.);
	//setting2 =  vec(0.,0., 0.);
	vec setting3 =  vec(1.,0.5, 1.0);
	vec setting4 =  vec(1.,0., 0.);
	vec camera = vec(3., 0., 30.);
	Material m1 = Material(vec(255, 217, 25),  setting1, 50.f,  0.05f);
	Material m2 = Material(vec(255, 255 ,0.),  setting2, 5000.0f, 0.f,   0.99f,  1.55f); //0.f, 1.55f);
	//m2 = Material(vec(255, 255 ,0.),  setting2, 60.0f, 0.0f,   1.f,  1.33f); 
	Material m3 = Material(vec(0., 0., 255),   setting3, 50.f,  0.1f, 0.f);
	Material m4 = Material(vec(0., 255., 255), setting4, 100.f, 0.05f);

	vec setting_for_plane1 =  vec(0.8,0.8, 0.5);
	Material m_for_plane1 = Material(vec(255., 0., 255),   setting_for_plane1, 50.f,  0.05f);
	Plane plane_chess =Plane(vec(4., 0., -1.), m_for_plane1);
	Scene my_scene(vec(255., 182., 193.), 4, plane_chess);//shine min_light refrac refr_sat


	my_scene.spheres.push_back(Sphere(vec(-3,    0,   -16), 2, m1));
	my_scene.spheres.push_back(Sphere(vec(0., -3, 0), 2, m2));
	my_scene.spheres.push_back(Sphere(vec(1.5, -0.5, -18), 3, m3));
	my_scene.spheres.push_back(Sphere(vec(7., 5., -18), 3, m4));
	my_scene.lamps.push_back(Light(vec(-30, 20,  20), 1.));//camera, 1.));//
	my_scene.lamps.push_back(Light(vec(10, -20,  40), 0.4));
	my_scene.lamps.push_back(Light(camera, 0.02));
	my_scene.print_scene(camera);
	return 0;
}