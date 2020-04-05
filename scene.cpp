#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <vector>
#include "omp.h"
#include "head.h"
#include <thread>
#include <string.h> 

#define PRETTY_YELLOW vec(255, 217, 25)
#define BLUE vec(0., 0., 255)
#define SKYBLUE vec(0., 255., 255)
#define LIGHT_GREEN vec(73,255,137)
#define LIGHT_PINK  vec(255., 182., 193.)
#define BLACK vec(0.,0.,0.)
#define WHITE vec(255., 255., 255.)


#define MIRROR Material(BLUE,  vec(1., 0.5, 1.), 50.f,  0.1f, 0.f)
#define REFRACT_GLASS  Material(BLACK,  vec(0., 1., 0.), 5000.0f, 0.f,   0.99f,  1.55f)
#define VELVET Material(SKYBLUE, vec(1., 0., 0.), 100.f, 0.05f)
#define GOLD  Material(PRETTY_YELLOW, vec(1., 1., 0.3), 50.f,  0.05f)
#define IRON Material(BLACK, vec(1., 1., 0.3), 50.f,  0.05f)

#define WIDTH 1280
#define HEIGHT 720

using namespace std;

void save_BMP(const char* fileName, unsigned char *buf, int width=WIDTH, int height=HEIGHT) {
    const int fileHeaderSize = 14;
    const int infoHeaderSize = 40;
    
    unsigned char bmp_file_header[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0, };
	unsigned char bmp_info_header[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, };
	unsigned char bmp_pad[3] = { 0, 0, 0, };
    
    int paddingSize = (4 - (width * 3) % 4) % 4;
    int size = fileHeaderSize + infoHeaderSize + (3 * width + paddingSize) * height;
    
	bmp_file_header[2]  = static_cast<unsigned char>(size      );
	bmp_file_header[3]  = static_cast<unsigned char>(size >>  8);
	bmp_file_header[4]  = static_cast<unsigned char>(size >> 16);
	bmp_file_header[5]  = static_cast<unsigned char>(size >> 24);

	bmp_info_header[4]  = static_cast<unsigned char>(width      );
	bmp_info_header[5]  = static_cast<unsigned char>(width >>  8);
	bmp_info_header[6]  = static_cast<unsigned char>(width >> 16);
	bmp_info_header[7]  = static_cast<unsigned char>(width >> 24);

	bmp_info_header[8]  = static_cast<unsigned char>(height      );
	bmp_info_header[9]  = static_cast<unsigned char>(height >>  8);
	bmp_info_header[10] = static_cast<unsigned char>(height >> 16);
	bmp_info_header[11] = static_cast<unsigned char>(height >> 24);
    
    FILE *file = fopen(fileName, "wb");
    
    if (file)
	{
		fwrite(bmp_file_header, 1, 14, file);
		fwrite(bmp_info_header, 1, 40, file);

		for (int i = 0; i < height; i++)
		{
			fwrite(buf + (width * i * 3), 3, width, file);
			fwrite(bmp_pad, 1, ((4 - (width * 3) % 4) % 4), file);
		}
		fclose(file);
	}
}

void save(const char* fileName, std::vector<vec>& frame) {
	std::ofstream ofs; 
    ofs.open(fileName);
    ofs << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";
    for (int i = 0; i < HEIGHT*WIDTH; i++) {
            ofs << (char)(frame[i][0]);
            ofs << (char)(frame[i][1]);
            ofs << (char)(frame[i][2]);
    }
    ofs.close();
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

struct Cube{
	Material material;
	vec vmin;
	vec vmax;
	Cube(const vec& c, const vec& c1, const Material& cc): vmin(c), vmax(c1), material(cc) {} 
};

struct Light {
	float brightness;
	vec location;
	Light(vec a, float b) : location(a), brightness(b) {}
};

struct Scene {
	std::vector<Sphere> spheres;
	std::vector<Light> lamps;
	std::vector<Cube> cubes;
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
	        
    				dist = sqrt(vco * vco);
    				float dist_hit = 2.f * sqrt(s.radius * s.radius - (s.center - point_norm)*(s.center - point_norm));
    				hit = camera + dir*(dist_hit);
    				WhereAmI = false;
    				return true;
   		        }
        	} 
        }
        return false;
	}

	bool intersection_plane(const vec& camera,const vec& dir, vec& hit, float& min_dist, vec& normalb, Material& material_pixel) {
		float dist_plane = std::numeric_limits<float>::max();
  		if (fabs(dir.x)>0.00001)  {
  			float d = (camera.x+plane.position.x)/(dir.x); 
		    if (d > 0) {
			    vec hit_plane = camera + dir*d;
			   if (d<min_dist) {
			        dist_plane = d;
			        hit = hit_plane;
			        normalb = vec(-1,0,0);
			        material_pixel = plane.material_table;
			        material_pixel.color = ((int(0.25 * hit.y - 30) + int(0.25 * hit.z)) & 1) ? WHITE : BLACK;
			        return true;
			    }
			}
		}
		return false;
	}
	bool intersection_side(vec camera, vec dir, Cube& cube, vec& hit, float& min_dist, vec& normalb) {

		double t_near = numeric_limits<float>::min(),
		t_far = numeric_limits<float>::max();
		float t1, t2;

		for (int i = 0; i < 3; i++) {
			if (abs(dir[i]) >= 0.0001) {
				t1 = (cube.vmin[i] - camera[i]) / dir[i];
				t2 = (cube.vmax[i] - camera[i]) / dir[i];

				if (t1 > t2)
					std::swap(t1, t2);
				if (t1 > t_near)
					t_near = t1;
				if (t2 < t_far)
					t_far = t2;
				if (t_near > t_far)
					return false;
				if (t_far < 0.0)
					return false;
			} 
			else {
				if (camera[i] < cube.vmin[i] || camera[i] > cube.vmax[i])
					return false;
			}
		} 
		if (t_near <= t_far && t_far >=0 && t_near < min_dist) {
			min_dist = t_near;
			hit = camera + min_dist * dir;
			vec tmp = vec(0., 0., 0.);
			for (int i = 0; i < 3; i++) {
				if (abs(hit[i] - cube.vmin[i]) < eps2) 
					tmp[i] = -1.;
				if (abs(hit[i] - cube.vmax[i]) < eps2)
					tmp[i] = 1.;
			}
			normalb = tmp;
			return true;
		}
		else return false;
	}

	bool intersection_cube(const vec& camera,const vec& dir, vec& hit, float& min_dist, vec& normalb, Material& material_pixel) {
		bool flag_cube = false;
		for(int i = 0; i < cubes.size(); i++) {
			if (intersection_side(camera, dir, cubes[i], hit, min_dist, normalb)) {
				flag_cube = true;
				material_pixel = cubes[i].material;
			}
		}
		return flag_cube;
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
		flag_sphere = flag_sphere || intersection_plane(camera, dir, hit, min_dist, normalb, material_pixel);
		flag_sphere = flag_sphere || intersection_cube(camera, dir, hit, min_dist, normalb, material_pixel);
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
		refraction_dir = WhereAmI ? refraction_dir - eps * normalb : refraction_dir + eps * normalb;

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
				vec refraction_dir = refract(material_pixel, dir, normalb);
				hit_shift = hit*normalb < 0 ? hit +  eps * normalb : hit -  eps * normalb;
				//refraction_dir = hit_shift * normalb ? refraction_dir - eps * normalb : refraction_dir + eps * normalb;
				refraction_color =  lighting_scene(hit_shift, refraction_dir, depth + 1);
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
			if ( (res[0] >= 256) || (res[1] >= 256) || (res[2] >= 256)) {
				res = vec(min(255.f, max(0.f, res[0])), min(255.f, max(0.f,res[1])), min(255.f, max(0.f,res[2])));
			}
			return res; 
		}
		return material_pixel.color;
	}

	void print_scene(vec camera, const char* fileName="out.bmp", bool BMP=true, int num_threads=8){
		float fov = M_PI / 6;
		std::vector<vec> frame(WIDTH * HEIGHT);
		unsigned char bmp_frame[HEIGHT][WIDTH][3];
		int i, j;
		omp_set_num_threads(10);
        #pragma omp parallel private(i, j)
        {


        #pragma omp for collapse(2) schedule(dynamic)

		for (i = 0; i < WIDTH; i++) {
			for (j = 0; j < HEIGHT; j++) {
				WhereAmI = true;
				float new_height = 2 * tan(fov / 2.);
				float new_width = new_height * (float)WIDTH / (float)HEIGHT;
				float step = new_height / (float)HEIGHT;
				vec res;



				for (float i_f = 0.33; i_f < 1.f; i_f = i_f + 0.34) {
					for (float j_f = 0.33; j_f < 1.f; j_f = j_f + 0.34) {
						float x = (j + j_f) * step - new_height / 2;
						float y = -((i + i_f) * step) + new_width / 2;
		            	x = -x;
		            	vec dir = vec(x, y, -1).normalize();
		            	res = res +  lighting_scene(camera, dir);
					}
				}
				res = (0.25) * res;
				if (!BMP) {
					frame[i+j*WIDTH] = res; 
				} else {
					bmp_frame[j][i][2] = (unsigned char) res.x; ///red
	                bmp_frame[j][i][1] = (unsigned char) res.y; ///green
	                bmp_frame[j][i][0] = (unsigned char) res.z; ///blue
				} 
			}
		}

    	}
		if (BMP) {
			save_BMP(fileName, (unsigned char *)bmp_frame);
		} else {
			save(fileName, frame);
		}
	}
};



int main(int argc, char* argv []){
	vec camera = vec(3., 0., 30.);

	vec setting_for_plane1 =  vec(0.8,0.8, 0.5);
	Material m_for_plane1 = Material(BLACK,   setting_for_plane1, 50.f,  0.05f);
	Plane plane_chess =Plane(vec(4., 0., -1.), m_for_plane1);
	Scene my_scene(LIGHT_PINK, 4, plane_chess);//shine min_light refrac refr_sat


	my_scene.spheres.push_back(Sphere(vec(-3,    0,   -16), 2, GOLD));
	my_scene.spheres.push_back(Sphere(vec(0., -2, 0), 2, REFRACT_GLASS));
	my_scene.spheres.push_back(Sphere(vec(1.5, -0.5, -18), 3, MIRROR));
	my_scene.spheres.push_back(Sphere(vec(7., 5., -18), 3, VELVET));
	my_scene.lamps.push_back(Light(vec(-30, 20,  20), 1.));//camera, 1.));//
	my_scene.lamps.push_back(Light(vec(10, -20,  40), 0.4));
	my_scene.lamps.push_back(Light(camera, 0.02));
	my_scene.cubes.push_back(Cube(vec(-4,    -4,   -4), vec(-2,    -2,   -2), IRON));

	std::vector<vec> frame(WIDTH * HEIGHT);
/*
	std::thread thr(&Scene::print_scene, &my_scene, camera, ref(frame));
	if (thr.joinable())
		thr.join();*/
	my_scene.print_scene(camera);
	
	for (int i = 0; i < argc; i++){
		cout<< argv[i] << wctype(argv[i]) << endl;
		if (!strcmp(argv[i], "./a"))
			cout << "eeeeeee" << endl;
	}

	return 0;
}