#include <cassert>
#include <iostream>
#include <cmath>

using namespace std;


struct vec {
	float x;
	float y;
	float z;

	vec(float _x = 0., float _y = 0., float _z = 0.) : x(_x), y(_y), z(_z) {}
	
	float norm() {
		return sqrt(x * x + y * y + z * z);
	}

	vec normalize() {
		float _n = (*this).norm();
		x /= _n;
		y /= _n;
		z /= _n;
		return (*this);
	}

	friend vec operator+(const vec& first, const vec& second);
	friend vec operator-(const vec& first, const vec& second);
	friend float operator*(const vec& first, const vec& second);
	friend vec operator*(const float& first, const vec& second);
	friend vec operator*(const float& first, const vec& second);
	friend vec operator*(const vec& first, const float& second);
	friend ostream& operator<<(ostream& out, const vec& first);

	float& operator[](const int i) {
		assert (i < 3);
		return i == 0 ? x : (i ==1 ? y : z);
	}
};

vec operator+(const vec& first, const vec& second) {
	return vec(first.x + second.x, first.y + second.y, first.z + second.z);
}

vec operator-(const vec& first, const vec& second) {
	return vec(first.x - second.x, first.y - second.y, first.z - second.z);
}

float operator*(const vec& first, const vec& second) {
	return (first.x * second.x + first.y * second.y + first.z * second.z);
}

vec operator*(const float& first, const vec& second) {
	return vec(first * second.x, first * second.y, first * second.z);
}

vec operator*(const vec& first, const float& second) {
	return vec(first.x * second, first.y * second, first.z * second);
}

ostream& operator<<(ostream& out, const vec& first) {
	out << first.x << ' ' << first.y << ' ' << first.z << endl;
	return out;
}