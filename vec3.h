#ifndef VEC3_H
#define VEC3_H

#include <cmath>

//Small vector library
// Represents a vector as 3 floats

struct vec3{
  float x,y,z;

  vec3(float x, float y, float z) : x(x), y(y), z(z) {}
  vec3() : x(0), y(0), z(0) {}

  //Clamp each component (used to clamp pixel colors)
  vec3 clampTo1(){
    return vec3(fmin(fmax(x,0),1), fmin(fmax(y,0),1), fmin(fmax(z,0),1));
  }

  //Compute vector length
  float length(){
    return sqrt(x*x+y*y+z*z);
  }

  //Create a unit-length vector
  vec3 normalized(){
    float len = sqrt(x*x+y*y+z*z);
    if (len < 0.0001f) return vec3(0,0,0);
    return vec3(x/len,y/len,z/len);
  }

};

//Multiply float and vector
inline vec3 operator*(float f, vec3 a){
  return vec3(a.x*f,a.y*f,a.z*f);
}

//Multiply vector and float
inline vec3 operator*(vec3 a, float f){
  return vec3(a.x*f,a.y*f,a.z*f);
}

//Component-wise multiplication (for materials and colors)
inline vec3 operator*(vec3 a, vec3 b){
  return vec3(a.x*b.x, a.y*b.y, a.z*b.z);
}

//Vector-vector dot product
inline float dot(vec3 a, vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

//Vector-vector cross product
inline vec3 cross(vec3 a, vec3 b){
  return vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y-a.y*b.x);
}

//Vector addition
inline vec3 operator+(vec3 a, vec3 b){
  return vec3(a.x+b.x, a.y+b.y, a.z+b.z);
}

//Vector subtraction
inline vec3 operator-(vec3 a, vec3 b){
  return vec3(a.x-b.x, a.y-b.y, a.z-b.z);
}

//Unary negation
inline vec3 operator-(vec3 a){
  return vec3(-a.x, -a.y, -a.z);
}

//Vector division by scalar
inline vec3 operator/(vec3 a, float f){
  return vec3(a.x/f, a.y/f, a.z/f);
}

#endif