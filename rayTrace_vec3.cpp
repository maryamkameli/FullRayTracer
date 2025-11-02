//CSCI 5607 HW3 - Rays & Files
//This HW has three steps:
// 1. Compile and run the program (the program takes a single command line argument)
// 2. Understand the code in this file (rayTrace_vec3.cpp), in particular be sure to understand:
//     -How ray-sphere intersection works
//     -How the rays are being generated
//     -The pipeline from rays, to intersection, to pixel color
//    After you finish this step, and understand the math, take the HW quiz on canvas
// 3. Update the file parse_vec3.h so that the function parseSceneFile() reads the passed in file
//     and sets the relevant global variables for the rest of the code to product to correct image

//To Compile: g++ -fsanitize=address -std=c++11 rayTrace_vec3.cpp

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // For fopen and sscanf
#define _USE_MATH_DEFINES 
#endif

//Images Lib includes:
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "image_lib.h" //Defines an image class and a color class

//#Vec3 Library
#include "vec3.h"

//High resolution timer
#include <chrono>
#include <cmath>
#include <algorithm>
#include <limits>

//Scene file parser
#include "parse_vec3.h"

const float EPSILON = 0.0001f;
const float INF = std::numeric_limits<float>::infinity();

struct Ray {
    vec3 origin;
    vec3 direction;
    
    Ray(vec3 o, vec3 d) : origin(o), direction(d) {}
};

struct HitInfo {
    bool hit;
    float t;
    vec3 point;
    vec3 normal;
    Material material;
    
    HitInfo() : hit(false), t(INF) {}
};

HitInfo raySphereIntersect(const Ray& ray, const Sphere& sphere) {
    HitInfo info;
    vec3 oc = ray.origin - sphere.center;
    
    float a = dot(ray.direction, ray.direction);
    float b = 2.0f * dot(ray.direction, oc);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) {
        return info;
    }
    
    float t0 = (-b - sqrt(discriminant))/(2.0f*a);
    float t1 = (-b + sqrt(discriminant))/(2.0f*a);
    
    float t = t0;
    if (t < EPSILON) {
        t = t1;
        if (t < EPSILON) {
            return info;
        }
    }
    
    info.hit = true;
    info.t = t;
    info.point = ray.origin + t * ray.direction;
    info.normal = (info.point - sphere.center).normalized();
    info.material = sphere.material;
    
    return info;
}

// Ray-Triangle Intersection
HitInfo rayTriangleIntersect(const Ray& ray, const Triangle& tri) {
    HitInfo info;
    
    vec3 v0 = vertices[tri.v0];
    vec3 v1 = vertices[tri.v1];
    vec3 v2 = vertices[tri.v2];
    
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    
    vec3 h = cross(ray.direction, edge2);
    float a = dot(edge1, h);
    
    // Ray is parallel to triangle
    if (a > -EPSILON && a < EPSILON) {
        return info;
    }
    
    float f = 1.0f / a;
    vec3 s = ray.origin - v0;
    float u = f * dot(s, h);
    
    if (u < 0.0f || u > 1.0f) {
        return info;
    }
    
    vec3 q = cross(s, edge1);
    float v = f * dot(ray.direction, q);
    
    if (v < 0.0f || u + v > 1.0f) {
        return info;
    }
    
    float t = f * dot(edge2, q);
    
    if (t > EPSILON) {
        info.hit = true;
        info.t = t;
        info.point = ray.origin + t * ray.direction;
        
        // smooth shading
        if (tri.n0 >= 0 && tri.n1 >= 0 && tri.n2 >= 0) {
            // Interpolate normals using barycentric coordinates
            float w = 1.0f - u - v;  // Barycentric coordinate for v0
            vec3 n0 = normals[tri.n0];
            vec3 n1 = normals[tri.n1];
            vec3 n2 = normals[tri.n2];
            info.normal = (w * n0 + u * n1 + v * n2).normalized();
        } else {
            info.normal = cross(edge1, edge2).normalized();
            
            if (dot(info.normal, ray.direction) > 0) {
                info.normal = -1.0f * info.normal;
            }
        }
        
        info.material = tri.material;
    }
    
    return info;
}

HitInfo intersectScene(const Ray& ray) {
    HitInfo closest;
    
    for (const auto& sphere : spheres) {
        HitInfo info = raySphereIntersect(ray, sphere);
        if (info.hit && info.t < closest.t) {
            closest = info;
        }
    }
    
    for (const auto& triangle : triangles) {
        HitInfo info = rayTriangleIntersect(ray, triangle);
        if (info.hit && info.t < closest.t) {
            closest = info;
        }
    }
    
    return closest;
}

bool inShadow(const vec3& point, const vec3& lightDir, float lightDist) {
    Ray shadowRay(point + EPSILON * lightDir, lightDir);
    HitInfo shadowHit = intersectScene(shadowRay);
    return shadowHit.hit && shadowHit.t < lightDist - EPSILON;
}

vec3 shade(const HitInfo& hit, const Ray& ray, int depth);

// Phong shading
vec3 computePhong(const HitInfo& hit, const vec3& viewDir) {
    // Ambient component
    vec3 color = vec3(
        hit.material.ambient.x * ambientLight.x,
        hit.material.ambient.y * ambientLight.y,
        hit.material.ambient.z * ambientLight.z
    );

    // Point lights
    for (const auto& light : pointLights) {
        vec3 lightDir = light.position - hit.point;
        float lightDist = lightDir.length();
        lightDir = lightDir.normalized();
        
        if (!inShadow(hit.point, lightDir, lightDist)) {
            // Diffuse
            float diff = std::max(0.0f, dot(hit.normal, lightDir));
            vec3 diffuse = vec3(
                hit.material.diffuse.x * light.color.x * diff,
                hit.material.diffuse.y * light.color.y * diff,
                hit.material.diffuse.z * light.color.z * diff
            );            
            float attenuation = 1.0f / (lightDist * lightDist);
            color = color + attenuation * diffuse;
        }
    }
    
    // Directional lights
    for (const auto& light : directionalLights) {
        vec3 lightDir = -1.0f * light.direction;
        
        if (!inShadow(hit.point, lightDir, INF)) {
            // Diffuse
            float diff = std::max(0.0f, dot(hit.normal, lightDir));
            vec3 diffuse = vec3(
                hit.material.diffuse.x * light.color.x * diff,
                hit.material.diffuse.y * light.color.y * diff,
                hit.material.diffuse.z * light.color.z * diff
            );            
            color = color + diffuse;
        }
    }
    
    // Spot lights
    for (const auto& light : spotLights) {
        vec3 lightDir = light.position - hit.point;
        float lightDist = lightDir.length();
        lightDir = lightDir.normalized();
        
        float theta = acos(dot(-1.0f * lightDir, light.direction)) * 180.0f / M_PI;
        
        if (theta <= light.angle2) {
            float intensity = 1.0f;
            if (theta > light.angle1) {
                intensity = (light.angle2 - theta) / (light.angle2 - light.angle1);
            }
            
            if (!inShadow(hit.point, lightDir, lightDist)) {
                // Diffuse
                float diff = std::max(0.0f, dot(hit.normal, lightDir));
                vec3 diffuse = intensity * vec3(
                    hit.material.diffuse.x * light.color.x * diff,
                    hit.material.diffuse.y * light.color.y * diff,
                    hit.material.diffuse.z * light.color.z * diff
                );                
                float attenuation = 1.0f / (lightDist * lightDist);
                color = color + attenuation * diffuse;
            }
        }
    }
    
    return color;
}

vec3 shade(const HitInfo& hit, const Ray& ray, int depth) {
    if (depth <= 0) return vec3(0,0,0);
    
    vec3 viewDir = -1.0f * ray.direction;
    vec3 color = computePhong(hit, viewDir);
    
    // Reflection 
    vec3 reflectCoeff = vec3(
        hit.material.transmissive.x + hit.material.specular.x,
        hit.material.transmissive.y + hit.material.specular.y,
        hit.material.transmissive.z + hit.material.specular.z
    );
    float reflectivity = reflectCoeff.x + reflectCoeff.y + reflectCoeff.z;
    
    if (reflectivity > 0.01f && depth > 0) {
        vec3 reflectDir = ray.direction - 2.0f * dot(ray.direction, hit.normal) * hit.normal;
        Ray reflectRay(hit.point + EPSILON * hit.normal, reflectDir);
        HitInfo reflectHit = intersectScene(reflectRay);
        
        if (reflectHit.hit) {
            vec3 reflectColor = shade(reflectHit, reflectRay, depth - 1);
            color = color + vec3(
                reflectCoeff.x * reflectColor.x,
                reflectCoeff.y * reflectColor.y,
                reflectCoeff.z * reflectColor.z
            );
        } else {
            // Reflect background
            color = color + vec3(
                reflectCoeff.x * backgroundColor.x,
                reflectCoeff.y * backgroundColor.y,
                reflectCoeff.z * backgroundColor.z
            );
        }
    }
    
    return color;
}

vec3 traceRay(const Ray& ray, int depth) {
    HitInfo hit = intersectScene(ray);
    
    if (hit.hit) {
        return shade(hit, ray, depth);
    }
    
    return backgroundColor;
}

int main(int argc, char** argv){
    if (argc != 2){
        std::cout << "Usage: " << argv[0] << " <scenefile>\n";
        std::cout << "Example: " << argv[0] << " triangle.txt\n";
        return 0;
    }
    std::string sceneFileName = argv[1];

    parseSceneFile(sceneFileName);

    float imgW = img_width, imgH = img_height;
    float halfW = imgW/2, halfH = imgH/2;
    float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));

    Image outputImg = Image(img_width, img_height);
    auto t_start = std::chrono::high_resolution_clock::now();
    
    printf("Rendering %dx%d image...\n", img_width, img_height);
    
    for (int i = 0; i < img_width; i++){
        if (i % 100 == 0) {
            printf("Progress: %d/%d\n", i, img_width);
        }
        
        for (int j = 0; j < img_height; j++){
            float u = -1.0f * (halfW - (imgW)*((i+0.5)/imgW));
            float v = (halfH - (imgH)*((j+0.5)/imgH));
            vec3 p = eye - d*forward + u*right + v*up;
            vec3 rayDir = (p - eye).normalized();

            Ray ray(eye, rayDir);
            vec3 color = traceRay(ray, maxDepth);
            
            color = color.clampTo1();

            outputImg.setPixel(i, j, Color(color.x, color.y, color.z));
        }
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    printf("Rendering took %.2f ms\n",
           std::chrono::duration<double, std::milli>(t_end-t_start).count());
        outputImg.write(imgName.c_str());
    
    return 0;
}