//CSCI 5607 Advanced Ray Tracer
//Features: Jittered/Adaptive Sampling, DoF, Area Lights, Ambient Occlusion, HDR/Tone Mapping

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 
#endif

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>

// Minimal vec3 implementation
struct vec3 {
    float x, y, z;
    vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
    
    float length() const { return sqrt(x*x + y*y + z*z); }
    vec3 normalized() const { float len = length(); return vec3(x/len, y/len, z/len); }
    vec3 clampTo1() const { return vec3(fmin(x,1), fmin(y,1), fmin(z,1)); }
    
    vec3 operator+(const vec3& v) const { return vec3(x+v.x, y+v.y, z+v.z); }
    vec3 operator-(const vec3& v) const { return vec3(x-v.x, y-v.y, z-v.z); }
    vec3 operator*(float f) const { return vec3(x*f, y*f, z*f); }
    vec3 operator/(float f) const { return vec3(x/f, y/f, z/f); }
};

inline vec3 operator*(float f, const vec3& v) { return v * f; }
inline float dot(const vec3& a, const vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline vec3 cross(const vec3& a, const vec3& b) { 
    return vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
}

// Random number generator
std::random_device rd;
std::mt19937 gen(42); // Fixed seed for reproducibility
std::uniform_real_distribution<float> dis(0.0f, 1.0f);

float random01() { return dis(gen); }

// STB Image Write (minimal declaration)
extern "C" int stbi_write_png(char const *filename, int w, int h, int comp, 
                               const void *data, int stride_in_bytes);

const float EPSILON = 0.001f;
const float INF = std::numeric_limits<float>::infinity();

// Configuration
struct RenderConfig {
    bool enableJitteredSampling = true;
    bool enableAdaptiveSampling = true;
    bool enableDepthOfField = false;
    bool enableAmbientOcclusion = true;
    bool enableHDR = true;
    
    int samplesPerPixel = 4;      // For jittered sampling
    int aoSamples = 16;           // Ambient occlusion samples
    float aoRadius = 0.5f;        // AO sampling radius
    float aoStrength = 0.8f;      // AO effect strength
    
    float aperture = 0.0f;        // DoF aperture size (0 = no DoF)
    float focalDistance = 5.0f;   // DoF focal distance
    
    float exposure = 1.0f;        // HDR exposure
    float gamma = 2.2f;           // Gamma correction
    bool enableBloom = true;
    float bloomThreshold = 1.0f;
    float bloomIntensity = 0.3f;
};

RenderConfig config;

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shininess;
    vec3 transmissive;
    float ior;
    vec3 emission;  // For area lights
    
    Material() : ambient(0.1,0.1,0.1), diffuse(0.7,0.7,0.7), 
                 specular(0.3,0.3,0.3), shininess(32),
                 transmissive(0,0,0), ior(1.0), emission(0,0,0) {}
};

struct Sphere {
    vec3 center;
    float radius;
    Material material;
    
    Sphere(vec3 c, float r, Material m) : center(c), radius(r), material(m) {}
};

struct AreaLight {
    vec3 center;
    vec3 u, v;  // Basis vectors defining the light plane
    float width, height;
    vec3 color;
    int samples;
    
    AreaLight(vec3 c, vec3 u, vec3 v, float w, float h, vec3 col, int s = 16)
        : center(c), u(u), v(v), width(w), height(h), color(col), samples(s) {}
    
    vec3 samplePoint() const {
        float su = (random01() - 0.5f) * width;
        float sv = (random01() - 0.5f) * height;
        return center + su * u + sv * v;
    }
};

struct Ray {
    vec3 origin;
    vec3 direction;
    Ray(vec3 o, vec3 d) : origin(o), direction(d.normalized()) {}
};

struct HitInfo {
    bool hit;
    float t;
    vec3 point;
    vec3 normal;
    Material material;
    
    HitInfo() : hit(false), t(INF) {}
};

// Scene data
std::vector<Sphere> spheres;
std::vector<AreaLight> areaLights;
vec3 backgroundColor(0.2f, 0.3f, 0.5f);
vec3 ambientLight(0.1f, 0.1f, 0.1f);

// Camera
vec3 eye(0, 0, 0);
vec3 forward(0, 0, -1);
vec3 up(0, 1, 0);
vec3 right(-1, 0, 0);
float halfAngleVFOV = 35;

// Image
int img_width = 800;
int img_height = 600;

// ============ Intersection Functions ============

HitInfo raySphereIntersect(const Ray& ray, const Sphere& sphere) {
    HitInfo info;
    vec3 oc = ray.origin - sphere.center;
    
    float a = dot(ray.direction, ray.direction);
    float b = 2.0f * dot(ray.direction, oc);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) return info;
    
    float t0 = (-b - sqrt(discriminant))/(2.0f*a);
    float t1 = (-b + sqrt(discriminant))/(2.0f*a);
    
    float t = t0;
    if (t < EPSILON) {
        t = t1;
        if (t < EPSILON) return info;
    }
    
    info.hit = true;
    info.t = t;
    info.point = ray.origin + t * ray.direction;
    info.normal = (info.point - sphere.center).normalized();
    info.material = sphere.material;
    
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
    return closest;
}

bool inShadow(const vec3& point, const vec3& lightDir, float lightDist) {
    Ray shadowRay(point + EPSILON * lightDir, lightDir);
    HitInfo shadowHit = intersectScene(shadowRay);
    return shadowHit.hit && shadowHit.t < lightDist - EPSILON;
}

// ============ Ambient Occlusion ============

vec3 randomHemisphereDirection(const vec3& normal) {
    // Cosine-weighted hemisphere sampling
    float u1 = random01();
    float u2 = random01();
    
    float r = sqrt(u1);
    float theta = 2.0f * M_PI * u2;
    
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(1.0f - u1);
    
    // Build orthonormal basis
    vec3 w = normal;
    vec3 u = (fabs(w.x) > 0.1f ? vec3(0,1,0) : vec3(1,0,0));
    u = cross(u, w).normalized();
    vec3 v = cross(w, u);
    
    return (x * u + y * v + z * w).normalized();
}

float computeAmbientOcclusion(const vec3& point, const vec3& normal) {
    if (!config.enableAmbientOcclusion) return 1.0f;
    
    int hitCount = 0;
    for (int i = 0; i < config.aoSamples; i++) {
        vec3 dir = randomHemisphereDirection(normal);
        Ray aoRay(point + EPSILON * normal, dir);
        HitInfo hit = intersectScene(aoRay);
        
        if (hit.hit && hit.t < config.aoRadius) {
            hitCount++;
        }
    }
    
    float occlusion = 1.0f - (float(hitCount) / config.aoSamples);
    return pow(occlusion, config.aoStrength);
}

// ============ Lighting ============

vec3 computeDirectLighting(const HitInfo& hit, const vec3& viewDir) {
    vec3 color = hit.material.ambient * ambientLight;
    
    // Area lights with soft shadows
    for (const auto& light : areaLights) {
        vec3 lightContribution(0, 0, 0);
        int visibleSamples = 0;
        
        for (int i = 0; i < light.samples; i++) {
            vec3 lightPoint = light.samplePoint();
            vec3 lightDir = lightPoint - hit.point;
            float lightDist = lightDir.length();
            lightDir = lightDir.normalized();
            
            if (!inShadow(hit.point, lightDir, lightDist)) {
                visibleSamples++;
                
                // Diffuse
                float diff = std::max(0.0f, dot(hit.normal, lightDir));
                vec3 diffuse = vec3(
                    hit.material.diffuse.x * light.color.x * diff,
                    hit.material.diffuse.y * light.color.y * diff,
                    hit.material.diffuse.z * light.color.z * diff
                );
                
                // Specular
                vec3 reflectDir = 2.0f * dot(lightDir, hit.normal) * hit.normal - lightDir;
                float spec = pow(std::max(0.0f, dot(reflectDir, viewDir)), hit.material.shininess);
                vec3 specular = vec3(
                    hit.material.specular.x * light.color.x * spec,
                    hit.material.specular.y * light.color.y * spec,
                    hit.material.specular.z * light.color.z * spec
                );
                
                float attenuation = 1.0f / (lightDist * lightDist);
                lightContribution = lightContribution + attenuation * (diffuse + specular);
            }
        }
        
        // Average over all samples (soft shadows)
        color = color + lightContribution * (1.0f / light.samples);
    }
    
    return color;
}

vec3 shade(const HitInfo& hit, const Ray& ray, int depth) {
    if (depth <= 0) return vec3(0,0,0);
    
    vec3 viewDir = ray.direction * -1.0f;
    
    // Direct lighting
    vec3 color = computeDirectLighting(hit, viewDir);
    
    // Ambient occlusion
    float ao = computeAmbientOcclusion(hit.point, hit.normal);
    color = color * ao;
    
    // Add emission (for area lights on spheres)
    color = color + hit.material.emission;
    
    // Simple reflection
    float reflectivity = (hit.material.specular.x + hit.material.specular.y + 
                          hit.material.specular.z) / 3.0f;
    if (reflectivity > 0.01f && depth > 0) {
        vec3 reflectDir = ray.direction - 2.0f * dot(ray.direction, hit.normal) * hit.normal;
        Ray reflectRay(hit.point + EPSILON * hit.normal, reflectDir);
        HitInfo reflectHit = intersectScene(reflectRay);
        
        if (reflectHit.hit) {
            vec3 reflectColor = shade(reflectHit, reflectRay, depth - 1);
            color = color + reflectivity * reflectColor;
        }
    }
    
    return color;
}

vec3 traceRay(const Ray& ray, int depth = 5) {
    HitInfo hit = intersectScene(ray);
    if (hit.hit) {
        return shade(hit, ray, depth);
    }
    return backgroundColor;
}

// ============ HDR & Tone Mapping ============

vec3 toneMap(const vec3& hdrColor) {
    if (!config.enableHDR) return hdrColor.clampTo1();
    
    // Reinhard tone mapping
    vec3 exposed = hdrColor * config.exposure;
    vec3 mapped(
        exposed.x / (1.0f + exposed.x),
        exposed.y / (1.0f + exposed.y),
        exposed.z / (1.0f + exposed.z)
    );
    
    // Gamma correction
    mapped = vec3(
        pow(mapped.x, 1.0f / config.gamma),
        pow(mapped.y, 1.0f / config.gamma),
        pow(mapped.z, 1.0f / config.gamma)
    );
    
    return mapped.clampTo1();
}

// ============ Sampling Strategies ============

vec3 jitteredSample(int i, int j, int sample, int totalSamples) {
    int sqrtSamples = (int)sqrt(totalSamples);
    int sx = sample % sqrtSamples;
    int sy = sample / sqrtSamples;
    
    float jitterX = (sx + random01()) / sqrtSamples;
    float jitterY = (sy + random01()) / sqrtSamples;
    
    return vec3(i + jitterX, j + jitterY, 0);
}

float colorVariance(const std::vector<vec3>& samples) {
    if (samples.size() < 2) return 0.0f;
    
    vec3 mean(0, 0, 0);
    for (const auto& s : samples) mean = mean + s;
    mean = mean / float(samples.size());
    
    float variance = 0.0f;
    for (const auto& s : samples) {
        vec3 diff = s - mean;
        variance += dot(diff, diff);
    }
    return variance / samples.size();
}

// ============ Depth of Field ============

Ray generateDOFRay(const vec3& pixelPoint, const vec3& focalPoint) {
    if (config.aperture <= 0.0f) {
        return Ray(eye, (pixelPoint - eye).normalized());
    }
    
    // Random point on aperture disk
    float angle = 2.0f * M_PI * random01();
    float radius = config.aperture * sqrt(random01());
    vec3 offset = right * (radius * cos(angle)) + up * (radius * sin(angle));
    
    vec3 newOrigin = eye + offset;
    vec3 newDir = (focalPoint - newOrigin).normalized();
    
    return Ray(newOrigin, newDir);
}

// ============ Main Rendering Loop ============

void setupScene() {
    // Create materials
    Material redMaterial;
    redMaterial.diffuse = vec3(0.8f, 0.2f, 0.2f);
    redMaterial.specular = vec3(0.3f, 0.3f, 0.3f);
    redMaterial.shininess = 32;
    
    Material greenMaterial;
    greenMaterial.diffuse = vec3(0.2f, 0.8f, 0.2f);
    greenMaterial.specular = vec3(0.3f, 0.3f, 0.3f);
    greenMaterial.shininess = 32;
    
    Material blueMaterial;
    blueMaterial.diffuse = vec3(0.2f, 0.2f, 0.8f);
    blueMaterial.specular = vec3(0.6f, 0.6f, 0.6f);
    blueMaterial.shininess = 64;
    
    Material groundMaterial;
    groundMaterial.diffuse = vec3(0.5f, 0.5f, 0.5f);
    groundMaterial.specular = vec3(0.1f, 0.1f, 0.1f);
    groundMaterial.shininess = 16;
    
    // Add spheres
    spheres.push_back(Sphere(vec3(-2, 0, -6), 1.0f, redMaterial));
    spheres.push_back(Sphere(vec3(0, 0, -6), 1.0f, greenMaterial));
    spheres.push_back(Sphere(vec3(2, 0, -6), 1.0f, blueMaterial));
    spheres.push_back(Sphere(vec3(0, -101, -6), 100.0f, groundMaterial)); // Ground
    
    // Add area light
    vec3 lightCenter(0, 4, -5);
    vec3 lightU = right;
    vec3 lightV = forward;
    areaLights.push_back(AreaLight(lightCenter, lightU, lightV, 2.0f, 2.0f, 
                                    vec3(1.5f, 1.5f, 1.5f), 16));
}

int main(int argc, char** argv) {
    std::cout << "=== Advanced Ray Tracer ===" << std::endl;
    std::cout << "Features: Jittered/Adaptive Sampling, DoF, Area Lights, AO, HDR" << std::endl;
    
    // Setup scene
    setupScene();
    
    // Setup camera
    eye = vec3(0, 1, 2);
    forward = vec3(0, 0, -1).normalized();
    up = vec3(0, 1, 0).normalized();
    right = cross(forward, up).normalized();
    up = cross(right, forward).normalized();
    
    // Configure rendering (can be set via command line or config file)
    config.enableJitteredSampling = true;
    config.enableAdaptiveSampling = false; // Set true for adaptive
    config.enableDepthOfField = false;     // Set true for DoF
    config.enableAmbientOcclusion = true;
    config.enableHDR = true;
    config.samplesPerPixel = 4;
    config.aoSamples = 16;
    config.aoRadius = 0.5f;
    
    std::vector<unsigned char> pixels(img_width * img_height * 3);
    std::vector<vec3> hdrBuffer(img_width * img_height);
    
    float imgW = img_width, imgH = img_height;
    float halfW = imgW/2, halfH = imgH/2;
    float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    int totalPixels = img_width * img_height;
    int processedPixels = 0;
    
    for (int j = 0; j < img_height; j++) {
        for (int i = 0; i < img_width; i++) {
            vec3 pixelColor(0, 0, 0);
            std::vector<vec3> samples;
            
            int numSamples = config.samplesPerPixel;
            
            for (int s = 0; s < numSamples; s++) {
                float u, v;
                
                if (config.enableJitteredSampling) {
                    vec3 jittered = jitteredSample(i, j, s, numSamples);
                    u = halfW - imgW * (jittered.x / imgW);
                    v = halfH - imgH * (jittered.y / imgH);
                } else {
                    u = halfW - imgW * ((i + 0.5f) / imgW);
                    v = halfH - imgH * ((j + 0.5f) / imgH);
                }
                
                vec3 pixelPoint = eye + d * forward * -1.0f + u * right + v * up;
                
                Ray ray(eye, (pixelPoint - eye).normalized());
                
                // Apply depth of field if enabled
                if (config.enableDepthOfField && config.aperture > 0) {
                    vec3 focalPoint = eye + config.focalDistance * ray.direction;
                    ray = generateDOFRay(pixelPoint, focalPoint);
                }
                
                vec3 sampleColor = traceRay(ray);
                samples.push_back(sampleColor);
                pixelColor = pixelColor + sampleColor;
            }
            
            // Adaptive sampling
            if (config.enableAdaptiveSampling && numSamples < 16) {
                float variance = colorVariance(samples);
                if (variance > 0.1f) { // High variance, need more samples
                    for (int s = 0; s < 4; s++) {
                        vec3 jittered = jitteredSample(i, j, numSamples + s, numSamples + 4);
                        float u = halfW - imgW * (jittered.x / imgW);
                        float v = halfH - imgH * (jittered.y / imgH);
                        vec3 pixelPoint = eye + d * forward * -1.0f + u * right + v * up;
                        Ray ray(eye, (pixelPoint - eye).normalized());
                        vec3 sampleColor = traceRay(ray);
                        pixelColor = pixelColor + sampleColor;
                        numSamples++;
                    }
                }
            }
            
            pixelColor = pixelColor / float(numSamples);
            hdrBuffer[j * img_width + i] = pixelColor;
            
            // Tone mapping
            vec3 finalColor = toneMap(pixelColor);
            
            int idx = (j * img_width + i) * 3;
            pixels[idx + 0] = (unsigned char)(finalColor.x * 255);
            pixels[idx + 1] = (unsigned char)(finalColor.y * 255);
            pixels[idx + 2] = (unsigned char)(finalColor.z * 255);
            
            processedPixels++;
            if (processedPixels % (totalPixels / 10) == 0) {
                std::cout << "Progress: " << (100 * processedPixels / totalPixels) << "%" << std::endl;
            }
        }
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    printf("Rendering took %.2f ms\n", 
           std::chrono::duration<double, std::milli>(t_end - t_start).count());
    
    std::string filename = "advanced_raytraced.png";
    std::cout << "Writing image to " << filename << std::endl;
    stbi_write_png(filename.c_str(), img_width, img_height, 3, 
                   pixels.data(), img_width * 3);
    
    std::cout << "\nFeatures used:" << std::endl;
    std::cout << "- Jittered Sampling: " << (config.enableJitteredSampling ? "ON" : "OFF") << std::endl;
    std::cout << "- Adaptive Sampling: " << (config.enableAdaptiveSampling ? "ON" : "OFF") << std::endl;
    std::cout << "- Depth of Field: " << (config.enableDepthOfField ? "ON" : "OFF") << std::endl;
    std::cout << "- Ambient Occlusion: " << (config.enableAmbientOcclusion ? "ON" : "OFF") << std::endl;
    std::cout << "- HDR/Tone Mapping: " << (config.enableHDR ? "ON" : "OFF") << std::endl;
    std::cout << "- Area Lights/Soft Shadows: ON" << std::endl;
    
    return 0;
}