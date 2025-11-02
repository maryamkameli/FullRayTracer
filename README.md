## Full Ray Tracer
This project is about rendering complex scenes with multiple geometry types, lighting models, and reflections. The journey involved solving several challenges.
All outputs match the provided reference images except the dragon. 

One of the chalenges I faced was the solid color for some of the challenging renderings. Some scenes (like ShadowTest.txt) would render as a single solid color. To solve this, I implemented complete smooth shading support. 


## Compile
g++ -O3 -std=c++11 rayTrace_vec3.cpp -o raytracer
flag -O3 optimized it to run so faster

## Run
./raytracer <text.file>



## Results
# Bottle
<p align="center">
  <img src="res/bottle.png" width="560">
</p>

# Dragon
<p align="center">
  <img src="res/dragon.png" width="560">
</p>

# Foo realistic
<p align="center">
  <img src="res/foo_r.png" width="560">
</p>

# Foo slow
<p align="center">
  <img src="res/foo-slow.png" width="560">
</p>

# Gear
<p align="center">
  <img src="res/gear.png" width="560">
</p>

# No Lable
<p align="center">
  <img src="res/noLable.png" width="560">
</p>

# Outdoor 
<p align="center">
  <img src="res/outdoor.png" width="560">
</p>

# Plant
<p align="center">
  <img src="res/plant.png" width="560">
</p>

# Reaching hand
<p align="center">
  <img src="res/reachingHand.png" width="560">
</p>

# Arm
<p align="center">
  <img src="res/arm.png" width="560">
</p>

# Shadow Test
<p align="center">
  <img src="res/ShadowTest.png" width="560">
</p>

# Triangle 
<p align="center">
  <img src="res/triangle.png" width="560">
</p>

## Bonus
<p align="center">
  <img src="res/watch_blue_and_gold.png" width="560">
</p>

# Easy mode
<p align="center">
  <img src="res/watch_easy_mode.png" width="560">
</p>


## Implemented Features 
| Feature | Status |
|---------|--------|
| Complex Lighting | Complete | 
| Spot Lights | Complete | 
| HDR Support | Partial |

# Bonus Implementation
Smooth Shading with Per-Vertex Normals

## References 
https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/shading-normals.html
https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html

# Holoween Plan
Spend Holoween at the MIT Visual Arts Center! 

<p align="center">
  <img src="mit.jpeg" width="280">
</p>

