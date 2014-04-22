//
//  Framework for a raytracer
//  File: scene.cpp
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Maarten Everts
//    Jasper van de Gronde
//
//  This framework is inspired by and uses code of the raytracer framework of 
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html 
//

#include "scene.h"
#include "material.h"
#include "triple.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

bool Shadows = true; // toggle shadows on/off

Color Scene::trace(const Ray &ray, int max_recurs)
{
    // Find hit object and distance
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    Object *obj = NULL;
    for (unsigned int i = 0; i < objects.size(); ++i) {
        Hit hit(objects[i]->intersect(ray));
        if (hit.t<min_hit.t) {
            min_hit = hit;
            obj = objects[i];
        }
    }

    // No hit? Return background color.
    if (!obj) return Color(0.0, 0.0, 0.0);

    Material *material = obj->material;            //the hit objects material
    Point hit = ray.at(min_hit.t);                 //the hit point
    Vector N = min_hit.N;                          //the normal at hit point
    Vector V = -ray.D;                             //the view vector


    /****************************************************
    * This is where you should insert the color
    * calculation (Phong model).
    *
    * Given: material, hit, N, V, lights[]
    * Sought: color
    *
    * Hints: (see triple.h)
    *        Triple.dot(Vector) dot product
    *        Vector+Vector      vector sum
    *        Vector-Vector      vector difference
    *        Point-Point        yields vector
    *        Vector.normalize() normalizes vector, returns length
    *        double*Color        scales each color component (r,g,b)
    *        Color*Color        dito
    *        pow(a,b)           a to the power of b
    ****************************************************/

	// Colour variables                 
	Color id, ia, is;
	Color cl, cr, cp; 
	cp = Color(1.0, 1.0, 1.0);	//cp is optional

	Color color;
	Color reflection;
	Color refraction;

	int max_shadow_ray = 15; // set the number of rays used for soft shadows

	for (int i = 0; i < lights.size(); i++) {
		// Phong's illumination Model Calculations
		Vector L = (lights[i]->position - hit).normalized();
		Vector R = (-L + 2*(L.dot(N))*N).normalized();			

		if (Shadows == true){
			Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
			Object *obj = NULL;

			double shade = 0;
			bool flag = true;

			// generate a number of rays that will strike the light source in random positions
			for (int numRays = 0; numRays < max_shadow_ray; numRays++) { 
						
				//calculate random point on sphere
				double rand2 = ((double)rand())/(RAND_MAX) * (2*3.14158265); //random value 0 to 2pi
				
				double r = 45;

				double x = r * cos(rand2);
				double y = r * sin(rand2);
				//double z = r * cos(rand2);
								
				//add the coordinate values to the original position of the light sphere
				Point randomPoint = Point(lights[i]->position.x + x, lights[i]->position.y + y, lights[i]->position.z);

				Vector shadowL = (randomPoint-hit).normalized();
			
				Ray shadow(hit, shadowL);
				Point hitj = shadow.at(pow(2.0,-32)); //jiggle
				shadow = Ray(hitj, shadowL);

				for (unsigned int j = 0; j < objects.size(); ++j) {
					Hit feeler(objects[j]->intersect(shadow));
					if (feeler.t<min_hit.t) {
						flag = false;
						break;
					}
				}	
				if (flag) {
					//add to average the values returned by the rays to determine a value of the shadow (number between 0 and 1)
					shade += 1;
				}				
			}
			shade = shade / max_shadow_ray; //divide to average the values

			//averaged values are used to weight the diffuse and specular components of the objects colour
			if (!obj) { // Shadow ray doesn't intersect another object
				cl = lights[i]->color;
				cr = material->color;
				id = cl*cr * max(0.0,(L.dot(N)))*(material->kd)* shade;
				ia = cl*cr * material->ka;
				is = cl * pow(max(0.0, (V.normalized().dot(R))), material->n) * material->ks * cp* shade;
				color += id + ia + is;
			}else{ //Shadow ray intersects an object
				cl = lights[i]->color;
				ia = cl * material->ka * material->color;
				color += ia;
			}			
		} else { // Shadows are toggled off
			cl = lights[i]->color;
			id = (cl*material->color) * max(0.0,(L.dot(N)))*(material->kd);
			ia = cl * material->ka * material->color;
			is = cl * pow(max(0.0, (V.normalized().dot(R))), material->n)*material->ks;
			color += id + ia + is;
		}
	}
	
	// reflection calculations and recursion
	if ((max_recurs > 0)&&((material->reflect) > 0)) {
		Vector ref = (ray.D-2*(ray.D.dot(N)*N)).normalized(); // direction of reflection

		Ray reflected_ray = Ray(hit, ref); // reflection ray
		Point hit_jiggle = reflected_ray.at(pow(2.0,-32));
		reflected_ray = Ray(hit_jiggle, ref);

		reflection = trace(reflected_ray, (max_recurs-1)) * material->reflect;
	}	
	
	// refraction calculations and recursion
	if ((max_recurs > 0)&&((material->refract) > 0)) {		
		double n = 1.00;
		double nt = material->eta;	

		Vector negV = -V;
		
		double discrim = (1.0-(pow(negV.dot(N),2.0)*pow(n,2.0))/(pow(nt,2.0)));
		Vector transmission = (n*(negV-N*(negV.dot(N)))/nt) - N*(sqrt(discrim));

		transmission.normalize();
		
		Ray refracted_ray = Ray(hit, transmission);
		Point hit_jiggle = refracted_ray.at(pow(2.0,-32));
		refracted_ray = Ray(hit_jiggle, transmission);

		Hit refract(obj->intersect(refracted_ray));

		if(refract.t) {
			Vector exit_N = -refract.N;
			double exit_discrim = (1.0 - (pow(nt,2.0)*(1-pow(transmission.dot(exit_N),2.0))))/(pow(n, 2.0));
			Vector exit_transmission = (nt*(transmission-exit_N*(transmission.dot(exit_N))))/(n)-exit_N*sqrt(exit_discrim);
			exit_transmission.normalize();

			refracted_ray = Ray(refracted_ray.at(refract.t), exit_transmission);
			Point hit_jiggle = refracted_ray.at(pow(2.0,-32));
			refracted_ray = Ray(hit_jiggle, exit_transmission);

			refraction = trace(refracted_ray, max_recurs-1)* material->refract;
		}
	}
	color += refraction;
	color += reflection;
    return color;
}

void Scene::render(Image &img)
{
    int w = img.width();
    int h = img.height();
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            Point pixel(x, h-1-y, 0);
            Ray ray(eye, (pixel-eye).normalized());
            Color col = trace(ray,5);
            col.clamp();
            img(x,y) = col;
        }
    }
}

void Scene::addObject(Object *o)
{
    objects.push_back(o);
}

void Scene::addLight(Light *l)
{
    lights.push_back(l);
}

void Scene::setEye(Triple e)
{
    eye = e;
}
