//
//  Framework for a raytracer
//  File: sphere.cpp
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

#include "sphere.h"
#include <iostream>
#include <math.h>

/************************** Sphere **********************************/

Hit Sphere::intersect(const Ray &ray)
{
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, position, r
    * Sought: intersects? if true: *t
    * 
    * Insert calculation of ray/sphere intersection here. 
    *
    * You have the sphere's center (C) and radius (r) as well as
    * the ray's origin (ray.O) and direction (ray.D).
    *
    * If the ray does not intersect the sphere, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/

	double t;
	
	double a = ray.D.dot(ray.D);
	double b = 2*(ray.D.dot(ray.O-position));
	double c = (ray.O-position).dot(ray.O-position)-(r*r);

	double discrim = ((b*b)-4*a*c);	//discriminant of a quadratic polynomial = b^2 - 4ac

	if (discrim < 0){				//two unequal complex roots
		return Hit::NO_HIT();
	} else {						//two real  roots
		double posT = (-b+sqrt((b*b)-4*a*c))/2*a;
		double negT = (-b-sqrt((b*b)-4*a*c))/2*a;
		
		t = min(posT, negT);

		if (posT < 0 && negT >= 0) {
			t = negT;
		}
		if (negT < 0 && posT >= 0) {
			t = posT;
		}
		if (negT < 0 && posT < 0) {
			return Hit::NO_HIT();
		}
	}

    /****************************************************
    * RT1.2: NORMAL CALCULATION
    *
    * Given: t, C, r
    * Sought: N
    * 
    * Insert calculation of the sphere's normal at the intersection point.
    ****************************************************/
	    
	Point intersect = ray.O + (t*ray.D);

	Vector N = (intersect-position).normalized();
    return Hit(t,N);
}
