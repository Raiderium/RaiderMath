module raider.math.misc;

import std.math;
import raider.math.vec;

T clamp(T)(T f, T min, T max)
{
	if(f>max) f = max;
	if(f<min) f = min;
	return f;
}

///Next higher power of two
size_t nhpo2(size_t x)
{
	size_t p = 1;
	while(p <= x) p <<= 1;
	return p;
}

///Is power of two
bool ispo2(size_t x)
{
	return x && !(x & (x - 1));
}

///Linear interpolation
T lerp(T)(double t, T a, T b)
{
	return a + (b - a)*t;
}

///Spherical linear interpolation
vec3 slerp(double t, vec3 a, vec3 b)
{
	double dot = clamp(a.dot(b), -1.0, 1.0);
	double theta = acos(dot)*t;
	vec3 r = b - a*dot;
	return a*cast(double)cos(theta) + r*cast(double)sin(theta);
}

///Normalised linear interpolation
vec3 nlerp(double t, vec3 a, vec3 b)
{
	vec3 result = lerp(t, a, b);
	result.normalize();
	return result;
}

///Cubic bezier interpolation
T bez(T)(double t, T p0, T p1, T p2, T p3)
{
	//De Casteljau's algorithm
	T l0, l1, l2;
	
	l0 = lerp(t, p0, p1);
	l1 = lerp(t, p1, p2);
	l2 = lerp(t, p2, p3);
	
	l0 = lerp(t, l0, l1);
	l1 = lerp(t, l1, l2);
	
	l0 = lerp(t, l0, l1);
	
	return l0;

	/* For the sake of curiosity, this is the polynomial form
	 * A really curious person might profile this against de casteljau

	T c0, c1, c2, c3;
	
	c0 = p0;
	c1 = 3*(p1 - p0);
	c2 = 3*(p0 - 2*p1 + p2);
	c3 = p3 - p0 + 3*(p1 - p2);

	T r = c0;
	r+=t*c1; t*=t;
	r+=t*c2; t*=t;
	r+=t*c3;
	return r;
	*/
}

///Use forward differencing to evaluate multiple points on a cubic bezier.
void bez(T)(ref T[] result, T p0, T p1, T p2, T p3)
{
	double t = 1.0 / result.length;
	double t2 = t*t;
	double t3 = t2*t;

	//Find polynomial coefficients
	T c1 = 3.0*(p1 - p0);
	T c2 = 3.0*(p0 - 2.0*p1 + p2);
	T c3 = p3 - p0 + 3.0*(p1 - p2);

	//Find initial forward differences
	T c3t3 = c3*t3;
	T c2t2 = c2*t2;
	T c3t36 = c3t3*6.0;

	T f0 = p0;
	T f1 = c3t3 + c2t2 + c1*t;
	T f2 = c3t36 + 2.0*c2t2;
	T f3 = c3t36;

	//Create curvy goodness
	foreach(ref T r; result)
	{
		r = f0;
		f0 += f1;
		f1 += f2;
		f2 += f3;
	}
}

///Cubic hermite interpolation
double her4(double t, double p0, double m0, double p1, double m1)
{
	double t12, t22, t32;
	double t13, t23, t33;

	t12 = t*t;				//1*t^2
	t22 = t12 + t12;		//2*t^2
	t32 = t22 + t12;		//3*t^2;

	t13 = t12 * t;			//1*t^3
	t23 = t13 + t13;		//2*t^3
	t33 = t23 + t13;		//3*t^3;
	
	//Basis functions
	double b0, b1, b2, b3;
	
	b0 = t23 - t32 + 1;     //For p0
	b1 = t13 - t22 + t;     //For m0
	b2 = t32 - t23;         //For p1
	b3 = t13 - t12;         //For m1

	return b0*p0 + b1*m0 + b2*p1 + b3*m1;
}
