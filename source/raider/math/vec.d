module raider.math.vec;

import std.conv;
import std.math;

import raider.math.mat;

alias Vec!(2, double) vec2;
alias Vec!(3, double) vec3;
alias Vec!(4, double) vec4;
alias Vec!(2, float) vec2f;
alias Vec!(3, float) vec3f;
alias Vec!(4, float) vec4f;
alias Vec!(2, int) vec2i;
alias Vec!(3, int) vec3i;
alias Vec!(4, int) vec4i;
alias Vec!(2, uint) vec2u;
alias Vec!(3, uint) vec3u;
alias Vec!(4, uint) vec4u;

package template isVecType(T)
{
	const isVecType = is(T == double) || is(T == float) || is(T == int) || is(T == uint);
}

private template fer(string stuff)
{
	const char[] fer = "foreach(x; 0..D) {" ~ stuff ~ "}";
}

struct Vec(int _D, _F) if( 2 <= _D && _D <= 4 && isVecType!_F)
{
	static const int D = _D;
	alias _F F;

	alias Vec!(D, F) Vec_;

	F f[D] = 0;

	this(T)(in T s)    if(isVecType!T) { f[] = cast(F)s; }
	this(T)(in T* o)   if(isVecType!T) { mixin(fer!("f[x] = cast(F) o[x];")); }
	this(T)(in T[] o)  if(isVecType!T) { mixin(fer!("if(x < o.length) f[x] = cast(F)o[x];")); }
	//this(T)(in T[D] o) if(isVecType!T) { mixin(fer!("f[x] = cast(F) o[x];")); }
	this(T)(in Vec!(D, T) o)           { mixin(fer!("f[x] = cast(F) o[x];")); }
	void opAssign(T)(in Vec!(D, T) o)  { mixin(fer!("f[x] = cast(F) o[x];")); }

	//TODO Why do method overloads need the same template parameters?
	static if(D == 2) this(T=int)(in F x, in F y) 					{ f = [x, y]; }
	static if(D == 3) this(T=int)(in F x, in F y, in F z) 			{ f = [x, y, z]; }
	static if(D == 4) this(T=int)(in F x, in F y, in F z, in F w) 	{ f = [x, y, z, w]; }

	const F opIndex(const size_t i) { return f[i]; }
	ref   F opIndex(const size_t i) { return f[i]; }

	F[] opSlice() { return f[]; }
	F[] opSlice(size_t x, size_t y) { return f[x..y]; }

	void opSliceAssign(in F[] t) { f[] = t[]; }
	void opSliceAssign(in F[] t, size_t x, size_t y) { f[x..y] = t[]; }

	@property F* ptr()
	{
		return f.ptr;
	}

	//this op= vec
	void opOpAssign(string op, T)(Vec!(D, T) o)
	{
		mixin(fer!("f[x]" ~ op ~ "= o.f[x];"));
	}

	//this op= scalar
	void opOpAssign(string op, T)(T o) if(isVecType!T)
	{
		mixin(fer!("f[x]" ~ op ~ "= o;"));
	}

	//op this
	Vec_ opUnary(string s)() if(s == "-")
	{
		Vec_ result;
		mixin(fer!("result[x] = -f[x];"));
		return result;
	}

	static if(isMatType!F)
	{
		//vec = this * matrix
		Vec_ opBinary(string op, T)(Mat!(D, T) o) if(op == "*")
		{
			Vec_ result;
			for(uint r=0; r<D; r++)
				for(uint c=0; c<D; c++)
					result[c] += f[r]*o[r][c];
			return result;
		}
	}

	//vec = this op vec
	Vec!(D, typeof(T + F)) opBinary(string op, T)(Vec!(D, T)o)
	{
		typeof(return) r;
		mixin(fer!("r[x] = f[x]" ~ op ~ " o.f[x];"));
		return r;
	}

	//vec = this op scalar
	Vec!(D, typeof(T + F)) opBinary(string op, T)(T o) if(isVecType!T)
	{
		typeof(return) r;
		mixin(fer!("r[x] = f[x]" ~ op ~ " o;"));
		return r;
	}

	///Returns a string containing the alias of the vec instantiation followed by its elements inside parentheses.
	//eg vec3f(0.0, 1.2, 9.8)
	string toString()
	{
		string result = "vec" ~ to!string(D);
		static if(is(F == double)) result ~= "(";
		static if(is(F == float)) result ~= "f(";
		static if(is(F == uint)) result ~= "u(";
		static if(is(F == int)) result ~= "i(";

		foreach(uint x, F ff; f)
		{
			result ~= to!string(ff);
			if(x != D-1) result ~= ", ";
		}
		result ~= ")";
		return result;
	}

	///Returns a copy with the absolute value of each element.
	Vec_ abs()
	{
		Vec_ r;
		mixin(fer!("r[x] = std.math.abs(f[x]);"));
		return r;
	}

	///Angle 'twixt this and another vec with same dimensions.
	double angle(T)(Vec!(D, T) o)
	{
		Vec!(D, double) v = o;
		Vec!(D, double) t = this;
		v.normalize();
		t.normalize();
		return acos(t.dot(v));
	}

	///Dot product with another vec with same dimensions.
	double dot(T)(Vec!(D, T) o)
	{
		double result = 0.0;
		mixin(fer!("result += f[x] * o[x];"));
		return result;
	}

	///Set length to 1.
	void normalize()
	{
		length = 1.0;
	}

	///Clamp length between a minimum and maximum value.
	void clampLength(double min, double max)
	{
		double length2 = dot(this);
		if(length2 != 0.0)
		{
			if(length2 > max*2) this *= max / sqrt(length2);
			if(length2 < min*2) this *= min / sqrt(length2);
		}
	}

	alias clampLength clampMagnitude;

	///Length of vector.
	@property double length()
	{
		return sqrt(dot(this));
	}

	///ditto
	@property void length(double value)
	{
		double length2 = dot(this);
		if(length2 != 0.0) this *= (value / sqrt(length2));
	}

	alias length magnitude;

	static if(D == 2 && !is(T == uint))
	{
		///Rotate a signed vec2 clockwise.
		void roll()
		{
			F temp = f[0];
			f[0] = f[1];
			f[1] = -temp;
		}
	}

	static if(D == 3 && (is(F == double) || is(F == float)))
	{
		/**
		 * Cross product.
		 * 
		 * With your right hand, point index finger forwards,
		 * middle finger to the left, and thumb upwards. 
		 * thumb = index.cross(middle)
		 * 
		 * In other words, this is a right cross. Powie!
		 * Only available on vec3 and vec3f.
		 */
		Vec!(D, F) cross(Vec!(D, F) o)
		{
			return Vec!(D, F)(
				f[1]*o[2] - f[2]*o[1],
				f[2]*o[0] - f[0]*o[2],
				f[0]*o[1] - f[1]*o[0]);
		}
	}

}

unittest
{
	//UUUUAAAAAGH I HAVE NO IDEA WHAT I'M DOING
	//TODO Proper unittests.

	bool eq(T)(T x, T y) { return feqrel(x, y) > 20; }

	//Constructors
	vec3 a3d = vec3();					assert(eq(a3d[0], 0.0) && eq(a3d[1], 0.0));

	vec2  a2d = vec2(-2.4, 0);			assert(eq(a2d[0], -2.4) && eq(a2d[1], 0.0));
	vec3i a3i = vec3i(-2, 0, 3);		assert(a3i[0] == -2 && a3i[1] == 0 && a3i[2] == 3);
	vec4f a4f = vec4f(-2.6, 0, 4, 2.3);	assert( eq(a4f[0], -2.6f) && 
												eq(a4f[1], 0.0f) && 
	       										eq(a4f[2], 4.0f) && 
	       										eq(a4f[3], 2.3f));
	//Casting constructors
	vec3 a = vec3(1.0, -2.3, 0.0);
	vec3i b = a; 						assert(b[0] == 1 && b[1] == -2 && b[2] == 0);

	//Casts
	b = a;								assert(b[0] == 1 && b[1] == -2 && b[2] == 0);

	//[] reference assignment
	b[2] = 2;							assert(b[2] == 2);

	//Binary assignment operators
	a *= 2.0;							assert(eq(a[0], 2.0) && eq(a[1], -4.6) && eq(a[2], 0.0));
	a /= b;								assert(eq(a[0], 2.0) && eq(a[1], 2.3) && eq(a[2], 0.0));

	//Binary operators
	b = a + a - a*4;					assert(b[0] == -4 && b[1] == -4 && b[2] == 0);

	//abs
	a = vec3(-2.0, -0.0, 1.0);
	a = a.abs();						assert(eq(a[0], 2.0) && eq(a[1], 0.0) && eq(a[2], 1.0));
}


