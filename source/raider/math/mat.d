module raider.math.mat;

import std.math;
import raider.math.vec;

alias Mat!(2, double) mat2;
alias Mat!(3, double) mat3;
alias Mat!(4, double) mat4;
alias Mat!(2, float) mat2f;
alias Mat!(3, float) mat3f;
alias Mat!(4, float) mat4f;

package template isMatType(T)
{
	const isMatType = is(T == double) || is(T == float);
}

struct Mat(int _D, _F) if(2 <= _D && _D <= 4 && isMatType!_F)
{
	static const int D = _D;
	alias _F F;

	alias Vec!(D, F) Vec_;
	alias Mat!(D, F) Mat_;

	union
	{
		F[D*D] f = 0;
		Vec_[D] v;
	}

	static if(D == 2)
	{
		this(T=int)(Vec_ v0, Vec_ v1)
		{ v[] = [v0, v1]; }

		this(T=int)(F f0, F f1, F f2, F f3)
		{ f[] = [f0, f1, f2, f3]; }

		static Mat_ identity = Mat_(1,0,
		                            0,1);
	}

	static if(D == 3)
	{
		this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2)
		{ v[] = [v0, v1, v2]; }

		this(T=int)(F f0, F f1, F f2, F f3, F f4, F f5, F f6, F f7, F f8)
		{ f[] = [f0, f1, f2, f3, f4, f5, f6, f7, f8]; }
		
		static Mat_ identity = Mat_(1,0,0,
		                            0,1,0,
		                            0,0,1);
	}

	static if(D == 4)
	{
		this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2, Vec_ v3, Vec_ v4)
		{ v[] = [v0, v1, v2, v3, v4]; }

		this(T=int)(F f0,  F f1,  F f2,  F f3,  
		            F f4,  F f5,  F f6,  F f7,  
		            F f8,  F f9,  F f10, F f11,
		            F f12, F f13, F f14, F f15)
		{ f[] = [f0,  f1,  f2,  f3, 
				 f4,  f5,  f6,  f7,
				 f8,  f9,  f10, f11,
				 f12, f13, f14, f15]; }
		
		static Mat_ identity = Mat_(1,0,0,0,
		                            0,1,0,0,
		                            0,0,1,0,
		                            0,0,0,1);

		//Construct affine orientation and translation matrix
		this(T=int)(auto ref in Mat!(3, T) o, Vec!(3, T) p)
		{
			f[0..3]   = o.f[0..3]; f[3]  = 0.0;
			f[4..7]   = o.f[3..6]; f[7]  = 0.0;
			f[8..11]  = o.f[6..9]; f[11] = 0.0;
			f[12..15] =   p[0..3]; f[15] = 1.0;
		}
	}

	//Construct submatrix, cutting out row r and column c.
	static if(D < 4)
	{
		this(T=int)(auto ref in Mat!(D+1, F) o, uint r, uint c)
		{
			assert(r < D);
			assert(c < D);

			foreach(ri; 0..D)
				foreach(ci; 0..D)
					v[ri][ci] = o[ri+(ri >= r)][ci+(ci >= c)];
		}
	}

	this(T)(const T* o) if(isMatType!T)
	{
		foreach(x; 0..D*D) f[x] = cast(F) o[x];
	}

	void opAssign(T)(auto ref in Mat!(D, T) o)
	{
		foreach(x; 0..D*D) f[x] = cast(F) o.f[x];
	}

	const Vec_ opIndex(const size_t i) { return v[i]; }
	ref   Vec_ opIndex(const size_t i) { return v[i]; }

	@property F* ptr() { return f.ptr; }

	//matrix = this * matrix
	Mat_ opBinary(string op)(auto ref in Mat_ o) if(op == "*")
	{
		Mat!(D, T) result;
		foreach(i; 0..D)
			foreach(j; 0..D)
				foreach(k; 0..D)
					result[i][j] += v[i][k] * o[k][j];

		return result;
	}

	//vec = this * vec
	Vec_ opBinary(string op)(auto ref in Vec_ o) if(op == "*")
	{
		Vec_ result;
		foreach(r; 0..D)
			foreach(c; 0..D)
				result[r] += o[c]*v[r][c];

		return result;
	}

	//matrix = this * scalar
	Mat_ opBinary(string op)(F o) if(op == "*")
	{
		Mat_ r = void;
		r.f[0..D*D] = f[0..D*D]*o;

		return r;
	}

	Mat_ transposed()
	{
		Mat_ t = void;
		foreach(r; 0..D)
			foreach(c; 0..D)
				t[c][r] = v[r][c];

		return t;
	}

	/**
	 * Compute the adjoint matrix (that is, the transposed cofactor matrix).
	 * 
	 * Inefficient for 4x4 inversion. Use invert().
	 */
	Mat_ adjoint()
	{
		//There is no 1x1 submatrix constructor. Implement 2x2 adjoint directly.
		static if(D == 2) return Mat_(f[3], -f[2], f[1], -f[0]);
		else
		{
			Mat_ adj = void;

			foreach(r; 0..D)
				foreach(c; 0..D)
					adj[c][r] = Mat!(D-1, F)(this, r, c).determinant * ((r+c)%2 ? -1 : 1);

			return adj;
		}
	}

	F determinant()
	{
		static if(D == 2) return f[0]*f[3] - f[1]*f[2];
		else
		{
			F r = 0.0;
			foreach(c; 0..D)
				r += v[0][c] * Mat!(D-1, F)(this, 0, c).determinant * (c%2 ? -1 : 1);

			return r;
		}
	}

	static if(D == 2)
	{
		void invert()
		{
			F det = determinant;
			if(det == 0) return;

			f[1] /= -det;
			f[2] /= -det;
			F temp = f[3]/det;
			f[3] = f[0]/det;
			f[0] = temp;
		}

		Mat_ rotation(F angle)
		{
			angle *= PI/180;
			F c = cos(angle);
			F s = sin(angle);

			return Mat_(c, s, -s, c);
		}
	}

	static if(D == 3)
	{
		void invert()
		{
			//m stands for minor
			F m0 = f[4]*f[8] - f[5]*f[7];
			F m1 = f[8]*f[3] - f[5]*f[6];
			F m2 = f[3]*f[7] - f[4]*f[6];
			F m3 = f[1]*f[8] - f[2]*f[7];
			F m4 = f[0]*f[8] - f[2]*f[6];
			F m5 = f[0]*f[7] - f[1]*f[6];
			F m6 = f[1]*f[5] - f[2]*f[4];
			F m7 = f[0]*f[5] - f[2]*f[3];
			F m8 = f[0]*f[4] - f[1]*f[3];
				 

			F det = 
				f[0] * m0 -
				f[1] * m1 +
				f[2] * m2;

			if(det == 0.0) return;

			double invdet = 1.0/det;

			//Notice this adds the sign to produce the cofactor matrix, 
			//transposes to produce the adjoint matrix,
			//and divides by the determinant to produce the inverse matrix.
			f[0] =  m0 * invdet;
			f[1] = -m3 * invdet;
			f[2] =  m6 * invdet;
			f[3] = -m1 * invdet;
			f[4] =  m4 * invdet;
			f[5] = -m7 * invdet;
			f[6] =  m2 * invdet;
			f[7] = -m5 * invdet;
			f[8] =  m8 * invdet;
		}

		Mat_ rotation(double angle, int axis)
		{
			assert(0 <= axis && axis < D);
			angle *= PI/180;

			F c = cos(angle);
			F s = sin(angle);

			switch(axis)
			{
				case 0: return Mat_(1, 0, 0, 
					                0, c, s, 
					                0,-s, c);
				case 1: return Mat_(c, 0,-s, 
					                0, 1, 0, 
					                s, 0, c);
				case 2: return Mat_(c, s, 0, 
					               -s, c, 0, 
					                0, 0, 1);
				default: assert(0);
			}
		}
	}

	static if(D == 4)
	{
		void invert()
		{
			/*
			 * This is an implementation of the Laplace Expansion Theorem which is described here:
			 * http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
			 * 
			 * Instead of finding the 4x4 determinant through the determinants of four 3x3 matrices
			 * and four 1x1 matrices (expansion by columns or rows), we find the determinants of 
			 * twelve 2x2 matrices (expansion by 2 columns or 2 rows). This requires fewer operations.
			 * 
			 * And, as it happens, all those 2x2 determinants may be reused while computing the 
			 * cofactor matrix, by using particular submatrices.
			 * 
			 * Finding determinant by...
			 * row expansion = 40 multiplications, 23 additions
			 * 2 row expansion = 30 mult, 17 add
			 * 
			 * Finding cofactors by...
			 * naive expansion = 160 mult, 80 add
			 * reusing determinants = 48 mult, 32 add
			 * 
			 * Finding inverse by...
			 * row expansion & naive expansion = 216 mult, 103 add
			 * 2 row expansion & reusing determinants = 94 mult, 49 add
			 */

			//s stands for submatrix
			F s0 = f[0] * f[5] - f[4] * f[1];
			F s1 = f[0] * f[6] - f[4] * f[2];
			F s2 = f[0] * f[7] - f[4] * f[3];
			F s3 = f[1] * f[6] - f[5] * f[2];
			F s4 = f[1] * f[7] - f[5] * f[3];
			F s5 = f[2] * f[7] - f[6] * f[3];

			//c stands for complementary submatrix
			F c5 = f[10] * f[15] - f[14] * f[11];
			F c4 = f[9]  * f[15] - f[13] * f[11];
			F c3 = f[9]  * f[14] - f[13] * f[10];
			F c2 = f[8]  * f[15] - f[12] * f[11];
			F c1 = f[8]  * f[14] - f[12] * f[10];
			F c0 = f[8]  * f[13] - f[12] * f[9];
			
			F det =
					s0 * c5 -
					s1 * c4 +
					s2 * c3 +
					s3 * c2 -
					s4 * c1 +
					s5 * c0;
			
			if(det == 0.0) return;
			
			F invdet = 1.0/det;
			
			Mat_ i = void;

			i.f[0]  = ( f[5]  * c5 - f[6]  * c4 + f[7]  * c3) * invdet;
			i.f[1]  = (-f[1]  * c5 + f[2]  * c4 - f[3]  * c3) * invdet;
			i.f[2]  = ( f[13] * s5 - f[14] * s4 + f[15] * s3) * invdet;
			i.f[3]  = (-f[9]  * s5 + f[10] * s4 - f[11] * s3) * invdet;
			
			i.f[4]  = (-f[4]  * c5 + f[6]  * c2 - f[7]  * c1) * invdet;
			i.f[5]  = ( f[0]  * c5 - f[2]  * c2 + f[3]  * c1) * invdet;
			i.f[6]  = (-f[12] * s5 + f[14] * s2 - f[15] * s1) * invdet;
			i.f[7]  = ( f[8]  * s5 - f[10] * s2 + f[11] * s1) * invdet;
			
			i.f[8]  = ( f[4]  * c4 - f[5]  * c2 + f[7]  * c0) * invdet;
			i.f[9]  = (-f[0]  * c4 + f[1]  * c2 - f[3]  * c0) * invdet;
			i.f[10] = ( f[12] * s4 - f[13] * s2 + f[15] * s0) * invdet;
			i.f[11] = (-f[8]  * s4 + f[9]  * s2 - f[11] * s0) * invdet;

			i.f[12] = (-f[4]  * c3 + f[5]  * c1 - f[6]  * c0) * invdet;
			i.f[13] = ( f[0]  * c3 - f[1]  * c1 + f[2]  * c0) * invdet;
			i.f[14] = (-f[12] * s3 + f[13] * s1 - f[14] * s0) * invdet;
			i.f[15] = ( f[8]  * s3 - f[9]  * s1 + f[10] * s0) * invdet;
			
			this = i;
		}

		void invertAffine()
		{
			//Invert rotation and scale
			Mat!(3, F) r = Mat!(3, F)(this, 3, 3);
			r.invert;

			//Invert translation
			Vec!(3, F) t = Vec!(3, F)(-f[12], -f[13], -f[14]);
			t = t * r;

			this = Mat_(r, t);
		}
	}
}

