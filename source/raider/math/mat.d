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

	static if(D == 2)
	{
		Vec_ v[D] = [Vec_(1,0), Vec_(0,1)];
		this(T=int)(Vec_ v0, Vec_ v1)
		{
			v[0] = v0;
			v[1] = v1;
		}
	}

	static if(D == 3)
	{
		Vec_ v[D] = [Vec_(1,0,0), Vec_(0,1,0), Vec_(0,0,1)];
		this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2)
		{
			v[0] = v0;
			v[1] = v1;
			v[2] = v2;
		}
	}

	static if(D == 4)
	{
		Vec_ v[D] = [Vec_(1,0,0,0), Vec_(0,1,0,0), Vec_(0,0,1,0), Vec_(0,0,0,1)];
		this(T=int)(Vec_ v0, Vec_ v1, Vec_ v2, Vec_ v3)
		{
			v[0] = v0;
			v[1] = v1;
			v[2] = v2;
			v[3] = v3;
		}

		//Construct affine rotation and translation matrix
		this(T=int)(Mat!(3, T) o, Vec!(3, T) p)
		{
			v[0][0..3] = o[0][]; v[0][3] = 0.0;
			v[1][0..3] = o[1][]; v[1][3] = 0.0;
			v[2][0..3] = o[2][]; v[2][3] = 0.0;
			v[3][0..3] = p[]; 	 v[3][3] = 1.0;
		}
	}

	//Construct submatrix, cutting out row r and column c.
	static if(D < 4)
	{
		this(T=int)(Mat!(D+1, F) o, uint r, uint c)
		{
			assert(r < D);
			assert(c < D);

			for(uint ri=0; ri<D; ri++)
			{
				for(uint ci=0; ci<D; ci++)
				{
					v[ri][ci] = o[ri+(ri >= r)][ci+(ci >= c)];
				}
			}
		}
	}

	this(T)(const T* o) if(isMatType!T)
	{
		for(uint x=0; x<D*D; x++)
			ptr[x] = cast(F) o[x];
	}

	void identify()
	{
		for(uint x=0; x<D; x++)
			for(uint y=0; y<D; y++)
				v[x][y] = x==y ? 1.0 : 0.0;
	}

	/*
	@property static Mat_ identity()
	{
		Mat_ m;
		for(uint x=0; x<D; x++) m[x][x] = 1.0;
		return m;
	}*/

	void opAssign(T)(Mat!(D, T) o)
	{
		for(uint x=0; x<D; x++) v[x] = o[x];
	}

	const Vec_ opIndex(const size_t i) { return v[i]; }
	ref   Vec_ opIndex(const size_t i) { return v[i]; }

	@property F* ptr() { return v[0].ptr; }

	//matrix = this * matrix
	Mat_ opBinary(string op)(Mat_ o) if(op == "*")
	{
		Mat!(D, T) result;
		for(uint i=0; i<D; i++)
			for(uint j=0; j<D; j++)
				for(uint k=0; k<D; k++)
					result[i][j] += v[i][k] * o[k][j];
		return result;
	}

	//vec = this * vec
	Vec_ opBinary(string op)(Vec_ o) if(op == "*")
	{
		Vec_ result;
		for(uint r=0; r<D; r++)
			for(uint c=0; c<D; c++)
				result[r] += o[c]*v[r][c];
		return result;
	}

	//matrix = this * scalar
	Mat_ opBinary(string op)(F o) if(op == "*")
	{
		Mat_ result;
		result.ptr[0..D*D] = ptr[0..D*D]*o;
		return result;
	}

	Mat_ transposed()
	{
		Mat_ t;
		for(uint r=0; r<D; r++)
			for(uint c=0; c<D; c++)
				t[c][r] = v[r][c];
		return t;
	}

	/**
	 * Compute the adjoint matrix (that is, the transposed cofactor matrix).
	 * 
	 * Inefficient for 4x4.
	 */
	Mat_ adjoint()
	{
		//There is no 1x1 submatrix constructor. Implement 2x2 adjoint directly.
		static if(D == 2) return Mat_(Vec_(v[1][1], -v[1][0]), Vec_(v[0][1], -v[0][0]));
		else
		{
			Mat_ adj;

			for(uint r=0; r<D; r++)
				for(uint c=0; c<D; c++)
					adj[c][r] = Mat!(D-1, F)(this, r, c).determinant * ((r+c)%2 ? -1 : 1);

			return adj;
		}
	}

	F determinant()
	{
		static if(D == 2) return v[0][0]*v[1][1] - v[0][1]*v[1][0];
		else
		{
			F result = 0.0;
			for(uint c=0; c<D; c++)
				result += v[0][c] * Mat!(D-1, F)(this, 0, c).determinant * (c%2 ? -1 : 1);
			return result;
		}
	}

	static if(D == 2)
	{
		void invert()
		{
			F det = determinant;
			if(det == 0) return;

			v[0][1] /= -det;
			v[1][0] /= -det;
			F temp = v[1][1]/det;
			v[1][1] = v[0][0]/det;
			v[0][0] = temp;
		}

		Mat_ rotation(F angle)
		{
			angle *= PI/180;
			F c = cos(angle);
			F s = sin(angle);

			return Mat_(Vec_(c, s), Vec_(-s, c));
		}
	}

	static if(D == 3)
	{
		void invert()
		{
			//m stands for minor
			F m0 = v[1][1]*v[2][2] - v[1][2]*v[2][1];
			F m1 = v[2][2]*v[1][0] - v[1][2]*v[2][0];
			F m2 = v[1][0]*v[2][1] - v[1][1]*v[2][0];
			F m3 = v[0][1]*v[2][2] - v[0][2]*v[2][1];
			F m4 = v[0][0]*v[2][2] - v[0][2]*v[2][0];
			F m5 = v[0][0]*v[2][1] - v[0][1]*v[2][0];
			F m6 = v[0][1]*v[1][2] - v[0][2]*v[1][1];
			F m7 = v[0][0]*v[1][2] - v[0][2]*v[1][0];
			F m8 = v[0][0]*v[1][1] - v[0][1]*v[1][0];
				 

			F det = 
				v[0][0] * m0 -
				v[0][1] * m1 +
				v[0][2] * m2;

			if(det == 0.0) return;

			double invdet = 1.0/det;

			Mat_ i;

			//Notice this adds the sign to produce the cofactor matrix, 
			//transposes to produce the adjoint matrix,
			//and divides by the determinant to produce the inverse matrix.
			i[0][0] =  m0 * invdet;
			i[0][1] = -m3 * invdet;
			i[0][2] =  m6 * invdet;
			i[1][0] = -m1 * invdet;
			i[1][1] =  m4 * invdet;
			i[1][2] = -m7 * invdet;
			i[2][0] =  m2 * invdet;
			i[2][1] = -m5 * invdet;
			i[2][2] =  m8 * invdet;

			this = i;
		}

		Mat_ rotation(double angle, int axis)
		{
			assert(0 <= axis && axis < D);
			angle *= PI/180;

			F c = cos(angle);
			F s = sin(angle);

			switch(axis)
			{
				case 0: return Mat_(Vec_( 1, 0, 0), 
					                Vec_( 0, c, s), 
					                Vec_( 0,-s, c));
				case 1: return Mat_(Vec_( c, 0,-s), 
					                Vec_( 0, 1, 0), 
					                Vec_( s, 0, c));
				case 2: return Mat_(Vec_( c, s, 0), 
					                Vec_(-s, c, 0), 
					                Vec_( 0, 0, 1));
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
			F s0 = v[0][0] * v[1][1] - v[1][0] * v[0][1];
			F s1 = v[0][0] * v[1][2] - v[1][0] * v[0][2];
			F s2 = v[0][0] * v[1][3] - v[1][0] * v[0][3];
			F s3 = v[0][1] * v[1][2] - v[1][1] * v[0][2];
			F s4 = v[0][1] * v[1][3] - v[1][1] * v[0][3];
			F s5 = v[0][2] * v[1][3] - v[1][2] * v[0][3];

			//c stands for complementary submatrix
			F c5 = v[2][2] * v[3][3] - v[3][2] * v[2][3];
			F c4 = v[2][1] * v[3][3] - v[3][1] * v[2][3];
			F c3 = v[2][1] * v[3][2] - v[3][1] * v[2][2];
			F c2 = v[2][0] * v[3][3] - v[3][0] * v[2][3];
			F c1 = v[2][0] * v[3][2] - v[3][0] * v[2][2];
			F c0 = v[2][0] * v[3][1] - v[3][0] * v[2][1];
			
			F det =
					s0 * c5 -
					s1 * c4 +
					s2 * c3 +
					s3 * c2 -
					s4 * c1 +
					s5 * c0;
			
			if(det == 0.0) return;
			
			F invdet = 1.0/det;
			
			Mat_ i;
			
			i[0][0] = ( v[1][1] * c5 - v[1][2] * c4 + v[1][3] * c3) * invdet;
			i[0][1] = (-v[0][1] * c5 + v[0][2] * c4 - v[0][3] * c3) * invdet;
			i[0][2] = ( v[3][1] * s5 - v[3][2] * s4 + v[3][3] * s3) * invdet;
			i[0][3] = (-v[2][1] * s5 + v[2][2] * s4 - v[2][3] * s3) * invdet;
			
			i[1][0] = (-v[1][0] * c5 + v[1][2] * c2 - v[1][3] * c1) * invdet;
			i[1][1] = ( v[0][0] * c5 - v[0][2] * c2 + v[0][3] * c1) * invdet;
			i[1][2] = (-v[3][0] * s5 + v[3][2] * s2 - v[3][3] * s1) * invdet;
			i[1][3] = ( v[2][0] * s5 - v[2][2] * s2 + v[2][3] * s1) * invdet;
			
			i[2][0] = ( v[1][0] * c4 - v[1][1] * c2 + v[1][3] * c0) * invdet;
			i[2][1] = (-v[0][0] * c4 + v[0][1] * c2 - v[0][3] * c0) * invdet;
			i[2][2] = ( v[3][0] * s4 - v[3][1] * s2 + v[3][3] * s0) * invdet;
			i[2][3] = (-v[2][0] * s4 + v[2][1] * s2 - v[2][3] * s0) * invdet;
			
			i[3][0] = (-v[1][0] * c3 + v[1][1] * c1 - v[1][2] * c0) * invdet;
			i[3][1] = ( v[0][0] * c3 - v[0][1] * c1 + v[0][2] * c0) * invdet;
			i[3][2] = (-v[3][0] * s3 + v[3][1] * s1 - v[3][2] * s0) * invdet;
			i[3][3] = ( v[2][0] * s3 - v[2][1] * s1 + v[2][2] * s0) * invdet;
			
			this = i;
		}

		void invertAffine()
		{
			//Invert rotation and scale
			Mat!(3, F) r = Mat!(3, F)(this, 3, 3);
			r.invert();

			//Invert translation
			Vec!(3, F) t = Vec!(3, F)(-v[3][0], -v[3][1], -v[3][2]);
			t = t * r;

			this = Mat_(r, t);
		}
	}
}

