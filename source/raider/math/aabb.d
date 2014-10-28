module raider.math.aabb;

import raider.math.vec;

alias Aabb!(2, double) aabb2;
alias Aabb!(3, double) aabb3;
alias Aabb!(2, float) aabb2f;
alias Aabb!(3, float) aabb3f;
alias Aabb!(2, int) aabb2i;
alias Aabb!(3, int) aabb3i;

package template isAabbType(T)
{
	const isAabbType = is(T == double) || is(T == float) || is(T == int);
}

private template fer(string stuff)
{
	const char[] fer = "foreach(x; 0..D) {" ~ stuff ~"}";
}

struct Aabb(int _D, _F) if(2 <= _D && _D <= 3 && isAabbType!_F)
{
	static const int D = _D;
	alias _F F;

	alias Aabb!(D, F) Aabb_;
	alias Vec!(D, F) Vec_;

	Vec_ f0, f1;

	this(T)(Vec!(D, T) v0, Vec!(D, T) v1) if(isAabbType!T) { f0 = v0; f1 = v1; }
}