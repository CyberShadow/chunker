module chunker.polynomials;

version(unittest) import std.format : format;
version(unittest) import std.stdio : stderr;
version = benchmark;
version (benchmark) enum N = 1;

/// Pol is a polynomial from F_2[X].
alias Pol = ulong;

/// Add returns x+y.
public Pol Add(/*this*/ Pol x, Pol y) {
	auto r = Pol(ulong(x) ^ ulong(y));
	return r;
}

unittest
{
	struct Test {
		private Pol x, y;
		private Pol sum;
	}
	Test[] tests = [
		{23, 16, 23 ^ 16},
		{0x9a7e30d1e855e0a0, 0x670102a1f4bcd414, 0xfd7f32701ce934b4},
		{0x9a7e30d1e855e0a0, 0x9a7e30d1e855e0a0, 0},
	];

	foreach (i, test; tests) {
		if (test.sum != test.x.Add(test.y)) {
			assert(false, format!"test %d failed: sum != x+y"( i));
		}

		if (test.sum != test.y.Add(test.x)) {
			assert(false, format!"test %d failed: sum != y+x"( i));
		}
	}
}


/// mulOverflows returns true if the multiplication would overflow ulong.
/// Code by Rob Pike, see
/// https://groups.google.com/d/msg/golang-nuts/h5oSN5t3Au4/KaNQREhZh0QJ
private bool mulOverflows(Pol a, Pol b) {
	if (a <= 1 || b <= 1) {
		return false;
	}
	auto c = a.mul(b);
	auto d = c.Div(b);
	if (d != a) {
		return true;
	}

	return false;
}

private Pol mul(/*this*/ Pol x, Pol y) {
	if (x == 0 || y == 0) {
		return 0;
	}

	Pol res;
	for (auto i = 0; i <= y.Deg(); i++) {
		if ((y & (1L << uint(i))) > 0) {
			res = res.Add(x << uint(i));
		}
	}

	return res;
}

/// Mul returns x*y. When an overflow occurs, Mul panics.
public Pol Mul(/*this*/ Pol x, Pol y) {
	if (mulOverflows(x, y)) {
		throw new Exception("multiplication would overflow ulong");
	}

	return x.mul(y);
}

version(unittest) private Pol parseBin(string s) {
	import std.conv : to;
	auto i = s.to!Pol(2);

	return Pol(i);
}

unittest
{
	struct Test {
		private Pol x, y;
		private Pol res;
	}
	Test[] tests = [
		{1, 2, 2},
		{
			parseBin("1101"),
			parseBin("10"),
			parseBin("11010"),
		},
		{
			parseBin("1101"),
			parseBin("11"),
			parseBin("10111"),
		},
		{
			0x40000000,
			0x40000000,
			0x1000000000000000,
		},
		{
			parseBin("1010"),
			parseBin("100100"),
			parseBin("101101000"),
		},
		{
			parseBin("100"),
			parseBin("11"),
			parseBin("1100"),
		},
		{
			parseBin("11"),
			parseBin("110101"),
			parseBin("1011111"),
		},
		{
			parseBin("10011"),
			parseBin("110101"),
			parseBin("1100001111"),
		},
	];

	foreach (i, test; tests) {
		auto m = test.x.Mul(test.y);
		if (test.res != m) {
			assert(false, format!"TestPolMul failed for test %d: %s * %s: want %s, got %s"(
				i, test.x, test.y, test.res, m));
		}
		m = test.y.Mul(test.x);
		if (test.res != test.y.Mul(test.x)) {
			assert(false, format!"TestPolMul failed for %d: %s * %s: want %s, got %s"(
				i, test.x, test.y, test.res, m));
		}
	}
}

unittest
{
	auto x = Pol(1L << 63);
	try
	{
		x.Mul(2);
		assert(false, "overflow test did not panic");
	}
	catch (Exception e)
	{
		// try to recover overflow error
		if (e.msg == "multiplication would overflow ulong")
			return;

		stderr.writefln!"invalid error raised: %s"(e);
		// re-raise error if not overflow
		throw e;
	}
}


/// Deg returns the degree of the polynomial x. If x is zero, -1 is returned.
public int Deg(/*this*/ Pol x) {
	// the degree of 0 is -1
	if (x == 0) {
		return -1;
	}

	// see https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog

	auto r = 0;
	if ((ulong(x) & 0xffffffff00000000) > 0) {
		x >>= 32;
		r |= 32;
	}

	if ((ulong(x) & 0xffff0000) > 0) {
		x >>= 16;
		r |= 16;
	}

	if ((ulong(x) & 0xff00) > 0) {
		x >>= 8;
		r |= 8;
	}

	if ((ulong(x) & 0xf0) > 0) {
		x >>= 4;
		r |= 4;
	}

	if ((ulong(x) & 0xc) > 0) {
		x >>= 2;
		r |= 2;
	}

	if ((ulong(x) & 0x2) > 0) {
		r |= 1;
	}

	return r;
}

unittest
{
	Pol x;
	if (x.Deg() != -1) {
		assert(false, format!"deg(0) is not -1: %s"( x.Deg()));
	}

	x = 1;
	if (x.Deg() != 0) {
		assert(false, format!"deg(1) is not 0: %s"( x.Deg()));
	}

	for (auto i = 0; i < 64; i++) {
		x = 1L << uint(i);
		if (x.Deg() != i) {
			assert(false, format!"deg(1<<%d) is not %d: %s"( i, i, x.Deg()));
		}
	}
}

unittest // benchmark
{
	auto f = Pol(0x3af4b284899);
	auto d = f.Deg();
	if (d != 41) {
		assert(false, format!"BenchmalPolDeg: Wrong degree %d returned, expected %d"
			(d, 41));
	}

	for (auto i = 0; i < N; i++) {
		f.Deg();
	}
}


/// String returns the coefficients in hex.
public string String(/*this*/ Pol x) {
	import std.conv : to;
	return "0x" ~ to!string(ulong(x), 16);
}

/// Expand returns the string representation of the polynomial x.
public string Expand(/*this*/ Pol x) {
	import std.format : format;

	if (x == 0) {
		return "0";
	}

	auto s = "";
	for (auto i = x.Deg(); i > 1; i--) {
		if ((x&(1L<<uint(i))) > 0) {
			s ~= format!"+x^%d"(i);
		}
	}

	if ((x&2) > 0) {
		s ~= "+x";
	}

	if ((x&1) > 0) {
		s ~= "+1";
	}

	return s[1 .. $];
}

unittest
{
	auto pol = Pol(0x3DA3358B4DC173);
	auto s = pol.Expand();
	if (s != "x^53+x^52+x^51+x^50+x^48+x^47+x^45+x^41+x^40+x^37+x^36+x^34+x^32+x^31+x^27+x^25+x^24+x^22+x^19+x^18+x^16+x^15+x^14+x^8+x^6+x^5+x^4+x+1") {
		assert(false, "wrong result");
	}
}


/// DivMod returns x / d = q, and remainder r,
/// see https://en.wikipedia.org/wiki/Division_algorithm
public Pol[2] DivMod(/*this*/ Pol x, Pol d) {
	if (x == 0) {
		return [0, 0];
	}

	if (d == 0) {
		assert(false, "division by zero");
	}

	auto D = d.Deg();
	auto diff = x.Deg() - D;
	if (diff < 0) {
		return [0, x];
	}

	Pol q;
	while (diff >= 0) {
		auto m = d << uint(diff);
		q |= (1L << uint(diff));
		x = x.Add(m);

		diff = x.Deg() - D;
	}

	return [q, x];
}

unittest // benchmark
{
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.DivMod(f);
	}
}


/// Div returns the integer division result x / d.
public Pol Div(/*this*/ Pol x, Pol d) {
	auto q = x.DivMod(d)[0];
	return q;
}

unittest
{
	struct Test {
		private Pol x, y;
		private Pol res;
	}

	Test[] tests = [
		{10, 50, 0},
		{0, 1, 0},
		{
			parseBin("101101000"), // 0x168
			parseBin("1010"),      // 0xa
			parseBin("100100"),    // 0x24
		},
		{2, 2, 1},
		{
			0x8000000000000000,
			0x8000000000000000,
			1,
		},
		{
			parseBin("1100"),
			parseBin("100"),
			parseBin("11"),
		},
		{
			parseBin("1100001111"),
			parseBin("10011"),
			parseBin("110101"),
		},
	];

	foreach (i, test; tests) {
		auto m = test.x.Div(test.y);
		if (test.res != m) {
			assert(false, format!"TestPolDiv failed for test %d: %s * %s: want %s, got %s"(
				i, test.x, test.y, test.res, m));
		}
	}
}

unittest // benchmark
{
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.Div(f);
	}
}


/// Mod returns the remainder of x / d
public Pol Mod(/*this*/ Pol x, Pol d) {
	auto r = x.DivMod(d)[1];
	return r;
}

unittest
{
	struct Test {
		private Pol x, y;
		private Pol res;
	}
	Test[] tests = [
		{10, 50, 10},
		{0, 1, 0},
		{
			parseBin("101101001"),
			parseBin("1010"),
			parseBin("1"),
		},
		{2, 2, 0},
		{
			0x8000000000000000,
			0x8000000000000000,
			0,
		},
		{
			parseBin("1100"),
			parseBin("100"),
			parseBin("0"),
		},
		{
			parseBin("1100001111"),
			parseBin("10011"),
			parseBin("0"),
		},
	];

	foreach (i, test; tests) {
		auto res = test.x.Mod(test.y);
		if (test.res != res) {
			assert(false, format!"test %d failed: want %s, got %s"( i, test.res, res));
		}
	}
}

unittest // benchmark
{
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.Mod(f);
	}
}


/// I really dislike having a function that does not terminate, so specify a
/// really large upper bound for finding a new irreducible polynomial, and
/// return an error when no irreducible polynomial has been found within
/// randPolMaxTries.
private enum randPolMaxTries = 1e6;

/// RandomPolynomial returns a new random irreducible polynomial
/// of degree 53 using the default System CSPRNG as source.
/// It is equivalent to calling DerivePolynomial(rand.Reader).
public Pol RandomPolynomial() {
	import std.random : rndGen;
	return DerivePolynomial(rndGen);
}

unittest
{
	RandomPolynomial();
}

unittest // benchmark
{
	for (auto i = 0; i < N; i++) {
		RandomPolynomial();
	}
}


/// DerivePolynomial returns an irreducible polynomial of degree 53
/// (largest prime number below 64-8) by reading bytes from source.
/// There are (2^53-2/53) irreducible polynomials of degree 53 in
/// F_2[X], c.f. Michael O. Rabin (1981): "Fingerprinting by Random
/// Polynomials", page 4. If no polynomial could be found in one
/// million tries, an error is returned.
public Pol DerivePolynomial(Random)(Random source) {
	for (auto i = 0; i < randPolMaxTries; i++) {
		Pol f;

		// choose polynomial at (pseudo)random
		import std.random : uniform;
		f = uniform!Pol(source);

		// mask away bits above bit 53
		f &= Pol((1L << 54) - 1);

		// set highest and lowest bit so that the degree is 53 and the
		// polynomial is not trivially reducible
		f |= (1L << 53) | 1;

		// test if f is irreducible
		if (f.Irreducible()) {
			return f;
		}
	}

	// If this is reached, we haven't found an irreducible polynomial in
	// randPolMaxTries. This error is very unlikely to occur.
	throw new Exception("unable to find new random irreducible polynomial");
}

/// GCD computes the Greatest Common Divisor x and f.
public Pol GCD(/*this*/ Pol x, Pol f) {
	if (f == 0) {
		return x;
	}

	if (x == 0) {
		return f;
	}

	if (x.Deg() < f.Deg()) {
		import std.algorithm.mutation : swap;
		swap(x, f);
	}

	return f.GCD(x.Mod(f));
}

unittest
{
	struct Test {
		private Pol f1;
		private Pol f2;
		private Pol gcd;
	}
	Test[] tests = [
		{10, 50, 2},
		{0, 1, 1},
		{
			parseBin("101101001"),
			parseBin("1010"),
			parseBin("1"),
		},
		{2, 2, 2},
		{
			parseBin("1010"),
			parseBin("11"),
			parseBin("11"),
		},
		{
			0x8000000000000000,
			0x8000000000000000,
			0x8000000000000000,
		},
		{
			parseBin("1100"),
			parseBin("101"),
			parseBin("11"),
		},
		{
			parseBin("1100001111"),
			parseBin("10011"),
			parseBin("10011"),
		},
		{
			0x3DA3358B4DC173,
			0x3DA3358B4DC173,
			0x3DA3358B4DC173,
		},
		{
			0x3DA3358B4DC173,
			0x230d2259defd,
			1,
		},
		{
			0x230d2259defd,
			0x51b492b3eff2,
			parseBin("10011"),
		},
	];

	foreach (i, test; tests) {
		auto gcd = test.f1.GCD(test.f2);
		if (test.gcd != gcd) {
			assert(false, format!"GCD test %d (%+s) failed: got %s, wanted %s"(
				i, test, gcd, test.gcd));
		}

		gcd = test.f2.GCD(test.f1);
		if (test.gcd != gcd) {
			assert(false, format!"GCD test %d (%+s) failed: got %s, wanted %s"(
				i, test, gcd, test.gcd));
		}
	}
}


/// Irreducible returns true iff x is irreducible over F_2. This function
/// uses Ben Or's reducibility test.
//
/// For details see "Tests and Constructions of Irreducible Polynomials over
/// Finite Fields".
public bool Irreducible(/*this*/ Pol x) {
	for (auto i = 1; i <= x.Deg()/2; i++) {
		if (x.GCD(qp(uint(i), x)) != 1) {
			return false;
		}
	}

	return true;
}

unittest
{
	struct Test {
		private Pol f;
		private bool irred;
	}
	Test[] tests = [
		{0x38f1e565e288df, false},
		{0x3DA3358B4DC173, true},
		{0x30a8295b9d5c91, false},
		{0x255f4350b962cb, false},
		{0x267f776110a235, false},
		{0x2f4dae10d41227, false},
		{0x2482734cacca49, true},
		{0x312daf4b284899, false},
		{0x29dfb6553d01d1, false},
		{0x3548245eb26257, false},
		{0x3199e7ef4211b3, false},
		{0x362f39017dae8b, false},
		{0x200d57aa6fdacb, false},
		{0x35e0a4efa1d275, false},
		{0x2ced55b026577f, false},
		{0x260b012010893d, false},
		{0x2df29cbcd59e9d, false},
		{0x3f2ac7488bd429, false},
		{0x3e5cb1711669fb, false},
		{0x226d8de57a9959, false},
		{0x3c8de80aaf5835, false},
		{0x2026a59efb219b, false},
		{0x39dfa4d13fb231, false},
		{0x3143d0464b3299, false},
	];

	foreach (_, test; tests) {
		if (test.f.Irreducible() != test.irred) {
			assert(false, format!"Irreducibility test for Polynomial %s failed: got %s, wanted %s"(
				test.f, test.f.Irreducible(), test.irred));
		}
	}

	// find first irreducible polynomial
	Pol pol;
	foreach (_, test; tests) {
		if (test.irred) {
			pol = test.f;
			break;
		}
	}

	for (auto i = 0; i < N; i++) {
		if (!pol.Irreducible()) {
			assert(false, format!"Irreducibility test for Polynomial %s failed"( pol));
		}
	}
}


/// MulMod computes x*f mod g
public Pol MulMod(/*this*/ Pol x, Pol f, Pol g) {
	if (x == 0 || f == 0) {
		return 0;
	}

	Pol res;
	for (auto i = 0; i <= f.Deg(); i++) {
		if ((f & (1L << uint(i))) > 0) {
			auto a = x;
			for (auto j = 0; j < i; j++) {
				a = a.Mul(2).Mod(g);
			}
			res = res.Add(a).Mod(g);
		}
	}

	return res;
}

unittest
{
	struct Test {
		private Pol f1;
		private Pol f2;
		private Pol g;
		private Pol mod;
	}
	Test[] tests = [
		{
			0x1230,
			0x230,
			0x55,
			0x22,
		},
		{
			0x0eae8c07dbbb3026,
			0xd5d6db9de04771de,
			0xdd2bda3b77c9,
			0x425ae8595b7a,
		},
	];

	foreach (i, test; tests) {
		auto mod = test.f1.MulMod(test.f2, test.g);
		if (mod != test.mod) {
			assert(false, format!"MulMod test %d (%+s) failed: got %s, wanted %s"(
				i, test, mod, test.mod));
		}
	}
}


/// qp computes the polynomial (x^(2^p)-x) mod g. This is needed for the
/// reducibility test.
private Pol qp(uint p, Pol g) {
	auto num = (1L << p);
	auto i = 1;

	// start with x
	auto res = Pol(2);

	while (i < num) {
		// repeatedly square res
		res = res.MulMod(res, g);
		i *= 2;
	}

	// add x
	return res.Add(2).Mod(g);
}

