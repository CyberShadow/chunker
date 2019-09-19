module chunker.polynomials;

/// Pol is a polynomial from `F_2[X]`.
struct Pol
{
	alias Base = ulong;
	Base value;

	/// Returns `this+y`.
	public Pol opBinary(string op)(Pol y) const
	if (op == "+")
	{
		return Pol(this.value ^ y.value);
	}

	unittest
	{
		struct Test
		{
			Pol x, y;
			Pol sum;
		}
		Test[] tests =
		[
			{Pol(23), Pol(16), Pol(23 ^ 16)},
			{Pol(0x9a7e30d1e855e0a0), Pol(0x670102a1f4bcd414), Pol(0xfd7f32701ce934b4)},
			{Pol(0x9a7e30d1e855e0a0), Pol(0x9a7e30d1e855e0a0), Pol(0)},
		];

		foreach (i, test; tests)
		{
			assert(test.sum == test.x + test.y,
				format!"test %d failed: sum != x+y"(i));

			assert(test.sum == test.y + test.x,
				format!"test %d failed: sum != y+x"(i));
		}
	}


	/// Returns true if the multiplication would overflow `Pol.Base`.
	/// Code by Rob Pike, see
	/// https://groups.google.com/d/msg/golang-nuts/h5oSN5t3Au4/KaNQREhZh0QJ
	private static bool mulOverflows(Pol a, Pol b)
	{
		if (a.value <= 1 || b.value <= 1)
			return false;
		auto c = mulImpl(a, b);
		auto d = c / b;
		if (d != a)
			return true;

		return false;
	}

	private static Pol mulImpl(Pol x, Pol y)
	{
		if (x.value == 0 || y.value == 0)
			return Pol(0);

		Pol res;
		foreach (i; 0 .. y.deg + 1)
			if ((y.value & (1L << uint(i))) > 0)
				res = res + Pol(x.value << uint(i));

		return res;
	}

	/// Returns `this*y`.
	/// When an overflow occurs, throws an exception.
	public Pol opBinary(string op)(Pol y) const
	if (op == "*")
	{
		if (mulOverflows(this, y))
			throw new Exception("multiplication would overflow ulong");

		return mulImpl(this, y);
	}

	version(unittest) static private Pol parseBin(string s)
	{
		import std.conv : to;
		return Pol(s.to!Base(2));
	}

	unittest
	{
		struct Test
		{
			Pol x, y;
			Pol res;
		}
		Test[] tests =
		[
			{Pol(1), Pol(2), Pol(2)},
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
				Pol(0x40000000),
				Pol(0x40000000),
				Pol(0x1000000000000000),
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

		foreach (i, test; tests)
		{
			auto m = test.x * test.y;
			assert(test.res == m,
				format!"TestPolMul failed for test %d: %s * %s: want %s, got %s"
				(i, test.x, test.y, test.res, m));

			m = test.y * test.x;
			assert(test.res == m,
				format!"TestPolMul failed for %d: %s * %s: want %s, got %s"
				(i, test.x, test.y, test.res, m));
		}
	}

	unittest
	{
		auto x = Pol(1L << 63);
		try
		{
			x * Pol(2);
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


	/// Returns the degree of the polynomial.
	/// If `this` is zero, -1 is returned.
	@property public int deg() const
	{
		Base x = this.value;

		// the degree of 0 is -1
		if (x == 0)
			return -1;

		// see https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog

		auto r = 0;
		if ((ulong(x) & 0xffffffff00000000) > 0)
		{
			x >>= 32;
			r |= 32;
		}

		if ((ulong(x) & 0xffff0000) > 0)
		{
			x >>= 16;
			r |= 16;
		}

		if ((ulong(x) & 0xff00) > 0)
		{
			x >>= 8;
			r |= 8;
		}

		if ((ulong(x) & 0xf0) > 0)
		{
			x >>= 4;
			r |= 4;
		}

		if ((ulong(x) & 0xc) > 0)
		{
			x >>= 2;
			r |= 2;
		}

		if ((ulong(x) & 0x2) > 0)
		{
			r |= 1;
		}

		return r;
	}

	unittest
	{
		Pol x;
		assert(x.deg == -1,
			format!"deg(0) is not -1: %s"(x.deg));

		x = Pol(1);
		assert(x.deg == 0,
			format!"deg(1) is not 0: %s"(x.deg));

		foreach (i; 0 .. 64)
		{
			x = Pol(1L << uint(i));
			assert(x.deg == i,
				format!"deg(1<<%d) is not %d: %s"(i, i, x.deg));
		}
	}

	version(benchmarkPolynomials) private static void benchmarkPolDeg()
	{
		auto f = Pol(0x3af4b284899);
		auto d = f.deg;
		assert(d == 41,
			format!"BenchmarkPolDeg: Wrong degree %d returned, expected %d"
			(d, 41));

		foreach (i; 0 .. Benchmark.N)
			f.deg;
	}


	/// Returns the coefficients in hex.
	public string toString() const
	{
		import std.conv : to;
		return "0x" ~ to!string(this.value, 16);
	}

	/// Returns the string representation of the polynomial x.
	public string expand() const
	{
		import std.format : format;

		if (value == 0)
			return "0";

		auto s = "";
		for (auto i = deg; i > 1; i--)
			if ((value & (1L<<uint(i))) > 0)
				s ~= format!"+x^%d"(i);

		if ((value & 2) > 0)
			s ~= "+x";

		if ((value & 1) > 0)
			s ~= "+1";

		return s[1 .. $];
	}

	unittest
	{
		auto pol = Pol(0x3DA3358B4DC173);
		auto s = pol.expand();
		assert(s == "x^53+x^52+x^51+x^50+x^48+x^47+x^45+x^41+x^40+x^37+x^36+x^34+x^32+x^31+x^27+x^25+x^24+x^22+x^19+x^18+x^16+x^15+x^14+x^8+x^6+x^5+x^4+x+1",
			"wrong result");
	}


	/// Returns quotient and remainder from division `[x / d, x % d]`,
	/// see https://en.wikipedia.org/wiki/Division_algorithm
	public static Pol[2] divMod(Pol x, Pol d)
	{
		if (x.value == 0)
			return [Pol(0), Pol(0)];

		if (d.value == 0)
			assert(false, "division by zero");

		auto d_deg = d.deg;
		auto diff = x.deg - d_deg;
		if (diff < 0)
			return [Pol(0), x];

		Pol q;
		while (diff >= 0)
		{
			auto m = d.value << uint(diff);
			q.value |= (1L << uint(diff));
			x = x + Pol(m);

			diff = x.deg - d_deg;
		}

		return [q, x];
	}

	version(benchmarkPolynomials) private static void benchmarkPolDivMod()
	{
		auto f = Pol(0x2482734cacca49);
		auto g = Pol(0x3af4b284899);

		foreach (i; 0 .. Benchmark.N)
			divMod(g, f);
	}


	/// Returns the integer division result `this / d`.
	public Pol opBinary(string op)(Pol d) const
	if (op == "/")
	{
		return divMod(this, d)[0];
	}

	unittest
	{
		struct Test
		{
			private Pol x, y;
			private Pol res;
		}

		Test[] tests =
		[
			{Pol(10), Pol(50), Pol(0)},
			{Pol(0), Pol(1), Pol(0)},
			{
				parseBin("101101000"), // 0x168
				parseBin("1010"),      // 0xa
				parseBin("100100"),    // 0x24
			},
			{Pol(2), Pol(2), Pol(1)},
			{
				Pol(0x8000000000000000),
				Pol(0x8000000000000000),
				Pol(1),
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

		foreach (i, test; tests)
		{
			auto m = test.x / test.y;
			assert(test.res == m,
				format!"TestPolDiv failed for test %d: %s * %s: want %s, got %s"
				(i, test.x, test.y, test.res, m));
		}
	}

	version(benchmarkPolynomials) private static void benchmarkPolDiv()
	{
		auto f = Pol(0x2482734cacca49);
		auto g = Pol(0x3af4b284899);

		foreach (i; 0 .. Benchmark.N)
			g / f;
	}


	/// Returns the remainder of `this / d`.
	public Pol opBinary(string op)(Pol d) const
	if (op == "%")
	{
		return divMod(this, d)[1];
	}

	unittest
	{
		struct Test
		{
			private Pol x, y;
			private Pol res;
		}
		Test[] tests =
		[
			{Pol(10), Pol(50), Pol(10)},
			{Pol(0), Pol(1), Pol(0)},
			{
				parseBin("101101001"),
				parseBin("1010"),
				parseBin("1"),
			},
			{Pol(2), Pol(2), Pol(0)},
			{
				Pol(0x8000000000000000),
				Pol(0x8000000000000000),
				Pol(0),
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

		foreach (i, test; tests)
		{
			auto res = test.x % test.y;
			assert(test.res == res,
				format!"test %d failed: want %s, got %s"(i, test.res, res));
		}
	}

	version(benchmarkPolynomials) private static void benchmarkPolMod()
	{
		auto f = Pol(0x2482734cacca49);
		auto g = Pol(0x3af4b284899);

		foreach (i; 0 .. Benchmark.N)
			g % f;
	}


	/// I really dislike having a function that does not terminate, so specify a
	/// really large upper bound for finding a new irreducible polynomial, and
	/// return an error when no irreducible polynomial has been found within
	/// randPolMaxTries.
	private enum randPolMaxTries = 1e6;

	/// Returns a new random irreducible polynomial of degree 53 using
	/// the default `rndGen` as source.  It is equivalent to calling
	/// `derive(rndGen)`.
	public static Pol getRandom()
	{
		import std.random : rndGen;
		return derive(rndGen);
	}

	unittest
	{
		getRandom();
	}

	version(benchmarkPolynomials) private static void benchmarkPolGetRandom()
	{
		foreach (i; 0 .. Benchmark.N)
			getRandom();
	}


	/// Returns an irreducible polynomial of degree 53
	/// (largest prime number below 64-8) by reading bytes from source.
	/// There are (2^53-2/53) irreducible polynomials of degree 53 in
	/// `F_2[X]`, c.f. Michael O. Rabin (1981): "Fingerprinting by Random
	/// Polynomials", page 4. If no polynomial could be found in one
	/// million tries, an error is returned.
	public static Pol derive(Random)(Random source)
	{
		foreach (i; 0 .. randPolMaxTries)
		{
			Pol f;

			// choose polynomial at (pseudo)random
			import std.random : uniform;
			f = Pol(uniform!Base(source));

			// mask away bits above bit 53
			f.value &= Base((1L << 54) - 1);

			// set highest and lowest bit so that the degree is 53 and the
			// polynomial is not trivially reducible
			f.value |= (1L << 53) | 1;

			// test if f is irreducible
			if (f.irreducible)
				return f;
		}

		// If this is reached, we haven't found an irreducible polynomial in
		// randPolMaxTries. This error is very unlikely to occur.
		throw new Exception("unable to find new random irreducible polynomial");
	}

	/// Computes the Greatest Common Divisor of `x` and `f`.
	public static Pol gcd(Pol x, Pol f)
	{
		if (f.value == 0)
			return x;

		if (x.value == 0)
			return f;

		if (x.deg < f.deg)
		{
			import std.algorithm.mutation : swap;
			swap(x, f);
		}

		return gcd(f, x % f);
	}

	unittest
	{
		struct Test
		{
			private Pol f1;
			private Pol f2;
			private Pol gcd;
		}
		Test[] tests =
		[
			{Pol(10), Pol(50), Pol(2)},
			{Pol(0), Pol(1), Pol(1)},
			{
				parseBin("101101001"),
				parseBin("1010"),
				parseBin("1"),
			},
			{Pol(2), Pol(2), Pol(2)},
			{
				parseBin("1010"),
				parseBin("11"),
				parseBin("11"),
			},
			{
				Pol(0x8000000000000000),
				Pol(0x8000000000000000),
				Pol(0x8000000000000000),
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
				Pol(0x3DA3358B4DC173),
				Pol(0x3DA3358B4DC173),
				Pol(0x3DA3358B4DC173),
			},
			{
				Pol(0x3DA3358B4DC173),
				Pol(0x230d2259defd),
				Pol(1),
			},
			{
				Pol(0x230d2259defd),
				Pol(0x51b492b3eff2),
				parseBin("10011"),
			},
		];

		foreach (i, test; tests)
		{
			auto r = gcd(test.f1, test.f2);
			assert(test.gcd == r,
				format!"GCD test %d (%+s) failed: got %s, wanted %s"
				(i, test, r, test.gcd));

			r = gcd(test.f2, test.f1);
			assert(test.gcd == r,
				format!"GCD test %d (%+s) failed: got %s, wanted %s"
				(i, test, r, test.gcd));
		}
	}


	/// Returns `true` iff `x` is irreducible over `F_2`.
	/// This function uses Ben Or's reducibility test.
	//
	/// For details see "Tests and Constructions of Irreducible
	/// Polynomials over Finite Fields".
	@property public bool irreducible() const
	{
		foreach (i; 1 .. this.deg/2 + 1)
			if (gcd(this, qp(uint(i), this)).value != 1)
				return false;

		return true;
	}

	private
	{
		struct IrredTest
		{
			private Pol f;
			private bool irred;
		}
		static immutable IrredTest[] irredTests =
		[
			{Pol(0x38f1e565e288df), false},
			{Pol(0x3DA3358B4DC173), true},
			{Pol(0x30a8295b9d5c91), false},
			{Pol(0x255f4350b962cb), false},
			{Pol(0x267f776110a235), false},
			{Pol(0x2f4dae10d41227), false},
			{Pol(0x2482734cacca49), true},
			{Pol(0x312daf4b284899), false},
			{Pol(0x29dfb6553d01d1), false},
			{Pol(0x3548245eb26257), false},
			{Pol(0x3199e7ef4211b3), false},
			{Pol(0x362f39017dae8b), false},
			{Pol(0x200d57aa6fdacb), false},
			{Pol(0x35e0a4efa1d275), false},
			{Pol(0x2ced55b026577f), false},
			{Pol(0x260b012010893d), false},
			{Pol(0x2df29cbcd59e9d), false},
			{Pol(0x3f2ac7488bd429), false},
			{Pol(0x3e5cb1711669fb), false},
			{Pol(0x226d8de57a9959), false},
			{Pol(0x3c8de80aaf5835), false},
			{Pol(0x2026a59efb219b), false},
			{Pol(0x39dfa4d13fb231), false},
			{Pol(0x3143d0464b3299), false},
		];
	}

	unittest
	{
		foreach (_, test; irredTests)
			assert(test.f.irreducible == test.irred,
				format!"Irreducibility test for Polynomial %s failed: got %s, wanted %s"
				(test.f, test.f.irreducible, test.irred));
	}

	version(benchmarkPolynomials) private static void benchmarkPolIrreducible()
	{
		// find first irreducible polynomial
		Pol pol;
		foreach (_, test; irredTests)
			if (test.irred)
			{
				pol = test.f;
				break;
			}

		foreach (i; 0 .. Benchmark.N)
			if (!pol.irreducible)
				assert(false, format!"Irreducibility test for Polynomial %s failed"(pol));
	}


	/// Computes `this*f mod g`.
	public Pol mulMod(Pol f, Pol g)
	{
		if (value == 0 || f.value == 0)
			return Pol(0);

		Pol res;
		foreach (i; 0 .. f.deg + 1)
			if ((f.value & (1L << uint(i))) > 0)
			{
				auto a = this;
				foreach (j; 0 .. i)
					a = (a * Pol(2)) % g;
				res = (res + a) % g;
			}

		return res;
	}

	unittest
	{
		struct Test
		{
			private Pol f1;
			private Pol f2;
			private Pol g;
			private Pol mod;
		}
		Test[] tests =
		[
			{
				Pol(0x1230),
				Pol(0x230),
				Pol(0x55),
				Pol(0x22),
			},
			{
				Pol(0x0eae8c07dbbb3026),
				Pol(0xd5d6db9de04771de),
				Pol(0xdd2bda3b77c9),
				Pol(0x425ae8595b7a),
			},
		];

		foreach (i, test; tests)
		{
			auto mod = test.f1.mulMod(test.f2, test.g);
			assert(mod == test.mod,
				format!"mulMod test %d (%+s) failed: got %s, wanted %s"
				(i, test, mod, test.mod));
		}
	}


	/// Computes the polynomial `(x^(2^p)-x) mod g`.
	/// This is needed for the reducibility test.
	static private Pol qp(uint p, Pol g)
	{
		auto num = (1L << p);
		auto i = 1;

		// start with x
		auto res = Pol(2);

		while (i < num)
		{
			// repeatedly square res
			res = res.mulMod(res, g);
			i *= 2;
		}

		// add x
		return (res + Pol(2)) % g;
	}
}

version (benchmarkPolynomials)
{
	import chunker.internal.benchmark;
	mixin BenchmarkThisModule;

	static foreach (name; __traits(allMembers, Pol))
		static if (name.length > 9 && name[0..9] == "benchmark")
			mixin(`private void ` ~ name ~ `() { Pol.` ~ name ~ `(); }`);
	version = test;
}
version(unittest) version = test;
version(test) import std.format : format;
version(test) import std.stdio : stderr;

