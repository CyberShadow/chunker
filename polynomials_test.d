module chunker.polynomials_test;

import chunker.polynomials;

version(unittest) import std.format : format;
version(unittest) import std.stdio : stderr;

private struct polAddTest {
	private Pol x, y;
	private Pol sum;
}
polAddTest[] polAddTests = [
	{23, 16, 23 ^ 16},
	{0x9a7e30d1e855e0a0, 0x670102a1f4bcd414, 0xfd7f32701ce934b4},
	{0x9a7e30d1e855e0a0, 0x9a7e30d1e855e0a0, 0},
];

@(`PolAdd`) unittest {
	foreach (i, test; polAddTests) {
		if (test.sum != test.x.Add(test.y)) {
			assert(false, format!"test %d failed: sum != x+y"( i));
		}

		if (test.sum != test.y.Add(test.x)) {
			assert(false, format!"test %d failed: sum != y+x"( i));
		}
	}
}

private Pol parseBin(string s) {
	import std.conv : to;
	auto i = s.to!Pol(2);

	return Pol(i);
}

private struct polMulTest {
	private Pol x, y;
	private Pol res;
};
polMulTest[] polMulTests = [
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

@(`PolMul`) unittest {
	foreach (i, test; polMulTests) {
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

@(`PolMulOverflow`) unittest {
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

private struct polDivTest {
	private Pol x, y;
	private Pol res;
}
polDivTest[] polDivTests = [
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

@(`PolDiv`) unittest {
	foreach (i, test; polDivTests) {
		auto m = test.x.Div(test.y);
		if (test.res != m) {
			assert(false, format!"TestPolDiv failed for test %d: %s * %s: want %s, got %s"(
				i, test.x, test.y, test.res, m));
		}
	}
}

@(`PolDeg`) unittest {
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

private struct polModTest {
	private Pol x, y;
	private Pol res;
}
polModTest[] polModTests = [
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

@(`PolModt`) unittest {
	foreach (i, test; polModTests) {
		auto res = test.x.Mod(test.y);
		if (test.res != res) {
			assert(false, format!"test %d failed: want %s, got %s"( i, test.res, res));
		}
	}
}

version = benchmark;
version (benchmark) enum N = 1;

version(benchmark) @(`BenchmarkPolDivMod`) unittest {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.DivMod(f);
	}
}

version(benchmark) @(`BenchmarkPolDiv`) unittest {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.Div(f);
	}
}

version(benchmark) @(`BenchmarkPolMod`) unittest {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < N; i++) {
		g.Mod(f);
	}
}

version(benchmark) @(`BenchmarkPolDeg`) unittest {
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

@(`RandomPolynomial`) unittest {
	RandomPolynomial();
}

version(benchmark) @(`BenchmarkRandomPolynomial`) unittest {
	for (auto i = 0; i < N; i++) {
		RandomPolynomial();
	}
}

@(`ExpandPolynomial`) unittest {
	auto pol = Pol(0x3DA3358B4DC173);
	auto s = pol.Expand();
	if (s != "x^53+x^52+x^51+x^50+x^48+x^47+x^45+x^41+x^40+x^37+x^36+x^34+x^32+x^31+x^27+x^25+x^24+x^22+x^19+x^18+x^16+x^15+x^14+x^8+x^6+x^5+x^4+x+1") {
		assert(false, "wrong result");
	}
}

private struct polIrredTest {
	private Pol f;
	private bool irred;
}
polIrredTest[] polIrredTests = [
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

@(`PolIrreducible`) unittest {
	foreach (_, test; polIrredTests) {
		if (test.f.Irreducible() != test.irred) {
			assert(false, format!"Irreducibility test for Polynomial %s failed: got %s, wanted %s"(
				test.f, test.f.Irreducible(), test.irred));
		}
	}
}

version(benchmark) @(`BenchmarkPolIrreducible`) unittest {
	// find first irreducible polynomial
	Pol pol;
	foreach (_, test; polIrredTests) {
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

private struct polGCDTest {
	private Pol f1;
	private Pol f2;
	private Pol gcd;
}
polGCDTest[] polGCDTests = [
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

@(`PolGCD`) unittest {
	foreach (i, test; polGCDTests) {
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

private struct polMulModTest {
	private Pol f1;
	private Pol f2;
	private Pol g;
	private Pol mod;
}
polMulModTest[] polMulModTests = [
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

@(`PolMulMod`) unittest {
	foreach (i, test; polMulModTests) {
		auto mod = test.f1.MulMod(test.f2, test.g);
		if (mod != test.mod) {
			assert(false, format!"MulMod test %d (%+s) failed: got %s, wanted %s"(
				i, test, mod, test.mod));
		}
	}
}
