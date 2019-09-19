module chunker.polynomials;

/// Pol is a polynomial from F_2[X].
alias Pol = ulong;

/// Add returns x+y.
public Pol Add(/*this*/ Pol x, Pol y) {
	auto r = Pol(ulong(x) ^ ulong(y));
	return r;
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
		if ((y & (1 << uint(i))) > 0) {
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
		if ((x&(1<<uint(i))) > 0) {
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
		q |= (1 << uint(diff));
		x = x.Add(m);

		diff = x.Deg() - D;
	}

	return [q, x];
}

/// Div returns the integer division result x / d.
public Pol Div(/*this*/ Pol x, Pol d) {
	auto q = x.DivMod(d)[0];
	return q;
}

/// Mod returns the remainder of x / d
public Pol Mod(/*this*/ Pol x, Pol d) {
	auto r = x.DivMod(d)[1];
	return r;
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

/// MulMod computes x*f mod g
public Pol MulMod(/*this*/ Pol x, Pol f, Pol g) {
	if (x == 0 || f == 0) {
		return 0;
	}

	Pol res;
	for (auto i = 0; i <= f.Deg(); i++) {
		if ((f & (1 << uint(i))) > 0) {
			auto a = x;
			for (auto j = 0; j < i; j++) {
				a = a.Mul(2).Mod(g);
			}
			res = res.Add(a).Mod(g);
		}
	}

	return res;
}

/// qp computes the polynomial (x^(2^p)-x) mod g. This is needed for the
/// reducibility test.
private Pol qp(uint p, Pol g) {
	auto num = (1 << p);
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

