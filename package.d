/**
Package chunker implements Content Defined Chunking (CDC) based on a rolling
Rabin Checksum.

Choosing a Random Irreducible Polynomial

The function RandomPolynomial() returns a new random polynomial of degree 53
for use with the chunker. The degree 53 is chosen because it is the largest
prime below 64-8 = 56, so that the top 8 bits of an ulong can be used for
optimising calculations in the chunker.

A random polynomial is chosen selecting 64 random bits, masking away bits
64..54 and setting bit 53 to one (otherwise the polynomial is not of the
desired degree) and bit 0 to one (otherwise the polynomial is trivially
reducible), so that 51 bits are chosen at random.

This process is repeated until Irreducible() returns true, then this
polynomials is returned. If this doesn't happen after 1 million tries, the
function returns an error. The probability for selecting an irreducible
polynomial at random is about 7.5% ( (2^53-2)/53 / 2^51), so the probability
that no irreducible polynomial has been found after 100 tries is lower than
0.04%.

Verifying Irreducible Polynomials

During development the results have been verified using the computational
discrete algebra system GAP, which can be obtained from the website at
http://www.gap-system.org/.

For filtering a given list of polynomials in hexadecimal coefficient notation,
the following script can be used:

	# create x over F_2 = GF(2);
	x := Indeterminate(GF(2), "x");

	# test if polynomial is irreducible, i.e. the number of factors is one;
	IrredPoly := function (poly);
		return (Length(Factors(poly)) = 1);
	end;;

	# create a polynomial in x from the hexadecimal representation of the;
	# coefficients;
	Hex2Poly := function (s);
		return ValuePol(CoefficientsQadic(IntHexString(s), 2), x);
	end;;

	# list of candidates, in hex;
	candidates := [ "3DA3358B4DC173" ];

	# create real polynomials;
	L := List(candidates, Hex2Poly);

	# filter and display the list of irreducible polynomials contained in L;
	Display(Filtered(L, x -> (IrredPoly(x))));

All irreducible polynomials from the list are written to the output.

Background Literature

An introduction to Rabin Fingerprints/Checksums can be found in the following articles:

Michael O. Rabin (1981): "Fingerprinting by Random Polynomials"
http://www.xmailserver.org/rabin.pdf

Ross N. Williams (1993): "A Painless Guide to CRC Error Detection Algorithms"
http://www.zlib.net/crc_v3.txt

Andrei Z. Broder (1993): "Some Applications of Rabin's Fingerprinting Method"
http://www.xmailserver.org/rabin_apps.pdf

Shuhong Gao and Daniel Panario (1997): "Tests and Constructions of Irreducible Polynomials over Finite Fields"
http://www.math.clemson.edu/~sgao/papers/GP97a.pdf

Andrew Kadatch, Bob Jenkins (2007): "Everything we know about CRC but afraid to forget"
http://crcutil.googlecode.com/files/crc-doc.1.0.pdf

License:

Copyright 2014 Alexander Neumann. All rights reserved.
Use of this source code is governed by a BSD-style
license that can be found in the LICENSE file.
*/

module chunker;

import chunker.polynomials;

import std.stdio : File;

private enum kiB = 1024;
private enum miB = 1024 * kiB;

/// WindowSize is the size of the sliding window.
private enum windowSize = 64;

/// MinSize is the default minimal size of a chunk.
public enum MinSize = 512 * kiB;
/// MaxSize is the default maximal size of a chunk.
public enum MaxSize = 8 * miB;

private enum chunkerBufSize = 512 * kiB;


private struct tables
{
	private Pol[256] out_;
	private Pol[256] mod;
}

/// cache precomputed tables, these are read-only anyway
private struct cache
{
static:
	tables[Pol] entries;
	Object mutex;
}

static this()
{
	cache.mutex = new Object;
}

/// Chunk is one content-dependent chunk of bytes whose end was cut when the
/// Rabin Fingerprint had the value stored in Cut.
public struct Chunk
{
	public uint Start;
	public uint Length;
	public ulong Cut;
	public ubyte[] Data;
}

private struct chunkerState
{
	private ubyte[windowSize] window;
	private int wpos;

	private ubyte[] buf;
	private uint bpos;
	private uint bmax;

	private uint start;
	private uint count;
	private uint pos;

	private uint pre; // wait for this many bytes before start calculating an new chunk

	private ulong digest;
}

private struct chunkerConfig
{
	public uint MinSize, MaxSize;

	private Pol pol;
	private uint polShift;
	private .tables tables;
	private bool tablesInitialized;
	private ulong splitmask;

	private File rd;
	private bool closed;
}

/// Chunker splits content with Rabin Fingerprints.
public struct Chunker
{
	.chunkerConfig config;
	.chunkerState state;
}

/// SetAverageBits allows to control the frequency of chunk discovery:
/// the lower averageBits, the higher amount of chunks will be identified.
/// The default value is 20 bits, so chunks will be of 1MiB size on average.
public void SetAverageBits(/*this*/ Chunker* c, int averageBits) {
	c.config.splitmask = (1 << ulong(averageBits)) - 1;
}

/// New returns a new Chunker based on polynomial p that reads from rd.
public Chunker* New(File rd, Pol pol) {
	return NewWithBoundaries(rd, pol, MinSize, MaxSize);
}

/// NewWithBoundaries returns a new Chunker based on polynomial p that reads from
/// rd and custom min and max size boundaries.
public Chunker* NewWithBoundaries(File rd, Pol pol, uint min, uint max) {
	Chunker v = {
		state: {
			buf: new ubyte[chunkerBufSize],
		},
		config: {
			pol:       pol,
			rd:        rd,
			MinSize:   min,
			MaxSize:   max,
			splitmask: (1 << 20) - 1, // aim to create chunks of 20 bits or about 1MiB on average.
		},
	};
	auto c = [v].ptr;

	c.reset();

	return c;
}

/// Reset reinitializes the chunker with a new reader and polynomial.
public void Reset(/*this*/ Chunker* c, File rd, Pol pol) {
	c.ResetWithBoundaries(rd, pol, MinSize, MaxSize);
}

/// ResetWithBoundaries reinitializes the chunker with a new reader, polynomial
/// and custom min and max size boundaries.
public void ResetWithBoundaries(/*this*/ Chunker* c, File rd, Pol pol, uint min, uint max) {
	Chunker v = {
		state: {
			buf: new ubyte[chunkerBufSize],
		},
		config: {
			pol:       pol,
			rd:        rd,
			MinSize:   min,
			MaxSize:   max,
			splitmask: (1 << 20) - 1, // aim to create chunks of 20 bits or about 1MiB on average.
		},
	};
	*c = v;

	c.reset();
}

private void reset(/*this*/ Chunker* c) {
	c.config.polShift = uint(c.config.pol.Deg() - 8);
	c.fillTables();

	for (auto i = 0; i < windowSize; i++) {
		c.state.window[i] = 0;
	}

	c.config.closed = false;
	c.state.digest = 0;
	c.state.wpos = 0;
	c.state.count = 0;
	c.state.digest = c.slide(c.state.digest, 1);
	c.state.start = c.state.pos;

	// do not start a new chunk unless at least MinSize bytes have been read
	c.state.pre = c.config.MinSize - windowSize;
}

/// fillTables calculates out_table and mod_table for optimization. This
/// implementation uses a cache in the global variable cache.
private void fillTables(/*this*/ Chunker* c) {
	// if polynomial hasn't been specified, do not compute anything for now
	if (c.config.pol == 0) {
		return;
	}

	c.config.tablesInitialized = true;

	// test if the tables are cached for this polynomial
	synchronized(cache.mutex)
	{
		if (auto t = c.config.pol in cache.entries) {
			c.config.tables = *t;
			return;
		}

		// calculate table for sliding out bytes. The byte to slide out is used as
		// the index for the table, the value contains the following:
		// out_table[b] = Hash(b || 0 ||        ...        || 0)
		//                          \ windowsize-1 zero bytes /
		// To slide out byte b_0 for window size w with known hash
		// H := H(b_0 || ... || b_w), it is sufficient to add out_table[b_0]:
		//    H(b_0 || ... || b_w) + H(b_0 || 0 || ... || 0)
		//  = H(b_0 + b_0 || b_1 + 0 || ... || b_w + 0)
		//  = H(    0     || b_1 || ...     || b_w)
		//
		// Afterwards a new byte can be shifted in.
		for (auto b = 0; b < 256; b++) {
			Pol h;

			h = appendByte(h, cast(ubyte)b, c.config.pol);
			for (auto i = 0; i < windowSize-1; i++) {
				h = appendByte(h, 0, c.config.pol);
			}
			c.config.tables.out_[b] = h;
		}

		// calculate table for reduction mod Polynomial
		auto k = c.config.pol.Deg();
		for (auto b = 0; b < 256; b++) {
			// mod_table[b] = A | B, where A = (b(x) * x^k mod pol) and  B = b(x) * x^k
			//
			// The 8 bits above deg(Polynomial) determine what happens next and so
			// these bits are used as a lookup to this table. The value is split in
			// two parts: Part A contains the result of the modulus operation, part
			// B is used to cancel out the 8 top bits so that one XOR operation is
			// enough to reduce modulo Polynomial
			c.config.tables.mod[b] = Pol(ulong(b)<<uint(k)).Mod(c.config.pol) | (Pol(b) << uint(k));
		}

		cache.entries[c.config.pol] = c.config.tables;
	}
}

/// Next returns the position and length of the next chunk of data. If an error
/// occurs while reading, the error is returned. Afterwards, the state of the
/// current chunk is undefined. When the last chunk has been returned, all
/// subsequent calls yield an io.EOF error.
public Chunk Next(/*this*/ Chunker* c, ubyte[] data) {
	data = data[];
	if (!c.config.tablesInitialized) {
		throw new Exception("tables for polynomial computation not initialized");
	}

	auto tabout = c.config.tables.out_;
	auto tabmod = c.config.tables.mod;
	auto polShift = c.config.polShift;
	auto minSize = c.config.MinSize;
	auto maxSize = c.config.MaxSize;
	auto buf = c.state.buf;
	while (true) {
		if (c.state.bpos >= c.state.bmax) {
			auto n = c.config.rd.rawRead(buf[]).length;

			// There are no more bytes to buffer, so this was the last
			// chunk.
			if (n == 0 && !c.config.closed) {
				c.config.closed = true;

				// return current chunk, if any bytes have been processed
				if (c.state.count > 0) {
					return Chunk(
						/*Start: */ c.state.start,
						/*Length:*/ c.state.count,
						/*Cut:   */ c.state.digest,
						/*Data:  */ data,
					);
				}
			}

			if (n == 0)
				return Chunk.init;

			c.state.bpos = 0;
			c.state.bmax = cast(uint)n;
		}

		// check if bytes have to be dismissed before starting a new chunk
		if (c.state.pre > 0) {
			auto n = c.state.bmax - c.state.bpos;
			if (c.state.pre > uint(n)) {
				c.state.pre -= uint(n);
				data ~= buf[c.state.bpos .. c.state.bmax];

				c.state.count += uint(n);
				c.state.pos += uint(n);
				c.state.bpos = c.state.bmax;

				continue;
			}

			data ~= buf[c.state.bpos .. c.state.bpos+c.state.pre];

			c.state.bpos += c.state.pre;
			c.state.count += c.state.pre;
			c.state.pos += c.state.pre;
			c.state.pre = 0;
		}

		auto add = c.state.count;
		auto digest = c.state.digest;
		auto win = c.state.window;
		auto wpos = c.state.wpos;
		foreach (_, b; buf[c.state.bpos .. c.state.bmax]) {
			// slide(b)
			auto out_ = win[wpos];
			win[wpos] = b;
			digest ^= ulong(tabout[out_]);
			wpos++;
			if (wpos >= windowSize) {
				wpos = 0;
			}

			// updateDigest
			auto index = cast(ubyte)(digest >> polShift);
			digest <<= 8;
			digest |= ulong(b);

			digest ^= ulong(tabmod[index]);
			// end manual inline

			add++;
			if (add < minSize) {
				continue;
			}

			if ((digest&c.config.splitmask) == 0 || add >= maxSize) {
				auto i = add - c.state.count - 1;
				data ~= c.state.buf[c.state.bpos .. c.state.bpos+uint(i)+1];
				c.state.count = add;
				c.state.pos += uint(i) + 1;
				c.state.bpos += uint(i) + 1;
				c.state.buf = buf;

				auto chunk = Chunk(
					/*Start: */ c.state.start,
					/*Length:*/ c.state.count,
					/*Cut:   */ digest,
					/*Data:  */ data,
				);

				c.reset();

				return chunk;
			}
		}
		c.state.digest = digest;
		c.state.window = win;
		c.state.wpos = wpos;

		auto steps = c.state.bmax - c.state.bpos;
		if (steps > 0) {
			data ~= c.state.buf[c.state.bpos .. c.state.bpos+steps];
		}
		c.state.count += steps;
		c.state.pos += steps;
		c.state.bpos = c.state.bmax;
	}
}

private ulong /*newDigest*/ updateDigest(ulong digest, uint polShift, tables tab, ubyte b) {
	auto index = digest >> polShift;
	digest <<= 8;
	digest |= ulong(b);

	digest ^= ulong(tab.mod[index]);
	return digest;
}

private ulong /*newDigest*/ slide(/*this*/ Chunker* c, ulong digest, ubyte b) {
	auto out_ = c.state.window[c.state.wpos];
	c.state.window[c.state.wpos] = b;
	digest ^= ulong(c.config.tables.out_[out_]);
	c.state.wpos = (c.state.wpos + 1) % windowSize;

	digest = updateDigest(digest, c.config.polShift, c.config.tables, b);
	return digest;
}

private Pol appendByte(Pol hash, ubyte b, Pol pol) {
	hash <<= 8;
	hash |= Pol(b);

	return hash.Mod(pol);
}

