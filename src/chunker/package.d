/**
Package chunker implements Content Defined Chunking (CDC) based on a rolling
Rabin Checksum.

Choosing a Random Irreducible Polynomial
========================================

The function `Pol.getRandom()` returns a new random polynomial of degree 53
for use with the chunker. The degree 53 is chosen because it is the largest
prime below 64-8 = 56, so that the top 8 bits of an ulong can be used for
optimising calculations in the chunker.

A random polynomial is chosen selecting 64 random bits, masking away bits
64..54 and setting bit 53 to one (otherwise the polynomial is not of the
desired degree) and bit 0 to one (otherwise the polynomial is trivially
reducible), so that 51 bits are chosen at random.

This process is repeated until `Pol.irreducible` returns `true`, then this
polynomials is returned. If this doesn't happen after 1 million tries, the
function returns an error. The probability for selecting an irreducible
polynomial at random is about 7.5% ( (2^53-2)/53 / 2^51), so the probability
that no irreducible polynomial has been found after 100 tries is lower than
0.04%.

Verifying Irreducible Polynomials
=================================

During development the results have been verified using the computational
discrete algebra system GAP, which can be obtained from the website at
http://www.gap-system.org/.

For filtering a given list of polynomials in hexadecimal coefficient notation,
the following script can be used:

```
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
```

All irreducible polynomials from the list are written to the output.

Background Literature
=====================

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

private enum size_t kiB = 1024;
private enum size_t miB = 1024 * kiB;

/// Size of the sliding window.
private enum size_t windowSize = 64;

/// Default minimal size of a chunk.
public enum size_t minSize = 512 * kiB;
/// Default maximal size of a chunk.
public enum size_t maxSize = 8 * miB;

private enum chunkerBufSize = 512 * kiB;


private struct Tables
{
	private Pol[256] out_;
	private Pol[256] mod;
}

// cache precomputed tables, these are read-only anyway
private struct cache
{
static:
	Tables[Pol] entries;
	Object mutex;
}

static this()
{
	cache.mutex = new Object;
}

/// Splits content with Rabin Fingerprints.
struct Chunker(R)
{
	/// One content-dependent chunk of bytes whose end was cut when
	/// the Rabin Fingerprint had the value stored in `cut`.
	public struct Chunk
	{
		/// Offset within the source data for the start of the chunk.
		public size_t start;
		/// Length of the chunk.
		public size_t length;
		/// Value of the rolling hash when the chunk was cut.
		public ulong cut;
		/// Contents of the chunk.
		public ubyte[] data;
	}

	private struct State
	{
		/// Current contents of the sliding window.
		private ubyte[windowSize] window;
		/// Sliding window cursor.
		private size_t wpos;

		/// Buffer used to receive and keep read data in.
		private ubyte[] buf;
		/// Current read position within `buf`.
		private size_t bpos;
		/// Number of bytes of data in `buf`.
		private size_t bmax;

		/// Start offset of the current chunk.
		private size_t start;
		/// Number of bytes within the current chunk.
		private size_t count;
		/// Current position within the stream.
		private size_t pos;

		/// Skip this many bytes before calculating a new chunk.
		private size_t pre;

		/// Current hash digest value.
		private ulong digest;
	}

	private struct Config
	{
		/// Minimum and maximum chunk sizes, as configured.
		public size_t minSize, maxSize;

		/// Polynomial to use when calculating the hash.
		private Pol pol;
		/// Bits to shift the digest when updating the hash.
		/// Calculated from the polynomial's degree.
		private uint polShift;
		/// Precomputed tables used for hashing.
		private Tables tables;
		/// Whether `tables` have already been computed.
		private bool tablesInitialized;
		/// Hash mask used to decide chunk boundaries.
		/// By default `(1 << 20) - 1`, or configured from `setAverageBits`.
		private ulong splitmask;

		/// Input data source.
		private R rd;
		/// Whether we've read all data from the source.
		private bool closed;
	}

	Config config;
	State state;

	/// Allows to control the frequency of chunk discovery: the lower
	/// `averageBits`, the higher amount of chunks will be identified.
	/// The default value is 20 bits, so chunks will be of 1MiB size
	/// on average.
	public void setAverageBits(uint averageBits)
	{
		config.splitmask = (1L << averageBits) - 1;
	}

	/// Constructs a new `Chunker` based on polynomial `pol` that
	/// reads from `rd`.
	public this(R rd, Pol pol)
	{
		this(rd, pol, minSize, maxSize);
	}

	/// Constructs a new `Chunker` based on polynomial `pol` that
	/// reads from `rd` and custom `min` and `max` size boundaries.
	public this(R rd, Pol pol, size_t min, size_t max)
	{
		resetWithBoundaries(rd, pol, min, max);
	}

	/// Reinitializes the `Chunker` with a new reader and polynomial.
	public void reset(R rd, Pol pol)
	{
		resetWithBoundaries(rd, pol, minSize, maxSize);
	}

	/// Reinitializes the `Chunker` with a new reader, polynomial and
	/// custom min and max size boundaries.
	public void resetWithBoundaries(R rd, Pol pol, size_t min, size_t max)
	{
		state = State();
		state.buf = new ubyte[chunkerBufSize];
		config = Config();
		config.pol = pol;
		config.rd = rd;
		config.minSize = min;
		config.maxSize = max;
		config.splitmask = (1 << 20) - 1; // aim to create chunks of 20 bits or about 1MiB on average.
		reset();
	}

	private void reset()
	{
		config.polShift = uint(config.pol.deg() - 8);
		fillTables();

		foreach (i; 0 .. windowSize)
			state.window[i] = 0;

		config.closed = false;
		state.digest = 0;
		state.wpos = 0;
		state.count = 0;
		slide(1);
		state.start = state.pos;

		// do not start a new chunk unless at least minSize bytes have been read
		state.pre = config.minSize - windowSize;
	}

	/// `fillTables` calculates `out_table` and `mod_table` for
	/// optimization. This implementation uses a cache in the global
	/// variable `cache`.
	private void fillTables()
	{
		// if polynomial hasn't been specified, do not compute anything for now
		if (config.pol.value == 0)
			return;

		config.tablesInitialized = true;

		// test if the tables are cached for this polynomial
		synchronized(cache.mutex)
		{
			if (auto t = config.pol in cache.entries)
			{
				config.tables = *t;
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
			foreach (b; 0 .. 256)
			{
				Pol h;

				h = appendByte(h, cast(ubyte)b, config.pol);
				foreach (i; 0 .. windowSize-1)
					h = appendByte(h, 0, config.pol);
				config.tables.out_[b] = h;
			}

			// calculate table for reduction mod Polynomial
			auto k = config.pol.deg();
			foreach (b; 0 .. 256)
			{
				// mod_table[b] = A | B, where A = (b(x) * x^k mod pol) and  B = b(x) * x^k
				//
				// The 8 bits above deg(Polynomial) determine what happens next and so
				// these bits are used as a lookup to this table. The value is split in
				// two parts: Part A contains the result of the modulus operation, part
				// B is used to cancel out the 8 top bits so that one XOR operation is
				// enough to reduce modulo Polynomial
				config.tables.mod[b] = Pol(
					(Pol(Pol.Base(b) << uint(k)) % config.pol).value |
					    (Pol.Base(b) << uint(k)));
			}

			cache.entries[config.pol] = config.tables;
		}
	}

	/// Returns the position and length of the next chunk of data.
	/// If an exception occurs while reading, it is propagated.
	/// Afterwards, the state of the current chunk is undefined.
	/// When the last chunk has been returned, all subsequent calls
	/// return `Chunk.init`.
	public Chunk next(ubyte[] data)
	{
		data = data[0..0];
		if (!config.tablesInitialized)
			throw new Exception("tables for polynomial computation not initialized");

		auto tabout = &config.tables.out_;
		auto tabmod = &config.tables.mod;
		auto polShift = config.polShift;
		auto minSize = config.minSize;
		auto maxSize = config.maxSize;
		auto buf = state.buf;
		while (true)
		{
			if (state.bpos >= state.bmax)
			{
				auto n = config.rd.rawRead(buf[]).length;

				// There are no more bytes to buffer, so this was the last
				// chunk.
				if (n == 0 && !config.closed)
				{
					config.closed = true;

					// return current chunk, if any bytes have been processed
					if (state.count > 0)
						return Chunk
						(
							/*Start: */ state.start,
							/*Length:*/ state.count,
							/*Cut:   */ state.digest,
							/*Data:  */ data,
						);
				}

				if (n == 0)
					return Chunk.init;

				state.bpos = 0;
				state.bmax = n;
			}

			// check if bytes have to be dismissed before starting a new chunk
			if (state.pre > 0)
			{
				auto n = state.bmax - state.bpos;
				if (state.pre > n)
				{
					state.pre -= n;
					data ~= buf[state.bpos .. state.bmax];

					state.count += n;
					state.pos += n;
					state.bpos = state.bmax;

					continue;
				}

				data ~= buf[state.bpos .. state.bpos+state.pre];

				state.bpos += state.pre;
				state.count += state.pre;
				state.pos += state.pre;
				state.pre = 0;
			}

			auto add = state.count;
			auto digest = state.digest;
			auto window = state.window;
			auto wpos = state.wpos;
			foreach (_, b; buf[state.bpos .. state.bmax])
			{
				.slide(window, wpos, digest, *tabout, *tabmod, polShift, b);

				add++;
				if (add < minSize)
					continue;

				if ((digest&config.splitmask) == 0 || add >= maxSize)
				{
					auto i = add - state.count - 1;
					data ~= state.buf[state.bpos .. state.bpos+i+1];
					state.count = add;
					state.pos += i + 1;
					state.bpos += i + 1;
					state.buf = buf;

					auto chunk = Chunk
					(
						/*Start: */ state.start,
						/*Length:*/ state.count,
						/*Cut:   */ digest,
						/*Data:  */ data,
					);

					reset();

					return chunk;
				}
			}
			state.digest = digest;
			state.window = window;
			state.wpos = wpos;

			auto steps = state.bmax - state.bpos;
			if (steps > 0)
				data ~= state.buf[state.bpos .. state.bpos+steps];
			state.count += steps;
			state.pos += steps;
			state.bpos = state.bmax;
		}
	}

	private void slide(ubyte b)
	{
		.slide(state.window, state.wpos, state.digest, config.tables.out_, config.tables.mod, config.polShift, b);
	}
}

private void slide(ref ubyte[windowSize] window, ref size_t wpos, ref ulong digest, ref Pol[256] tabout, ref Pol[256] tabmod, uint polShift, ubyte b)
{
	auto out_ = window[wpos];
	window[wpos] = b;
	digest ^= ulong(tabout[out_].value);
	wpos++;
	if (wpos >= windowSize)
		wpos = 0;

	digest = updateDigest(digest, polShift, tabmod, b);
}

private ulong /*newDigest*/ updateDigest(ulong digest, uint polShift, ref Pol[256] tabmod, ubyte b)
{
	pragma(inline, true);
	auto index = cast(ubyte)(digest >> polShift);
	digest <<= 8;
	digest |= ulong(b);

	digest ^= ulong(tabmod[index].value);
	return digest;
}

private Pol appendByte(Pol hash, ubyte b, Pol pol)
{
	hash.value <<= 8;
	hash.value |= Pol.Base(b);

	return hash % pol;
}

// -----------------------------------------------------------------------------
package: // shared with example

private import chunker.internal.gorng;

package ubyte[] getRandom(int seed, int count)
{
	import std.random : Random, uniform;
	auto buf = new ubyte[count];

	RNGSource rnd;
	.seed(&rnd, seed);
	for (auto i = 0; i < count; i += 4)
	{
		auto r = int63(&rnd) >> 31;
		buf[i] = cast(ubyte)(r);
		buf[i+1] = cast(ubyte)(r >> 8);
		buf[i+2] = cast(ubyte)(r >> 16);
		buf[i+3] = cast(ubyte)(r >> 24);
	}

	return buf;
}

package ubyte[32] hashData(ubyte[] d)
{
	import std.digest.sha;
	return sha256Of(d);
}

// Temporary D shim
import std.stdio : File;
package File bufFile(ubyte[] buf)
{
	File("temp.bin", "wb").rawWrite(buf);
	return File("temp.bin", "rb");
}

// -----------------------------------------------------------------------------
version(unittest) version = test;
version(benchmarkChunker) version = test;
version(test):
private:

private template parseDigest(string s)
{
	import std.conv : hexString;
	enum ubyte[] parseDigest = cast(ubyte[])hexString!s;
}

private struct TestChunk
{
	public size_t length;
	public ulong cutFP;
	public ubyte[] digest;
}

/// polynomial used for all the tests below
private enum testPol = Pol(0x3DA3358B4DC173);

/// created for 32MB of random data out of math/rand's Uint32() seeded by
/// constant 23
//
/// chunking configuration:
/// window size 64, avg chunksize 1<<20, min chunksize 1<<19, max chunksize 1<<23
/// polynom 0x3DA3358B4DC173
TestChunk[] chunks1 =
[
	{2163460, 0x000b98d4cdf00000, parseDigest!"4b94cb2cf293855ea43bf766731c74969b91aa6bf3c078719aabdd19860d590d"},
	{643703, 0x000d4e8364d00000, parseDigest!"5727a63c0964f365ab8ed2ccf604912f2ea7be29759a2b53ede4d6841e397407"},
	{1528956, 0x0015a25c2ef00000, parseDigest!"a73759636a1e7a2758767791c69e81b69fb49236c6929e5d1b654e06e37674ba"},
	{1955808, 0x00102a8242e00000, parseDigest!"c955fb059409b25f07e5ae09defbbc2aadf117c97a3724e06ad4abd2787e6824"},
	{2222372, 0x00045da878000000, parseDigest!"6ba5e9f7e1b310722be3627716cf469be941f7f3e39a4c3bcefea492ec31ee56"},
	{2538687, 0x00198a8179900000, parseDigest!"8687937412f654b5cfe4a82b08f28393a0c040f77c6f95e26742c2fc4254bfde"},
	{609606, 0x001d4e8d17100000, parseDigest!"5da820742ff5feb3369112938d3095785487456f65a8efc4b96dac4be7ebb259"},
	{1205738, 0x000a7204dd600000, parseDigest!"cc70d8fad5472beb031b1aca356bcab86c7368f40faa24fe5f8922c6c268c299"},
	{959742, 0x00183e71e1400000, parseDigest!"4065bdd778f95676c92b38ac265d361f81bff17d76e5d9452cf985a2ea5a4e39"},
	{4036109, 0x001fec043c700000, parseDigest!"b9cf166e75200eb4993fc9b6e22300a6790c75e6b0fc8f3f29b68a752d42f275"},
	{1525894, 0x000b1574b1500000, parseDigest!"2f238180e4ca1f7520a05f3d6059233926341090f9236ce677690c1823eccab3"},
	{1352720, 0x00018965f2e00000, parseDigest!"afd12f13286a3901430de816e62b85cc62468c059295ce5888b76b3af9028d84"},
	{811884, 0x00155628aa100000, parseDigest!"42d0cdb1ee7c48e552705d18e061abb70ae7957027db8ae8db37ec756472a70a"},
	{1282314, 0x001909a0a1400000, parseDigest!"819721c2457426eb4f4c7565050c44c32076a56fa9b4515a1c7796441730eb58"},
	{1318021, 0x001cceb980000000, parseDigest!"842eb53543db55bacac5e25cb91e43cc2e310fe5f9acc1aee86bdf5e91389374"},
	{948640, 0x0011f7a470a00000, parseDigest!"b8e36bf7019bb96ac3fb7867659d2167d9d3b3148c09fe0de45850b8fe577185"},
	{645464, 0x00030ce2d9400000, parseDigest!"5584bd27982191c3329f01ed846bfd266e96548dfa87018f745c33cfc240211d"},
	{533758, 0x0004435c53c00000, parseDigest!"4da778a25b72a9a0d53529eccfe2e5865a789116cb1800f470d8df685a8ab05d"},
	{1128303, 0x0000c48517800000, parseDigest!"08c6b0b38095b348d80300f0be4c5184d2744a17147c2cba5cc4315abf4c048f"},
	{800374, 0x000968473f900000, parseDigest!"820284d2c8fd243429674c996d8eb8d3450cbc32421f43113e980f516282c7bf"},
	{2453512, 0x001e197c92600000, parseDigest!"5fa870ed107c67704258e5e50abe67509fb73562caf77caa843b5f243425d853"},
	{2651975, 0x000ae6c868000000, parseDigest!"181347d2bbec32bef77ad5e9001e6af80f6abcf3576549384d334ee00c1988d8"},
	{237392, 0x0000000000000001, parseDigest!"fcd567f5d866357a8e299fd5b2359bb2c8157c30395229c4e9b0a353944a7978"},
];

/// test if nullbytes are correctly split, even if length is a multiple of minSize.
TestChunk[] chunks2 =
[
	{minSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{minSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{minSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{minSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
];

/// the same as chunks1, but avg chunksize is 1<<19
TestChunk[] chunks3 =
[
	{1491586, 0x00023e586ea80000, parseDigest!"4c008237df602048039287427171cef568a6cb965d1b5ca28dc80504a24bb061"},
	{671874, 0x000b98d4cdf00000, parseDigest!"fa8a42321b90c3d4ce9dd850562b2fd0c0fe4bdd26cf01a24f22046a224225d3"},
	{643703, 0x000d4e8364d00000, parseDigest!"5727a63c0964f365ab8ed2ccf604912f2ea7be29759a2b53ede4d6841e397407"},
	{1284146, 0x0012b527e4780000, parseDigest!"16d04cafecbeae9eaedd49da14c7ad7cdc2b1cc8569e5c16c32c9fb045aa899a"},
	{823366, 0x000d1d6752180000, parseDigest!"48662c118514817825ad4761e8e2e5f28f9bd8281b07e95dcafc6d02e0aa45c3"},
	{810134, 0x0016071b6e180000, parseDigest!"f629581aa05562f97f2c359890734c8574c5575da32f9289c5ba70bfd05f3f46"},
	{567118, 0x00102a8242e00000, parseDigest!"d4f0797c56c60d01bac33bfd49957a4816b6c067fc155b026de8a214cab4d70a"},
	{821315, 0x001b3e42c8180000, parseDigest!"8ebd0fd5db0293bd19140da936eb8b1bbd3cd6ffbec487385b956790014751ca"},
	{1401057, 0x00045da878000000, parseDigest!"001360af59adf4871ef138cfa2bb49007e86edaf5ac2d6f0b3d3014510991848"},
	{2311122, 0x0005cbd885380000, parseDigest!"8276d489b566086d9da95dc5c5fe6fc7d72646dd3308ced6b5b6ddb8595f0aa1"},
	{608723, 0x001cfcd86f280000, parseDigest!"518db33ba6a79d4f3720946f3785c05b9611082586d47ea58390fc2f6de9449e"},
	{980456, 0x0013edb7a7f80000, parseDigest!"0121b1690738395e15fecba1410cd0bf13fde02225160cad148829f77e7b6c99"},
	{1140278, 0x0001f9f017e80000, parseDigest!"28ca7c74804b5075d4f5eeb11f0845d99f62e8ea3a42b9a05c7bd5f2fca619dd"},
	{2015542, 0x00097bf5d8180000, parseDigest!"6fe8291f427d48650a5f0f944305d3a2dbc649bd401d2655fc0bdd42e890ca5a"},
	{904752, 0x000e1863eff80000, parseDigest!"62af1f1eb3f588d18aff28473303cc4731fc3cafcc52ce818fee3c4c2820854d"},
	{713072, 0x001f3bb1b9b80000, parseDigest!"4bda9dc2e3031d004d87a5cc93fe5207c4b0843186481b8f31597dc6ffa1496c"},
	{675937, 0x001fec043c700000, parseDigest!"5299c8c5acec1b90bb020cd75718aab5e12abb9bf66291465fd10e6a823a8b4a"},
	{1525894, 0x000b1574b1500000, parseDigest!"2f238180e4ca1f7520a05f3d6059233926341090f9236ce677690c1823eccab3"},
	{1352720, 0x00018965f2e00000, parseDigest!"afd12f13286a3901430de816e62b85cc62468c059295ce5888b76b3af9028d84"},
	{811884, 0x00155628aa100000, parseDigest!"42d0cdb1ee7c48e552705d18e061abb70ae7957027db8ae8db37ec756472a70a"},
	{1282314, 0x001909a0a1400000, parseDigest!"819721c2457426eb4f4c7565050c44c32076a56fa9b4515a1c7796441730eb58"},
	{1093738, 0x0017f5d048880000, parseDigest!"5dddfa7a241b68f65d267744bdb082ee865f3c2f0d8b946ea0ee47868a01bbff"},
	{962003, 0x000b921f7ef80000, parseDigest!"0cb5c9ebba196b441c715c8d805f6e7143a81cd5b0d2c65c6aacf59ca9124af9"},
	{856384, 0x00030ce2d9400000, parseDigest!"7734b206d46f3f387e8661e81edf5b1a91ea681867beb5831c18aaa86632d7fb"},
	{533758, 0x0004435c53c00000, parseDigest!"4da778a25b72a9a0d53529eccfe2e5865a789116cb1800f470d8df685a8ab05d"},
	{1128303, 0x0000c48517800000, parseDigest!"08c6b0b38095b348d80300f0be4c5184d2744a17147c2cba5cc4315abf4c048f"},
	{800374, 0x000968473f900000, parseDigest!"820284d2c8fd243429674c996d8eb8d3450cbc32421f43113e980f516282c7bf"},
	{2453512, 0x001e197c92600000, parseDigest!"5fa870ed107c67704258e5e50abe67509fb73562caf77caa843b5f243425d853"},
	{665901, 0x00118c842cb80000, parseDigest!"deceec26163842fdef6560311c69bf8a9871a56e16d719e2c4b7e4d668ceb61f"},
	{1986074, 0x000ae6c868000000, parseDigest!"64cd64bf3c3bc389eb20df8310f0427d1c36ab2eaaf09e346bfa7f0453fc1a18"},
	{237392, 0x0000000000000001, parseDigest!"fcd567f5d866357a8e299fd5b2359bb2c8157c30395229c4e9b0a353944a7978"},
];

import std.format;

private Chunker!R.Chunk[] testWithData(R)(ref Chunker!R chnker, TestChunk[] testChunks, bool checkDigest)
{
	Chunker!R.Chunk[] chunks;

	size_t pos = 0;
	foreach (i, chunk; testChunks)
	{
		auto c = chnker.next(null);

		assert(c.start == pos,
			format!"Start for chunk %d does not match: expected %d, got %d"
			(i, pos, c.start));

		assert(c.length == chunk.length,
			format!"Length for chunk %d does not match: expected %d, got %d"
			(i, chunk.length, c.length));

		assert(c.cut == chunk.cutFP,
			format!"Cut fingerprint for chunk %d/%d does not match: expected %016x, got %016x"
			(i, testChunks.length-1, chunk.cutFP, c.cut));

		if (checkDigest)
		{
			auto digest = hashData(c.data);
			assert(chunk.digest == digest,
				format!"Digest fingerprint for chunk %d/%d does not match: expected %(%02x%), got %(%02x%)"
				(i, testChunks.length-1, chunk.digest, digest));
		}

		pos += c.length;
		chunks ~= c;
	}

	auto c = chnker.next(null);
	if (c !is typeof(c).init)
		assert(false, "Wrong error returned after last chunk");

	assert(chunks.length == testChunks.length,
		"Amounts of test and resulting chunks do not match");

	return chunks;
}

import std.array : replicate;

@(`Chunker`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = Chunker!File(bufFile(buf), testPol);
	testWithData!File(ch, chunks1, true);

	// setup nullbyte data source
	buf = replicate([ubyte(0)], chunks2.length*minSize);
	ch = Chunker!File(bufFile(buf), testPol);

	testWithData(ch, chunks2, true);
}

@(`ChunkerWithCustomAverageBits`) unittest
{
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = Chunker!File(bufFile(buf), testPol);

	// sligthly decrease averageBits to get more chunks
	ch.setAverageBits(19);
	testWithData(ch, chunks3, true);
}

@(`ChunkerReset`) unittest
{
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = Chunker!File(bufFile(buf), testPol);
	testWithData(ch, chunks1, true);

	ch.reset(bufFile(buf), testPol);
	testWithData(ch, chunks1, true);
}

import std.stdio : stderr;

@(`ChunkerWithRandomPolynomial`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	// generate a new random polynomial
	import std.datetime.systime : Clock, SysTime;
	auto start = Clock.currTime();
	auto p = Pol.getRandom();
	stderr.writefln!"generating random polynomial took %s"(Clock.currTime() - start);

	start = Clock.currTime();
	auto ch = Chunker!File(bufFile(buf), p);
	stderr.writefln!"creating chunker took %s"(Clock.currTime() - start);

	// make sure that first chunk is different
	auto c = ch.next(null);

	assert(c.cut != chunks1[0].cutFP,
		"Cut point is the same");

	assert(c.length != chunks1[0].length,
		"Length is the same");

	assert(hashData(c.data) != chunks1[0].digest,
		"Digest is the same");
}

@(`ChunkerWithoutHash`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	auto ch = Chunker!File(bufFile(buf), testPol);
	auto chunks = testWithData(ch, chunks1, false);

	// test reader
	foreach (i, c; chunks)
	{
		assert(c.data.length == chunks1[i].length,
			format!"reader returned wrong number of bytes: expected %d, got %d"
			(chunks1[i].length, c.data.length));

		assert(buf[c.start .. c.start+c.length] == c.data,
			format!"invalid data for chunk returned: expected %(%02x%), got %(%02x%)"
			(buf[c.start .. c.start+c.length], c.data));
	}

	// setup nullbyte data source
	buf = replicate([ubyte(0)], chunks2.length*minSize);
	ch = Chunker!File(bufFile(buf), testPol);

	testWithData(ch, chunks2, false);
}

version (benchmarkChunker)
{
	import chunker.internal.benchmark;
	mixin BenchmarkThisModule;

	void _benchmarkChunker(bool checkDigest)
	{
		auto size = 32 * 1024 * 1024;
		auto rd = bufFile(getRandom(23, size));
		auto ch = Chunker!File(rd, testPol);
		auto buf = new ubyte[maxSize];

		// b.SetBytes(long(size));

		int chunks;
		Benchmark.benchmark({
			chunks = 0;

			rd.seek(0);

			ch.reset(rd, testPol);

			auto cur = 0;
			while (true)
			{
				auto chunk = ch.next(buf);

				if (chunk is typeof(chunk).init)
				{
					break;
				}

				assert(chunk.length == chunks1[cur].length,
					format!"wrong chunk length, want %d, got %d"
					(chunks1[cur].length, chunk.length));

				assert(chunk.cut == chunks1[cur].cutFP,
					format!"wrong cut fingerprint, want 0x%x, got 0x%x"
					(chunks1[cur].cutFP, chunk.cut));

				if (checkDigest)
				{
					auto h = hashData(chunk.data);
					assert(h == chunks1[cur].digest,
						format!"wrong digest, want %(%02x%), got %(%02x%)"
						(chunks1[cur].digest, h));
				}

				chunks++;
				cur++;
			}
		});

		stderr.writefln!"%d chunks, average chunk size: %d bytes"(chunks, size/chunks);
	}

	void benchmarkChunkerWithSHA256()
	{
		_benchmarkChunker(true);
	}

	void benchmarkChunker()
	{
		_benchmarkChunker(false);
	}

	void benchmarkNewChunker()
	{
		auto p = Pol.getRandom();

		Benchmark.benchmark({
			Chunker!File(bufFile(null), p);
		});
	}
}
