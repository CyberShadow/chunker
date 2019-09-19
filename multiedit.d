// ----------------------------------------------------------- chunker.d
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
	private tables tables;
	private bool tablesInitialized;
	private ulong splitmask;

	private io.Reader rd;
	private bool closed;
}

/// Chunker splits content with Rabin Fingerprints.
public struct Chunker
{
	chunkerConfig chunkerConfig;
	chunkerState chunkerState;
}

/// SetAverageBits allows to control the frequency of chunk discovery:
/// the lower averageBits, the higher amount of chunks will be identified.
/// The default value is 20 bits, so chunks will be of 1MiB size on average.
public void SetAverageBits(/*this*/ Chunker* c, int averageBits) {
	c.splitmask = (1 << ulong(averageBits)) - 1;
}

/// New returns a new Chunker based on polynomial p that reads from rd.
public Chunker* New(io.Reader rd, Pol pol) {
	return NewWithBoundaries(rd, pol, MinSize, MaxSize);
}

/// NewWithBoundaries returns a new Chunker based on polynomial p that reads from
/// rd and custom min and max size boundaries.
public Chunker* NewWithBoundaries(io.Reader rd, Pol pol, uint min, uint max) {
	Chunker c = {
		chunkerState: {
			buf: new ubyte[chunkerBufSize],
		},
		chunkerConfig: {
			pol:       pol,
			rd:        rd,
			MinSize:   min,
			MaxSize:   max,
			splitmask: (1 << 20) - 1, // aim to create chunks of 20 bits or about 1MiB on average.
		},
	};

	c.reset();

	return [c].ptr;
}

/// Reset reinitializes the chunker with a new reader and polynomial.
public void Reset(/*this*/ Chunker* c, io.Reader rd, Pol pol) {
	c.ResetWithBoundaries(rd, pol, MinSize, MaxSize);
}

/// ResetWithBoundaries reinitializes the chunker with a new reader, polynomial
/// and custom min and max size boundaries.
public void ResetWithBoundaries(/*this*/ Chunker* c, io.Reader rd, Pol pol, uint min, uint max) {
	Chunker v = {
		chunkerState: {
			buf: new ubyte[chunkerBufSize],
		},
		chunkerConfig: {
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
	c.polShift = uint(c.pol.Deg() - 8);
	c.fillTables();

	for (auto i = 0; i < windowSize; i++) {
		c.window[i] = 0;
	}

	c.closed = false;
	c.digest = 0;
	c.wpos = 0;
	c.count = 0;
	c.digest = c.slide(c.digest, 1);
	c.start = c.pos;

	// do not start a new chunk unless at least MinSize bytes have been read
	c.pre = c.MinSize - windowSize;
}

/// fillTables calculates out_table and mod_table for optimization. This
/// implementation uses a cache in the global variable cache.
private void fillTables(/*this*/ Chunker* c) {
	// if polynomial hasn't been specified, do not compute anything for now
	if (c.pol == 0) {
		return;
	}

	c.tablesInitialized = true;

	// test if the tables are cached for this polynomial
	synchronized(cache.mutex)
	{
		if (auto t = c.pol in cache.entries) {
			c.tables = *t;
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

			h = appendByte(h, ubyte(b), c.pol);
			for (auto i = 0; i < windowSize-1; i++) {
				h = appendByte(h, 0, c.pol);
			}
			c.tables.out_[b] = h;
		}

		// calculate table for reduction mod Polynomial
		auto k = c.pol.Deg();
		for (auto b = 0; b < 256; b++) {
			// mod_table[b] = A | B, where A = (b(x) * x^k mod pol) and  B = b(x) * x^k
			//
			// The 8 bits above deg(Polynomial) determine what happens next and so
			// these bits are used as a lookup to this table. The value is split in
			// two parts: Part A contains the result of the modulus operation, part
			// B is used to cancel out the 8 top bits so that one XOR operation is
			// enough to reduce modulo Polynomial
			c.tables.mod[b] = Pol(ulong(b)<<uint(k)).Mod(c.pol) | (Pol(b) << uint(k));
		}

		cache.entries[c.pol] = c.tables;
	}
}

/// Next returns the position and length of the next chunk of data. If an error
/// occurs while reading, the error is returned. Afterwards, the state of the
/// current chunk is undefined. When the last chunk has been returned, all
/// subsequent calls yield an io.EOF error.
public void Next(/*this*/ Chunker* c, ubyte[] data) (Chunk, error) {
	data = data[];
	if (!c.tablesInitialized) {
		throw new Exception("tables for polynomial computation not initialized");
	}

	auto tabout = c.tables.out_;
	auto tabmod = c.tables.mod;
	auto polShift = c.polShift;
	auto minSize = c.MinSize;
	auto maxSize = c.MaxSize;
	auto buf = c.buf;
	while (true) {
		if (c.bpos >= c.bmax) {
			n, err := io.ReadFull(c.rd, buf[]);

			if (err == io.ErrUnexpectedEOF) {
				err = nil;
			}

			// io.ReadFull only returns io.EOF when no bytes could be read. If
			// this is the case and we're in this branch, there are no more
			// bytes to buffer, so this was the last chunk. If a different
			// error has occurred, return that error and abandon the current
			// chunk.
			if (err == io.EOF && !c.closed) {
				c.closed = true;

				// return current chunk, if any bytes have been processed
				if (c.count > 0) {
					return Chunk{
						Start:  c.start,
						Length: c.count,
						Cut:    c.digest,
						Data:   data,
					}, nil;
				}
			}

			if (err != nil) {
				return Chunk{}, err;
			}

			c.bpos = 0;
			c.bmax = uint(n);
		}

		// check if bytes have to be dismissed before starting a new chunk
		if (c.pre > 0) {
			auto n = c.bmax - c.bpos;
			if (c.pre > uint(n)) {
				c.pre -= uint(n);
				data ~= buf[c.bpos .. c.bmax];

				c.count += uint(n);
				c.pos += uint(n);
				c.bpos = c.bmax;

				continue;
			}

			data ~= buf[c.bpos .. c.bpos+c.pre];

			c.bpos += c.pre;
			c.count += c.pre;
			c.pos += c.pre;
			c.pre = 0;
		}

		auto add = c.count;
		auto digest = c.digest;
		auto win = c.window;
		auto wpos = c.wpos;
		foreach (_, b; buf[c.bpos .. c.bmax]) {
			// slide(b)
			auto out_ = win[wpos];
			win[wpos] = b;
			digest ^= ulong(tabout[out_]);
			wpos++;
			if (wpos >= windowSize) {
				wpos = 0;
			}

			// updateDigest
			auto index = ubyte(digest >> polShift);
			digest <<= 8;
			digest |= ulong(b);

			digest ^= ulong(tabmod[index]);
			// end manual inline

			add++;
			if (add < minSize) {
				continue;
			}

			if ((digest&c.splitmask) == 0 || add >= maxSize) {
				auto i = add - c.count - 1;
				data ~= c.buf[c.bpos .. c.bpos+uint(i)+1];
				c.count = add;
				c.pos += uint(i) + 1;
				c.bpos += uint(i) + 1;
				c.buf = buf;

				auto chunk = Chunk{
					Start:  c.start,
					Length: c.count,
					Cut:    digest,
					Data:   data,
				}

				c.reset();

				return chunk, nil;
			}
		}
		c.digest = digest;
		c.window = win;
		c.wpos = wpos;

		auto steps = c.bmax - c.bpos;
		if (steps > 0) {
			data ~= c.buf[c.bpos .. c.bpos+steps];
		}
		c.count += steps;
		c.pos += steps;
		c.bpos = c.bmax;
	}
}

private (ulong newDigest) updateDigest(ulong digest, uint polShift, tables tab, ubyte b) {
	auto index = digest >> polShift;
	digest <<= 8;
	digest |= ulong(b);

	digest ^= ulong(tab.mod[index]);
	return digest;
}

private void slide(/*this*/ Chunker* c, ulong digest, ubyte b) (ulong newDigest) {
	auto out_ = c.window[c.wpos];
	c.window[c.wpos] = b;
	digest ^= ulong(c.tables.out_[out_]);
	c.wpos = (c.wpos + 1) % windowSize;

	digest = updateDigest(digest, c.polShift, c.tables, b);
	return digest;
}

private Pol appendByte(Pol hash, ubyte b, Pol pol) {
	hash <<= 8;
	hash |= Pol(b);

	return hash.Mod(pol);
}

// ----------------------------------------------------------- chunker_test.d
module chunker.chunker_test;

private ubyte[] parseDigest(string s) {
	d, err := hex.DecodeString(s);
	if (err != nil) {
		panic(err);
	}

	return d;
}

private struct chunk
{
	public uint Length;
	public ulong CutFP;
	public ubyte[] Digest;
}

/// polynomial used for all the tests below
private enum testPol = Pol(0x3DA3358B4DC173);

/// created for 32MB of random data out of math/rand's Uint32() seeded by
/// constant 23
//
/// chunking configuration:
/// window size 64, avg chunksize 1<<20, min chunksize 1<<19, max chunksize 1<<23
/// polynom 0x3DA3358B4DC173
var chunks1 = chunk[]{
	chunk{2163460, 0x000b98d4cdf00000, parseDigest("4b94cb2cf293855ea43bf766731c74969b91aa6bf3c078719aabdd19860d590d")},
	chunk{643703, 0x000d4e8364d00000, parseDigest("5727a63c0964f365ab8ed2ccf604912f2ea7be29759a2b53ede4d6841e397407")},
	chunk{1528956, 0x0015a25c2ef00000, parseDigest("a73759636a1e7a2758767791c69e81b69fb49236c6929e5d1b654e06e37674ba")},
	chunk{1955808, 0x00102a8242e00000, parseDigest("c955fb059409b25f07e5ae09defbbc2aadf117c97a3724e06ad4abd2787e6824")},
	chunk{2222372, 0x00045da878000000, parseDigest("6ba5e9f7e1b310722be3627716cf469be941f7f3e39a4c3bcefea492ec31ee56")},
	chunk{2538687, 0x00198a8179900000, parseDigest("8687937412f654b5cfe4a82b08f28393a0c040f77c6f95e26742c2fc4254bfde")},
	chunk{609606, 0x001d4e8d17100000, parseDigest("5da820742ff5feb3369112938d3095785487456f65a8efc4b96dac4be7ebb259")},
	chunk{1205738, 0x000a7204dd600000, parseDigest("cc70d8fad5472beb031b1aca356bcab86c7368f40faa24fe5f8922c6c268c299")},
	chunk{959742, 0x00183e71e1400000, parseDigest("4065bdd778f95676c92b38ac265d361f81bff17d76e5d9452cf985a2ea5a4e39")},
	chunk{4036109, 0x001fec043c700000, parseDigest("b9cf166e75200eb4993fc9b6e22300a6790c75e6b0fc8f3f29b68a752d42f275")},
	chunk{1525894, 0x000b1574b1500000, parseDigest("2f238180e4ca1f7520a05f3d6059233926341090f9236ce677690c1823eccab3")},
	chunk{1352720, 0x00018965f2e00000, parseDigest("afd12f13286a3901430de816e62b85cc62468c059295ce5888b76b3af9028d84")},
	chunk{811884, 0x00155628aa100000, parseDigest("42d0cdb1ee7c48e552705d18e061abb70ae7957027db8ae8db37ec756472a70a")},
	chunk{1282314, 0x001909a0a1400000, parseDigest("819721c2457426eb4f4c7565050c44c32076a56fa9b4515a1c7796441730eb58")},
	chunk{1318021, 0x001cceb980000000, parseDigest("842eb53543db55bacac5e25cb91e43cc2e310fe5f9acc1aee86bdf5e91389374")},
	chunk{948640, 0x0011f7a470a00000, parseDigest("b8e36bf7019bb96ac3fb7867659d2167d9d3b3148c09fe0de45850b8fe577185")},
	chunk{645464, 0x00030ce2d9400000, parseDigest("5584bd27982191c3329f01ed846bfd266e96548dfa87018f745c33cfc240211d")},
	chunk{533758, 0x0004435c53c00000, parseDigest("4da778a25b72a9a0d53529eccfe2e5865a789116cb1800f470d8df685a8ab05d")},
	chunk{1128303, 0x0000c48517800000, parseDigest("08c6b0b38095b348d80300f0be4c5184d2744a17147c2cba5cc4315abf4c048f")},
	chunk{800374, 0x000968473f900000, parseDigest("820284d2c8fd243429674c996d8eb8d3450cbc32421f43113e980f516282c7bf")},
	chunk{2453512, 0x001e197c92600000, parseDigest("5fa870ed107c67704258e5e50abe67509fb73562caf77caa843b5f243425d853")},
	chunk{2651975, 0x000ae6c868000000, parseDigest("181347d2bbec32bef77ad5e9001e6af80f6abcf3576549384d334ee00c1988d8")},
	chunk{237392, 0x0000000000000001, parseDigest("fcd567f5d866357a8e299fd5b2359bb2c8157c30395229c4e9b0a353944a7978")},
}

/// test if nullbytes are correctly split, even if length is a multiple of MinSize.
var chunks2 = chunk[]{
	chunk{MinSize, 0, parseDigest("07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541")},
	chunk{MinSize, 0, parseDigest("07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541")},
	chunk{MinSize, 0, parseDigest("07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541")},
	chunk{MinSize, 0, parseDigest("07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541")},
}

/// the same as chunks1, but avg chunksize is 1<<19
var chunks3 = chunk[]{
	chunk{1491586, 0x00023e586ea80000, parseDigest("4c008237df602048039287427171cef568a6cb965d1b5ca28dc80504a24bb061")},
	chunk{671874, 0x000b98d4cdf00000, parseDigest("fa8a42321b90c3d4ce9dd850562b2fd0c0fe4bdd26cf01a24f22046a224225d3")},
	chunk{643703, 0x000d4e8364d00000, parseDigest("5727a63c0964f365ab8ed2ccf604912f2ea7be29759a2b53ede4d6841e397407")},
	chunk{1284146, 0x0012b527e4780000, parseDigest("16d04cafecbeae9eaedd49da14c7ad7cdc2b1cc8569e5c16c32c9fb045aa899a")},
	chunk{823366, 0x000d1d6752180000, parseDigest("48662c118514817825ad4761e8e2e5f28f9bd8281b07e95dcafc6d02e0aa45c3")},
	chunk{810134, 0x0016071b6e180000, parseDigest("f629581aa05562f97f2c359890734c8574c5575da32f9289c5ba70bfd05f3f46")},
	chunk{567118, 0x00102a8242e00000, parseDigest("d4f0797c56c60d01bac33bfd49957a4816b6c067fc155b026de8a214cab4d70a")},
	chunk{821315, 0x001b3e42c8180000, parseDigest("8ebd0fd5db0293bd19140da936eb8b1bbd3cd6ffbec487385b956790014751ca")},
	chunk{1401057, 0x00045da878000000, parseDigest("001360af59adf4871ef138cfa2bb49007e86edaf5ac2d6f0b3d3014510991848")},
	chunk{2311122, 0x0005cbd885380000, parseDigest("8276d489b566086d9da95dc5c5fe6fc7d72646dd3308ced6b5b6ddb8595f0aa1")},
	chunk{608723, 0x001cfcd86f280000, parseDigest("518db33ba6a79d4f3720946f3785c05b9611082586d47ea58390fc2f6de9449e")},
	chunk{980456, 0x0013edb7a7f80000, parseDigest("0121b1690738395e15fecba1410cd0bf13fde02225160cad148829f77e7b6c99")},
	chunk{1140278, 0x0001f9f017e80000, parseDigest("28ca7c74804b5075d4f5eeb11f0845d99f62e8ea3a42b9a05c7bd5f2fca619dd")},
	chunk{2015542, 0x00097bf5d8180000, parseDigest("6fe8291f427d48650a5f0f944305d3a2dbc649bd401d2655fc0bdd42e890ca5a")},
	chunk{904752, 0x000e1863eff80000, parseDigest("62af1f1eb3f588d18aff28473303cc4731fc3cafcc52ce818fee3c4c2820854d")},
	chunk{713072, 0x001f3bb1b9b80000, parseDigest("4bda9dc2e3031d004d87a5cc93fe5207c4b0843186481b8f31597dc6ffa1496c")},
	chunk{675937, 0x001fec043c700000, parseDigest("5299c8c5acec1b90bb020cd75718aab5e12abb9bf66291465fd10e6a823a8b4a")},
	chunk{1525894, 0x000b1574b1500000, parseDigest("2f238180e4ca1f7520a05f3d6059233926341090f9236ce677690c1823eccab3")},
	chunk{1352720, 0x00018965f2e00000, parseDigest("afd12f13286a3901430de816e62b85cc62468c059295ce5888b76b3af9028d84")},
	chunk{811884, 0x00155628aa100000, parseDigest("42d0cdb1ee7c48e552705d18e061abb70ae7957027db8ae8db37ec756472a70a")},
	chunk{1282314, 0x001909a0a1400000, parseDigest("819721c2457426eb4f4c7565050c44c32076a56fa9b4515a1c7796441730eb58")},
	chunk{1093738, 0x0017f5d048880000, parseDigest("5dddfa7a241b68f65d267744bdb082ee865f3c2f0d8b946ea0ee47868a01bbff")},
	chunk{962003, 0x000b921f7ef80000, parseDigest("0cb5c9ebba196b441c715c8d805f6e7143a81cd5b0d2c65c6aacf59ca9124af9")},
	chunk{856384, 0x00030ce2d9400000, parseDigest("7734b206d46f3f387e8661e81edf5b1a91ea681867beb5831c18aaa86632d7fb")},
	chunk{533758, 0x0004435c53c00000, parseDigest("4da778a25b72a9a0d53529eccfe2e5865a789116cb1800f470d8df685a8ab05d")},
	chunk{1128303, 0x0000c48517800000, parseDigest("08c6b0b38095b348d80300f0be4c5184d2744a17147c2cba5cc4315abf4c048f")},
	chunk{800374, 0x000968473f900000, parseDigest("820284d2c8fd243429674c996d8eb8d3450cbc32421f43113e980f516282c7bf")},
	chunk{2453512, 0x001e197c92600000, parseDigest("5fa870ed107c67704258e5e50abe67509fb73562caf77caa843b5f243425d853")},
	chunk{665901, 0x00118c842cb80000, parseDigest("deceec26163842fdef6560311c69bf8a9871a56e16d719e2c4b7e4d668ceb61f")},
	chunk{1986074, 0x000ae6c868000000, parseDigest("64cd64bf3c3bc389eb20df8310f0427d1c36ab2eaaf09e346bfa7f0453fc1a18")},
	chunk{237392, 0x0000000000000001, parseDigest("fcd567f5d866357a8e299fd5b2359bb2c8157c30395229c4e9b0a353944a7978")},
}

private Chunk[] testWithData(testing*.T t, Chunker* chnker, chunk[] testChunks, bool checkDigest) {
	auto chunks = Chunk[]{};

	auto pos = uint(0);
	foreach (i, chunk; testChunks) {
		c, err := chnker.Next(nil);

		if (err != nil) {
			t.Fatalf("Error returned with chunk %d: %v", i, err);
		}

		if (c.Start != pos) {
			t.Fatalf("Start for chunk %d does not match: expected %d, got %d",
				i, pos, c.Start);
		}

		if (c.Length != chunk.Length) {
			t.Fatalf("Length for chunk %d does not match: expected %d, got %d",
				i, chunk.Length, c.Length);
		}

		if (c.Cut != chunk.CutFP) {
			t.Fatalf("Cut fingerprint for chunk %d/%d does not match: expected %016x, got %016x",
				i, len(chunks)-1, chunk.CutFP, c.Cut);
		}

		if (checkDigest) {
			auto digest = hashData(c.Data);
			if (!bytes.Equal(chunk.Digest, digest)) {
				t.Fatalf("Digest fingerprint for chunk %d/%d does not match: expected %02x, got %02x",
					i, len(chunks)-1, chunk.Digest, digest);
			}
		}

		pos += c.Length;
		chunks ~= c;
	}

	_, err := chnker.Next(nil);
	if (err != io.EOF) {
		t.Fatal("Wrong error returned after last chunk");
	}

	if (len(chunks) != len(testChunks)) {
		t.Fatal("Amounts of test and resulting chunks do not match");
	}

	return chunks;
}

private ubyte[] getRandom(long seed, int count) {
	auto buf = new ubyte[count];

	auto rnd = rand.New(rand.NewSource(seed));
	for (auto i = 0; i < count; i += 4) {
		auto r = rnd.Uint32();
		buf[i] = ubyte(r);
		buf[i+1] = ubyte(r >> 8);
		buf[i+2] = ubyte(r >> 16);
		buf[i+3] = ubyte(r >> 24);
	}

	return buf;
}

private ubyte[] hashData(ubyte[] d) {
	auto h = sha256.New();
	h.Write(d);
	return h.Sum(nil);
}

func TestChunker(testing*.T t) {
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bytes.NewReader(buf), testPol);
	testWithData(t, ch, chunks1, true);

	// setup nullbyte data source
	buf = bytes.Repeat(ubyte[]{0}, len(chunks2)MinSize*);
	ch = New(bytes.NewReader(buf), testPol);

	testWithData(t, ch, chunks2, true);
}

func TestChunkerWithCustomAverageBits(testing*.T t) {
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bytes.NewReader(buf), testPol);

	// sligthly decrease averageBits to get more chunks
	ch.SetAverageBits(19);
	testWithData(t, ch, chunks3, true);
}

func TestChunkerReset(testing*.T t) {
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bytes.NewReader(buf), testPol);
	testWithData(t, ch, chunks1, true);

	ch.Reset(bytes.NewReader(buf), testPol);
	testWithData(t, ch, chunks1, true);
}

func TestChunkerWithRandomPolynomial(testing*.T t) {
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	// generate a new random polynomial
	auto start = time.Now();
	p, err := RandomPolynomial();
	if (err != nil) {
		t.Fatal(err);
	}
	t.Logf("generating random polynomial took %v", time.Since(start));

	start = time.Now();
	auto ch = New(bytes.NewReader(buf), p);
	t.Logf("creating chunker took %v", time.Since(start));

	// make sure that first chunk is different
	c, err := ch.Next(nil);
	if (err != nil) {
		t.Fatal(err.Error());
	}

	if (c.Cut == chunks1[0].CutFP) {
		t.Fatal("Cut point is the same");
	}

	if (c.Length == chunks1[0].Length) {
		t.Fatal("Length is the same");
	}

	if (bytes.Equal(hashData(c.Data), chunks1[0].Digest)) {
		t.Fatal("Digest is the same");
	}
}

func TestChunkerWithoutHash(testing*.T t) {
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	auto ch = New(bytes.NewReader(buf), testPol);
	auto chunks = testWithData(t, ch, chunks1, false);

	// test reader
	foreach (i, c; chunks) {
		if (uint(len(c.Data)) != chunks1[i].Length) {
			t.Fatalf("reader returned wrong number of bytes: expected %d, got %d",
				chunks1[i].Length, len(c.Data));
		}

		if (!bytes.Equal(buf[c.Start .. c.Start+c.Length], c.Data)) {
			t.Fatalf("invalid data for chunk returned: expected %02x, got %02x",
				buf[c.Start .. c.Start+c.Length], c.Data);
		}
	}

	// setup nullbyte data source
	buf = bytes.Repeat(ubyte[]{0}, len(chunks2)MinSize*);
	ch = New(bytes.NewReader(buf), testPol);

	testWithData(t, ch, chunks2, false);
}

func benchmarkChunker(testing*.B b, bool checkDigest) {
	auto size = 32 * 1024 * 1024;
	auto rd = bytes.NewReader(getRandom(23, size));
	auto ch = New(rd, testPol);
	auto buf = new ubyte[MaxSize];

	b.ResetTimer();
	b.SetBytes(long(size));

	int chunks;
	for (auto i = 0; i < b.N; i++) {
		chunks = 0;

		_, err := rd.Seek(0, 0);
		if (err != nil) {
			b.Fatalf("Seek() return error %v", err);
		}

		ch.Reset(rd, testPol);

		auto cur = 0;
		while (true) {
			chunk, err := ch.Next(buf);

			if (err == io.EOF) {
				break;
			}

			if (err != nil) {
				b.Fatalf("Unexpected error occurred: %v", err);
			}

			if (chunk.Length != chunks1[cur].Length) {
				b.Errorf("wrong chunk length, want %d, got %d",
					chunks1[cur].Length, chunk.Length);
			}

			if (chunk.Cut != chunks1[cur].CutFP) {
				b.Errorf("wrong cut fingerprint, want 0x%x, got 0x%x",
					chunks1[cur].CutFP, chunk.Cut);
			}

			if (checkDigest) {
				auto h = hashData(chunk.Data);
				if (!bytes.Equal(h, chunks1[cur].Digest)) {
					b.Errorf("wrong digest, want %x, got %x",
						chunks1[cur].Digest, h);
				}
			}

			chunks++;
			cur++;
		}
	}

	b.Logf("%d chunks, average chunk size: %d bytes", chunks, size/chunks)
}

func BenchmarkChunkerWithSHA256(testing*.B b) {
	benchmarkChunker(b, true);
}

func BenchmarkChunker(testing*.B b) {
	benchmarkChunker(b, false);
}

func BenchmarkNewChunker(testing*.B b) {
	p, err := RandomPolynomial();
	if (err != nil) {
		b.Fatal(err);
	}

	b.ResetTimer();

	for (auto i = 0; i < b.N; i++) {
		New(bytes.NewBuffer(nil), p);
	}
}

// ----------------------------------------------------------- example_test.d
module chunker.example_test;

func ExampleChunker() {
	// generate 32MiB of deterministic pseudo-random data
	auto data = getRandom(23, 32*1024*1024);

	// create a chunker
	auto chunker = New(bytes.NewReader(data), Pol(0x3DA3358B4DC173));

	// reuse this buffer
	auto buf = new ubyte[8*1024*1024];

	for (auto i = 0; i < 5; i++) {
		chunk, err := chunker.Next(buf);
		if (err == io.EOF) {
			break;
		}

		if (err != nil) {
			panic(err);
		}

		fmt.Printf("%d %02x\n", chunk.Length, sha256.Sum256(chunk.Data));
	}

	// Output:
	// 2163460 4b94cb2cf293855ea43bf766731c74969b91aa6bf3c078719aabdd19860d590d
	// 643703 5727a63c0964f365ab8ed2ccf604912f2ea7be29759a2b53ede4d6841e397407
	// 1528956 a73759636a1e7a2758767791c69e81b69fb49236c6929e5d1b654e06e37674ba
	// 1955808 c955fb059409b25f07e5ae09defbbc2aadf117c97a3724e06ad4abd2787e6824
	// 2222372 6ba5e9f7e1b310722be3627716cf469be941f7f3e39a4c3bcefea492ec31ee56
}

// ----------------------------------------------------------- polynomials.d
module chunker.polynomials;

/// Pol is a polynomial from F_2[X].
type Pol ulong

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
		panic("multiplication would overflow ulong");
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
	if (ulong(x)&0xffffffff00000000 > 0) {
		x >>= 32;
		r |= 32;
	}

	if (ulong(x)&0xffff0000 > 0) {
		x >>= 16;
		r |= 16;
	}

	if (ulong(x)&0xff00 > 0) {
		x >>= 8;
		r |= 8;
	}

	if (ulong(x)&0xf0 > 0) {
		x >>= 4;
		r |= 4;
	}

	if (ulong(x)&0xc > 0) {
		x >>= 2;
		r |= 2;
	}

	if (ulong(x)&0x2 > 0) {
		r |= 1;
	}

	return r;
}

/// String returns the coefficients in hex.
public string String(/*this*/ Pol x) {
	return "0x" + strconv.FormatUint(ulong(x), 16);
}

/// Expand returns the string representation of the polynomial x.
public string Expand(/*this*/ Pol x) {
	if (x == 0) {
		return "0";
	}

	auto s = "";
	for (auto i = x.Deg(); i > 1; i--) {
		if (x&(1<<uint(i)) > 0) {
			s += fmt.Sprintf("+x^%d", i);
		}
	}

	if (x&2 > 0) {
		s += "+x";
	}

	if (x&1 > 0) {
		s += "+1";
	}

	return s[1 .. $];
}

/// DivMod returns x / d = q, and remainder r,
/// see https://en.wikipedia.org/wiki/Division_algorithm
public void DivMod(/*this*/ Pol x, Pol d) (Pol, Pol) {
	if (x == 0) {
		return 0, 0;
	}

	if (d == 0) {
		panic("division by zero");
	}

	auto D = d.Deg();
	auto diff = x.Deg() - D;
	if (diff < 0) {
		return 0, x;
	}

	Pol q;
	while (diff >= 0) {
		auto m = d << uint(diff);
		q |= (1 << uint(diff));
		x = x.Add(m);

		diff = x.Deg() - D;
	}

	return q, x;
}

/// Div returns the integer division result x / d.
public Pol Div(/*this*/ Pol x, Pol d) {
	q, _ := x.DivMod(d);
	return q;
}

/// Mod returns the remainder of x / d
public Pol Mod(/*this*/ Pol x, Pol d) {
	_, r := x.DivMod(d);
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
public (Pol, error) RandomPolynomial() {
	return DerivePolynomial(rand.Reader);
}

/// DerivePolynomial returns an irreducible polynomial of degree 53
/// (largest prime number below 64-8) by reading bytes from source.
/// There are (2^53-2/53) irreducible polynomials of degree 53 in
/// F_2[X], c.f. Michael O. Rabin (1981): "Fingerprinting by Random
/// Polynomials", page 4. If no polynomial could be found in one
/// million tries, an error is returned.
public (Pol, error) DerivePolynomial(io.Reader source) {
	for (auto i = 0; i < randPolMaxTries; i++) {
		Pol f;

		// choose polynomial at (pseudo)random
		auto err = binary.Read(source, binary.LittleEndian, &f);
		if (err != nil) {
			return 0, err;
		}

		// mask away bits above bit 53
		f &= Pol((1 << 54) - 1);

		// set highest and lowest bit so that the degree is 53 and the
		// polynomial is not trivially reducible
		f |= (1 << 53) | 1;

		// test if f is irreducible
		if (f.Irreducible()) {
			return f, nil;
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
		x, f = f, x;
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

/// MarshalJSON returns the JSON representation of the Pol.
public void MarshalJSON(/*this*/ Pol x) (ubyte[], error) {
	auto buf = strconv.AppendUint(ubyte[]{'"'}, ulong(x), 16);
	buf ~= '"';
	return buf, nil;
}

/// UnmarshalJSON parses a Pol from the JSON data.
public error UnmarshalJSON(/*this*/ Pol* x, ubyte[] data) {
	if (len(data) < 2) {
		throw new Exception("invalid string for polynomial");
	}
	n, err := strconv.ParseUint(string(data[1 .. len(data)-1]), 16, 64);
	if (err != nil) {
		return err;
	}
	x* = Pol(n);

	return nil;
}

// ----------------------------------------------------------- polynomials_test.d
module chunker.polynomials_test;

var polAddTests = struct[] {
	private Pol x, y;
	private Pol sum;
}{
	{23, 16, 23 ^ 16},
	{0x9a7e30d1e855e0a0, 0x670102a1f4bcd414, 0xfd7f32701ce934b4},
	{0x9a7e30d1e855e0a0, 0x9a7e30d1e855e0a0, 0},
}

func TestPolAdd(testing*.T t) {
	foreach (i, test; polAddTests) {
		if (test.sum != test.x.Add(test.y)) {
			t.Errorf("test %d failed: sum != x+y", i);
		}

		if (test.sum != test.y.Add(test.x)) {
			t.Errorf("test %d failed: sum != y+x", i);
		}
	}
}

private Pol parseBin(string s) {
	i, err := strconv.ParseUint(s, 2, 64);
	if (err != nil) {
		panic(err);
	}

	return Pol(i);
}

var polMulTests = struct[] {
	private Pol x, y;
	private Pol res;
}{
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
}

func TestPolMul(testing*.T t) {
	foreach (i, test; polMulTests) {
		auto m = test.x.Mul(test.y);
		if (test.res != m) {
			t.Errorf("TestPolMul failed for test %d: %v * %v: want %v, got %v",
				i, test.x, test.y, test.res, m);
		}
		m = test.y.Mul(test.x);
		if (test.res != test.y.Mul(test.x)) {
			t.Errorf("TestPolMul failed for %d: %v * %v: want %v, got %v",
				i, test.x, test.y, test.res, m);
		}
	}
}

func TestPolMulOverflow(testing*.T t) {
	defer func() {
		// try to recover overflow error
		auto err = recover();

		if (e, ok := err.(string); ok && e == "multiplication would overflow ulong") {
			return;
		}

		t.Logf("invalid error raised: %v", err);
		// re-raise error if not overflow
		panic(err);
	}();

	auto x = Pol(1 << 63);
	x.Mul(2);
	t.Fatal("overflow test did not panic");
}

var polDivTests = struct[] {
	private Pol x, y;
	private Pol res;
}{
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
}

func TestPolDiv(testing*.T t) {
	foreach (i, test; polDivTests) {
		auto m = test.x.Div(test.y);
		if (test.res != m) {
			t.Errorf("TestPolDiv failed for test %d: %v * %v: want %v, got %v",
				i, test.x, test.y, test.res, m);
		}
	}
}

func TestPolDeg(testing*.T t) {
	Pol x;
	if (x.Deg() != -1) {
		t.Errorf("deg(0) is not -1: %v", x.Deg());
	}

	x = 1;
	if (x.Deg() != 0) {
		t.Errorf("deg(1) is not 0: %v", x.Deg());
	}

	for (auto i = 0; i < 64; i++) {
		x = 1 << uint(i);
		if (x.Deg() != i) {
			t.Errorf("deg(1<<%d) is not %d: %v", i, i, x.Deg());
		}
	}
}

var polModTests = struct[] {
	private Pol x, y;
	private Pol res;
}{
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
}

func TestPolModt(testing*.T t) {
	foreach (i, test; polModTests) {
		auto res = test.x.Mod(test.y);
		if (test.res != res) {
			t.Errorf("test %d failed: want %v, got %v", i, test.res, res);
		}
	}
}

func BenchmarkPolDivMod(testing*.B t) {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < t.N; i++) {
		g.DivMod(f);
	}
}

func BenchmarkPolDiv(testing*.B t) {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < t.N; i++) {
		g.Div(f);
	}
}

func BenchmarkPolMod(testing*.B t) {
	auto f = Pol(0x2482734cacca49);
	auto g = Pol(0x3af4b284899);

	for (auto i = 0; i < t.N; i++) {
		g.Mod(f);
	}
}

func BenchmarkPolDeg(testing*.B t) {
	auto f = Pol(0x3af4b284899);
	auto d = f.Deg();
	if (d != 41) {
		t.Fatalf("BenchmalPolDeg: Wrong degree %d returned, expected %d",
			d, 41);
	}

	for (auto i = 0; i < t.N; i++) {
		f.Deg();
	}
}

func TestRandomPolynomial(testing*.T t) {
	_, err := RandomPolynomial();
	if (err != nil) {
		t.Fatal(err);
	}
}

func BenchmarkRandomPolynomial(testing*.B t) {
	for (auto i = 0; i < t.N; i++) {
		_, err := RandomPolynomial();
		if (err != nil) {
			t.Fatal(err);
		}
	}
}

func TestExpandPolynomial(testing*.T t) {
	auto pol = Pol(0x3DA3358B4DC173);
	auto s = pol.Expand();
	if (s != "x^53+x^52+x^51+x^50+x^48+x^47+x^45+x^41+x^40+x^37+x^36+x^34+x^32+x^31+x^27+x^25+x^24+x^22+x^19+x^18+x^16+x^15+x^14+x^8+x^6+x^5+x^4+x+1") {
		t.Fatal("wrong result");
	}
}

var polIrredTests = struct[] {
	private Pol f;
	private bool irred;
}{
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
}

func TestPolIrreducible(testing*.T t) {
	foreach (_, test; polIrredTests) {
		if (test.f.Irreducible() != test.irred) {
			t.Errorf("Irreducibility test for Polynomial %v failed: got %v, wanted %v",
				test.f, test.f.Irreducible(), test.irred);
		}
	}
}

func BenchmarkPolIrreducible(testing*.B b) {
	// find first irreducible polynomial
	Pol pol;
	foreach (_, test; polIrredTests) {
		if (test.irred) {
			pol = test.f;
			break;
		}
	}

	for (auto i = 0; i < b.N; i++) {
		if (!pol.Irreducible()) {
			b.Errorf("Irreducibility test for Polynomial %v failed", pol);
		}
	}
}

var polGCDTests = struct[] {
	private Pol f1;
	private Pol f2;
	private Pol gcd;
}{
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
}

func TestPolGCD(testing*.T t) {
	foreach (i, test; polGCDTests) {
		auto gcd = test.f1.GCD(test.f2);
		if (test.gcd != gcd) {
			t.Errorf("GCD test %d (%+v) failed: got %v, wanted %v",
				i, test, gcd, test.gcd);
		}

		gcd = test.f2.GCD(test.f1);
		if (test.gcd != gcd) {
			t.Errorf("GCD test %d (%+v) failed: got %v, wanted %v",
				i, test, gcd, test.gcd);
		}
	}
}

var polMulModTests = struct[] {
	private Pol f1;
	private Pol f2;
	private Pol g;
	private Pol mod;
}{
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
}

func TestPolMulMod(testing*.T t) {
	foreach (i, test; polMulModTests) {
		auto mod = test.f1.MulMod(test.f2, test.g);
		if (mod != test.mod) {
			t.Errorf("MulMod test %d (%+v) failed: got %v, wanted %v",
				i, test, mod, test.mod);
		}
	}
}
