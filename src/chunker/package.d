/**
Thin package implements Content Defined Chunking (CDC) based on a rolling
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
import chunker.rabin;

// Array support
import std.range.primitives : empty, front, popFront;

private enum size_t kiB = 1024;
private enum size_t miB = 1024 * kiB;

/// Aim to create chunks of 20 bits or about 1MiB on average.
public enum size_t chunkDefaultAverageBits = 20;
/// Default minimal size of a chunk.
public enum size_t chunkDefaultMinSize = 512 * kiB;
/// Default maximal size of a chunk.
public enum size_t chunkDefaultMaxSize = 8 * miB;


/// Splits content with Rabin Fingerprints.
struct Chunker(R)
{
	/// One content-dependent chunk of bytes whose end was cut when
	/// the Rabin Fingerprint had the value stored in `cut`.
	struct Chunk
	{
		/// Contents of the chunk.
		ubyte[] data;
		/// Value of the rolling hash when the chunk was cut.
		ulong cut;
	}

	private // internal state
	{
		/// Current hash internal state.
		RabinHash hash;

		/// Buffer used to receive and keep read data in.
		ubyte[] buf;

		/// Buffer used to copy chunk data to.
		ubyte[] cbuf;
		/// Number of bytes within the current chunk.
		size_t count;
	}

	private // configuration
	{
		/// Minimum and maximum chunk sizes, as configured.
		size_t minSize, maxSize;

		/// Hash mask used to decide chunk boundaries.
		/// By default `(1 << 20) - 1`, or configured from `averageBits`.
		ulong splitmask;

		/// Input data source.
		R source;
	}

	/// Probably use the `byCDChunk` convenience function instead of this constructor directly.
	/// Note that `Chunker` will `popFront` the input range during construction.
	public this(R source, Pol pol, uint averageBits, size_t minSize, size_t maxSize, ubyte[] cbuf)
	{
		this.cbuf = cbuf ? cbuf : new ubyte[maxSize];
		this.source = source;
		this.minSize = minSize;
		this.maxSize = maxSize;
		this.splitmask = (1L << averageBits) - 1;
		this.hash = RabinHash(pol);
		if (this.source.empty)
			empty = true;
		else
		{
			this.buf = this.source.front;
			popFront();
		}
	}

	Chunk front; /// Current chunk.
	bool empty; /// `true` when there are no more chunks (input range is empty).

	/// Updates `front` to contain the position and data of the next
	/// chunk. Note that `Chunker` reuses the chunk buffer, so the
	/// previous contents of `front.data` will be overwritten.
	public void popFront()
	{
		assert(!empty);

		hash.start();
		hash.slide(1);
		count = 0;
		auto minSize = this.minSize;
		auto maxSize = this.maxSize;
		// do not start a new chunk unless at least minSize bytes have been read
		auto pre = minSize - RabinHash.windowSize;
		auto buf = this.buf;
		scope(success) this.buf = buf;

		void copyBytes(size_t numBytes)
		{
			auto bytes = buf[0 .. numBytes];
			buf = buf[numBytes .. $];
			auto newLen = count + bytes.length;
			if (cbuf.length < newLen)
				cbuf.length = newLen;
			cbuf[count .. newLen] = bytes[];
			count += numBytes;
		}

		while (true)
		{
			if (buf.length == 0)
			{
				if (source.empty)
				{
					// Last chunk has been read.
					empty = true;
					return;
				}

				source.popFront();

				if (source.empty)
				{
					// There are no more bytes to buffer, so this was
					// the last chunk. Return current chunk, if any
					// bytes have been processed.
					if (count > 0)
						break;
					else
					{
						empty = true;
						return;
					}
				}

				buf = source.front;
			}

			// check if bytes have to be dismissed before starting a new chunk
			if (pre > 0)
			{
				if (pre > buf.length)
				{
					pre -= buf.length;
					copyBytes(buf.length);
					continue;
				}

				copyBytes(pre);
				pre = 0;
			}

			if (count < minSize)
			{
				auto warmUp = minSize - count;
				if (warmUp > buf.length)
					warmUp = buf.length;
				hash.put(buf[0 .. warmUp]);
				copyBytes(warmUp);
			}

			auto toWrite = buf.length;
			if (count + toWrite > maxSize)
				toWrite = maxSize - count;
			auto written = hash.putUntil(buf[0 .. toWrite], splitmask);
			copyBytes(written);
			if (buf.length) // stopped early due to maxSize or mask match
				break;
		}
		front = Chunk(cbuf[0 .. count], hash.peek());
	}
}

/// Constructs a new `Chunker` based on polynomial `pol` that
/// reads from `source`.
/// Params:
///   source = An input range of chunks of bytes, such as that
///     returned by `File.byChunk` or `ubyte[].chunks` or simply
///     `ubyte[].only`.
///   pol = An irreducible polynomial from `F_2[X]`. Use
///     `Pol.getRandom()` to generate a random one.
///   averageBits = Allows to control the frequency of chunk
///     discovery: the lower `averageBits`, the higher amount of
///     chunks will be identified.  The default value is 20 bits,
///     so chunks will be of 1MiB size on average.
///   minSize = Minimum size of emitted chunks. Chunk boundaries that
///     occur less than `minSize` bytes from the chunk start are
///     ignored.
///   maxSize = Maximum size of emitted chunks. If a chunk boundary is
///     not encountered after `maxSize` bytes from the chunk start, it
///     is forcibly split at that point.
///   cbuf = A buffer to store chunk data. When `null` (default), a
///     new buffer is allocated on construction of length `maxSize`.
/// Returns:
///   An instance of `Chunker`, an input range of `Chunker.Chunk`,
///   which contains the chunk data, and the fingerprint value when it
///   was cut.
Chunker!R byCDChunk(R)(
	R source,
	Pol pol,
	uint averageBits = chunkDefaultAverageBits,
	size_t minSize = chunkDefaultMinSize,
	size_t maxSize = chunkDefaultMaxSize,
	ubyte[] cbuf = null)
{
	return Chunker!R(source, pol, averageBits, minSize, maxSize, cbuf);
}

/// ditto
Chunker!R byCDChunk(R)(R source, Pol pol, ubyte[] cbuf)
{
	return Chunker!R(source, pol, chunkDefaultAverageBits, chunkDefaultMinSize, chunkDefaultMaxSize, cbuf);
}

// -----------------------------------------------------------------------------
version(unittest) version = test;
version(benchmarkChunker) version = test;
version(test):
private:

import chunker.internal.helpers;

enum chunkerBufSize = 512 * kiB;

template parseDigest(string s)
{
	import std.conv : hexString;
	enum ubyte[] parseDigest = cast(ubyte[])hexString!s;
}

struct TestChunk
{
	public size_t length;
	public ulong cutFP;
	public ubyte[] digest;
}

/// polynomial used for all the tests below
enum testPol = Pol(0x3DA3358B4DC173);

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
	{chunkDefaultMinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{chunkDefaultMinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{chunkDefaultMinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{chunkDefaultMinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
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

/// same as chunks1, but with a much bounded min/max size
TestChunk[] chunks4 =
[
	{1310720, 0x0008ec095c2a5e99, parseDigest!"5997b608c2a8b7d827bdfbd69dcffa86984c3a754aa4be92292f9f2ba0af4c22"},
	{852740, 0x000b98d4cdf00000, parseDigest!"8df8723fbcb4639c94cd831a3d8c8a6ff3ff746fe77769ee2244d1ef69bbe31e"},
	{1310720, 0x0008271396402fa4, parseDigest!"7c29efffbda477bbe56274573c3ba1b7e616eb76f3b9c92168a1f214aa557ba6"},
	{861939, 0x0015a25c2ef00000, parseDigest!"5080a1c4d3b538dce8c04c7738a043185d8c76e3251ae69e483b4b17d94dd7b9"},
	{1310720, 0x000fe23afea2e3be, parseDigest!"95507aa5fb821925e92ece579b5b50df0a6b311bd517937c4ae1a39e14cb79ea"},
	{1310720, 0x000e7cc6e6b093f0, parseDigest!"d2efa31fd050ca04232401b84a06fbde30bcc62f2108a5c0f52cca319952ab90"},
	{1310720, 0x00090b61ebb67fc7, parseDigest!"3a7e29b6a323b111ff5860953de004dfb0f69e5ccadb8831176cb0a6e3563262"},
	{1310720, 0x0008998f21ddfe08, parseDigest!"d25815c6d82e29d8d59fdde4bd85804526cd45b027106e9526d2bfdba6e8c7ae"},
	{1310720, 0x0010ebc94474a0db, parseDigest!"1ef10b00b9894102517e4a2d201877c9eba662dcbd22b2259f44c2c88436cbae"},
	{1310720, 0x00180e6341e54abb, parseDigest!"d7a9ebce0573e5a2f6d688a495044766b2323a7a8b73ad2a9f1bdead52429d83"},
	{1310720, 0x0010fa125ffa7553, parseDigest!"9803d879cf7cc8221eae87c2bf1ac9e6a4328929dfe32e5cb461663186402531"},
	{1310720, 0x000884f3046bc00f, parseDigest!"0d35b3a8949b673821f004eb845af80ae34f082dd4c20b7715c42aff21aeabad"},
	{1310720, 0x00154b752483d77e, parseDigest!"fb464beea91aee4eb91b7e5da58cc51f28a9ee74a61bc5ea40e703603e860fd8"},
	{1310720, 0x00068be82cf57932, parseDigest!"ab4ebe4015214c4f449782ded96011d3ac1be3daa682a1df2d1905104008f1b5"},
	{1310720, 0x000a1c2e12d6c196, parseDigest!"9198f8799fc64d2f6da6365a42f32c7b2d31ebfb23ccd957b705cfc24c69bc5e"},
	{1310720, 0x000f73c53cb2ae82, parseDigest!"2d05e2f054b3ce480e7ffa4490db250ed2126fe36b8d767fc51c6e2d3baff87b"},
	{915082, 0x000f55a1fea00000, parseDigest!"bcc88a08794bf3f194afa0c351836846eed95acfe7de41dc32d7f63c8fa02a4c"},
	{1310720, 0x0016d908ad8f8bbc, parseDigest!"214f5fdf4136b08499ae944d6507ae97ae0eee5e051e336182a18d0881251d96"},
	{1310720, 0x001740bda9456bd2, parseDigest!"1fcdad3287449b74a61ac626599dc9b0efb7bf01b599931f141f0ec340547c06"},
	{1310720, 0x000c1606d9a0d99f, parseDigest!"d88982bf999af432cfcc94bace5000308cbb1cc6cd2126286d6304ebd83c4765"},
	{837117, 0x00030ce2d9400000, parseDigest!"bf6c025df617eac82f643026d1d0295159a63a4e8dea41ed0341bd0469cbd7d3"},
	{1310720, 0x0008fb9882a96b1f, parseDigest!"b5a2a7d0dfaa4d66e6c693afe77037cd42a21deac68f21fdafdfdeeea881e6a1"},
	{1151715, 0x000968473f900000, parseDigest!"31882d126190910a0260abe5d10cb31b180e3e2d0d2c3e45981a78e1675738a5"},
	{1310720, 0x0015cadcf5b2d7ac, parseDigest!"7d878ceb8e6ef8a3adc8a892a6c32d8afc93b9ea6b8a55d2271230847fa4b925"},
	{1142792, 0x001e197c92600000, parseDigest!"96095ef256bb29934c024dbd0e3eab9c9dee2aa3308a7f4d77f20ae712ea57f1"},
	{1310720, 0x0005408e6716fd12, parseDigest!"2825054e95b14ceb9f6c1947e5f5da0b34ce898850e05390adfcd976a271c107"},
	{1310720, 0x000fb439ba516157, parseDigest!"eb4d3e5686f4cc810faa8bb356dbdd3242cb9acfd87fa7c2229c44f3f33aa714"},
	{267927, 0x0000000000000001, parseDigest!"19f1aa1f5c3b49452675de394497fe5943cd00ee5a61291c6c0ca92b01ca2312"},
];

import std.format : format;
import std.digest.sha : sha256Of;

Chunker!R.Chunk[] testWithData(R)(ref Chunker!R chnker, TestChunk[] testChunks, bool checkDigest, bool copyData = false)
{
	Chunker!R.Chunk[] chunks;

	foreach (i, chunk; testChunks)
	{
		auto c = chnker.front;

		assert(c.data.length == chunk.length,
			format!"Length for chunk %d does not match: expected %d, got %d"
			(i, chunk.length, c.data.length));

		assert(c.cut == chunk.cutFP,
			format!"Cut fingerprint for chunk %d/%d does not match: expected %016x, got %016x"
			(i, testChunks.length, chunk.cutFP, c.cut));

		if (checkDigest)
		{
			auto digest = sha256Of(c.data);
			assert(chunk.digest == digest,
				format!"Digest fingerprint for chunk %d/%d does not match: expected %(%02x%), got %(%02x%)"
				(i, testChunks.length, chunk.digest, digest));
		}

		if (copyData)
			c.data = c.data.dup;
		chunks ~= c;
		chnker.popFront();
	}

	if (!chnker.empty)
		assert(false, "Wrong error returned after last chunk");

	assert(chunks.length == testChunks.length,
		"Amounts of test and resulting chunks do not match");

	return chunks;
}

import std.array : replicate;
import std.range : chunks;

@(`Chunker`) unittest
{
	// setup data source
	auto data = getRandom(23, 32*1024*1024);
	auto ch = data.chunks(chunkerBufSize).byCDChunk(testPol);
	testWithData(ch, chunks1, true);

	// setup nullbyte data source
	data = replicate([ubyte(0)], chunks2.length*chunkDefaultMinSize);
	ch = data.chunks(chunkerBufSize).byCDChunk(testPol);

	testWithData(ch, chunks2, true);
}

@(`ChunkerWithArrayInput`) unittest
{
	// Test with all data in one array/chunk
	auto data = getRandom(23, 32*1024*1024);
	auto ch = [data].byCDChunk(testPol);
	testWithData(ch, chunks1, true);

	data = replicate([ubyte(0)], chunks2.length*chunkDefaultMinSize);
	ch = [data].byCDChunk(testPol);
	testWithData(ch, chunks2, true);
}

@(`ChunkerWithFileInput`) unittest
{
	// Test with File.byChunk
	import std.stdio : File;
	import std.file : remove;
	auto data = getRandom(23, 32*1024*1024);
	File("test.bin", "wb").rawWrite(data);
	scope(exit) remove("test.bin");
	auto ch = File("test.bin", "rb").byChunk(chunkerBufSize).byCDChunk(testPol);
	testWithData(ch, chunks1, true);

	data = replicate([ubyte(0)], chunks2.length*chunkDefaultMinSize);
	File("test.bin", "wb").rawWrite(data);
	ch = File("test.bin", "rb").byChunk(chunkerBufSize).byCDChunk(testPol);
	testWithData(ch, chunks2, true);
}

@(`ChunkerWithCustomAverageBits`) unittest
{
	auto data = getRandom(23, 32*1024*1024);

	// sligthly decrease averageBits to get more chunks
	auto ch = data.chunks(chunkerBufSize).byCDChunk(testPol, 19);

	testWithData(ch, chunks3, true);
}

@(`ChunkerWithCustomMinMax`) unittest
{
	auto data = getRandom(23, 32*1024*1024);
	auto ch = data.chunks(chunkerBufSize).byCDChunk(testPol,
		20,
		(1 << 20) - (1 << 18),
		(1 << 20) + (1 << 18));

	testWithData(ch, chunks4, true);
}

// Ensure the min/max bounds are strictly followed.
// (Check for off-by-one errors.)
@(`ChunkerMinMaxBounds`) unittest
{
	auto data = getRandom(23, 64*1024);
	auto ch = data.chunks(chunkerBufSize).byCDChunk(testPol,
		7,
		RabinHash.windowSize * 2 - 2,
		RabinHash.windowSize * 2 + 2);

	size_t[] sizes;
	foreach (chunk; ch)
		if (chunk.data.length != RabinHash.windowSize * 2 + 2)
			sizes ~= chunk.data.length;

	assert(sizes == [
			126, 129, 126, 127, 128, 128, 126, 126, 127, 129,
			129, 128, 127, 128, 126, 129, 126, 129, 128, 127, 67]);
}

import std.stdio : stderr;

@(`ChunkerWithRandomPolynomial`) unittest
{
	// setup data source
	auto data = getRandom(23, 32*1024*1024);

	// generate a new random polynomial
	import std.datetime.systime : Clock, SysTime;
	auto start = Clock.currTime();
	auto p = Pol.getRandom();
	stderr.writefln!"generating random polynomial took %s"(Clock.currTime() - start);

	start = Clock.currTime();
	auto ch = data.chunks(chunkerBufSize).byCDChunk(p);
	stderr.writefln!"creating chunker took %s"(Clock.currTime() - start);

	// make sure that first chunk is different
	auto c = ch.front;

	assert(c.cut != chunks1[0].cutFP,
		"Cut point is the same");

	assert(c.data.length != chunks1[0].length,
		"Length is the same");

	assert(sha256Of(c.data) != chunks1[0].digest,
		"Digest is the same");
}

@(`ChunkerWithoutHash`) unittest
{
	// setup data source
	auto data = getRandom(23, 32*1024*1024);

	auto ch = data.chunks(chunkerBufSize).byCDChunk(testPol);
	auto chunks = testWithData(ch, chunks1, false, true);

	// test reader
	size_t pos = 0;
	foreach (i, c; chunks)
	{
		assert(c.data.length == chunks1[i].length,
			format!"reader returned wrong number of bytes: expected %d, got %d"
			(chunks1[i].length, c.data.length));

		assert(data[pos .. pos+c.data.length] == c.data,
			format!"invalid data for chunk returned: expected %(%02x%), got %(%02x%)"
			(data[pos .. pos+c.data.length], c.data));
		pos += c.data.length;
	}

	// setup nullbyte data source
	data = replicate([ubyte(0)], chunks2.length*chunkDefaultMinSize);
	ch = data.chunks(chunkerBufSize).byCDChunk(testPol);

	testWithData(ch, chunks2, false);
}

version (benchmarkChunker)
{
	import chunker.internal.benchmark;
	mixin BenchmarkThisModule;

	void _benchmarkChunker(bool checkDigest)
	{
		auto size = 32 * 1024 * 1024;
		auto buf = new ubyte[chunkDefaultMaxSize];
		auto data = [getRandom(23, size)];

		// b.SetBytes(long(size));

		int chunks;
		Benchmark.benchmark({
			chunks = 0;

			auto ch = data.byCDChunk(testPol, buf);

			auto cur = 0;
			while (!ch.empty)
			{
				auto chunk = ch.front;

				assert(chunk.data.length == chunks1[cur].length,
					format!"wrong chunk length, want %d, got %d"
					(chunks1[cur].length, chunk.data.length));

				assert(chunk.cut == chunks1[cur].cutFP,
					format!"wrong cut fingerprint, want 0x%x, got 0x%x"
					(chunks1[cur].cutFP, chunk.cut));

				if (checkDigest)
				{
					auto h = sha256Of(chunk.data);
					assert(h == chunks1[cur].digest,
						format!"wrong digest, want %(%02x%), got %(%02x%)"
						(chunks1[cur].digest, h));
				}

				chunks++;
				cur++;
				ch.popFront();
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
		ubyte[][] buf;

		Benchmark.benchmark({
			buf.byCDChunk(p);
		});
	}
}
