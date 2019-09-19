module chunker.chunker_test;

import std.stdio : File;

import chunker;
import chunker.polynomials;

private template parseDigest(string s)
{
	import std.conv : hexString;
	enum ubyte[] parseDigest = cast(ubyte[])hexString!s;
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
chunk[] chunks1 =
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

/// test if nullbytes are correctly split, even if length is a multiple of MinSize.
chunk[] chunks2 =
[
	{MinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{MinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{MinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
	{MinSize, 0, parseDigest!"07854d2fef297a06ba81685e660c332de36d5d18d546927d30daad6d7fda1541"},
];

/// the same as chunks1, but avg chunksize is 1<<19
chunk[] chunks3 =
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

version(unittest) import std.format;

version(unittest) private Chunk[] testWithData(Chunker* chnker, chunk[] testChunks, bool checkDigest)
{
	Chunk[] chunks;

	auto pos = uint(0);
	foreach (i, chunk; testChunks)
	{
		auto c = chnker.Next(null);

		if (c.Start != pos)
		{
			assert(false, format!"Start for chunk %d does not match: expected %d, got %d"
				(i, pos, c.Start));
		}

		if (c.Length != chunk.Length)
		{
			assert(false, format!"Length for chunk %d does not match: expected %d, got %d"
				(i, chunk.Length, c.Length));
		}

		if (c.Cut != chunk.CutFP)
		{
			assert(false, format!"Cut fingerprint for chunk %d/%d does not match: expected %016x, got %016x"
				(i, testChunks.length-1, chunk.CutFP, c.Cut));
		}

		if (checkDigest)
		{
			auto digest = hashData(c.Data);
			if (!(chunk.Digest == digest))
			{
				assert(false, format!"Digest fingerprint for chunk %d/%d does not match: expected %(%02x%), got %(%02x%)"
					(i, testChunks.length-1, chunk.Digest, digest));
			}
		}

		pos += c.Length;
		chunks ~= c;
	}

	auto c = chnker.Next(null);
	if (c !is Chunk.init)
	{
		assert(false, "Wrong error returned after last chunk");
	}

	if (chunks.length != testChunks.length)
	{
		assert(false, "Amounts of test and resulting chunks do not match");
	}

	return chunks;
}

private import chunker.gorng;

package ubyte[] getRandom(int seed, int count)
{
	import std.random : Random, uniform;
	auto buf = new ubyte[count];

	rngSource rnd;
	Seed(&rnd, seed);
	for (auto i = 0; i < count; i += 4)
	{
		auto r = Int63(&rnd) >> 31;
		buf[i] = cast(ubyte)(r);
		buf[i+1] = cast(ubyte)(r >> 8);
		buf[i+2] = cast(ubyte)(r >> 16);
		buf[i+3] = cast(ubyte)(r >> 24);
	}

	return buf;
}

version(unittest) private ubyte[32] hashData(ubyte[] d)
{
	import std.digest.sha;
	return sha256Of(d);
}

version(unittest) import std.array : replicate;

// Temporary D shim
package File bufFile(ubyte[] buf)
{
	File("temp.bin", "wb").rawWrite(buf);
	return File("temp.bin", "rb");
}

@(`Chunker`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bufFile(buf), testPol);
	testWithData(ch, chunks1, true);

	// setup nullbyte data source
	buf = replicate([ubyte(0)], chunks2.length*MinSize);
	ch = New(bufFile(buf), testPol);

	testWithData(ch, chunks2, true);
}

@(`ChunkerWithCustomAverageBits`) unittest
{
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bufFile(buf), testPol);

	// sligthly decrease averageBits to get more chunks
	ch.SetAverageBits(19);
	testWithData(ch, chunks3, true);
}

@(`ChunkerReset`) unittest
{
	auto buf = getRandom(23, 32*1024*1024);
	auto ch = New(bufFile(buf), testPol);
	testWithData(ch, chunks1, true);

	ch.Reset(bufFile(buf), testPol);
	testWithData(ch, chunks1, true);
}

version(unittest) import std.stdio : stderr;

@(`ChunkerWithRandomPolynomial`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	// generate a new random polynomial
	import std.datetime.systime : Clock, SysTime;
	auto start = Clock.currTime();
	auto p = RandomPolynomial();
	stderr.writefln!"generating random polynomial took %s"(Clock.currTime() - start);

	start = Clock.currTime();
	auto ch = New(bufFile(buf), p);
	stderr.writefln!"creating chunker took %s"(Clock.currTime() - start);

	// make sure that first chunk is different
	auto c = ch.Next(null);

	if (c.Cut == chunks1[0].CutFP)
	{
		assert(false, "Cut point is the same");
	}

	if (c.Length == chunks1[0].Length)
	{
		assert(false, "Length is the same");
	}

	if (hashData(c.Data) == chunks1[0].Digest)
	{
		assert(false, "Digest is the same");
	}
}

@(`ChunkerWithoutHash`) unittest
{
	// setup data source
	auto buf = getRandom(23, 32*1024*1024);

	auto ch = New(bufFile(buf), testPol);
	auto chunks = testWithData(ch, chunks1, false);

	// test reader
	foreach (i, c; chunks)
	{
		if (c.Data.length != chunks1[i].Length)
		{
			assert(false, format!"reader returned wrong number of bytes: expected %d, got %d"
				(chunks1[i].Length, c.Data.length));
		}

		if (!(buf[c.Start .. c.Start+c.Length] == c.Data))
		{
			assert(false, format!"invalid data for chunk returned: expected %(%02x%), got %(%02x%)"
				(buf[c.Start .. c.Start+c.Length], c.Data));
		}
	}

	// setup nullbyte data source
	buf = replicate([ubyte(0)], chunks2.length*MinSize);
	ch = New(bufFile(buf), testPol);

	testWithData(ch, chunks2, false);
}

version = benchmark;
version (benchmark) enum N = 1;

version (benchmark) void benchmarkChunker(bool checkDigest)
{
	auto size = 32 * 1024 * 1024;
	auto rd = bufFile(getRandom(23, size));
	auto ch = New(rd, testPol);
	auto buf = new ubyte[MaxSize];

	// b.ResetTimer();
	// b.SetBytes(long(size));

	int chunks;
	for (auto i = 0; i < N; i++)
	{
		chunks = 0;

		rd.seek(0);

		ch.Reset(rd, testPol);

		auto cur = 0;
		while (true)
		{
			auto chunk = ch.Next(buf);

			if (chunk is Chunk.init)
			{
				break;
			}

			if (chunk.Length != chunks1[cur].Length)
			{
				assert(false, format!"wrong chunk length, want %d, got %d"
				(
					chunks1[cur].Length, chunk.Length));
			}

			if (chunk.Cut != chunks1[cur].CutFP)
			{
				assert(false, format!"wrong cut fingerprint, want 0x%x, got 0x%x"
				(
					chunks1[cur].CutFP, chunk.Cut));
			}

			if (checkDigest)
			{
				auto h = hashData(chunk.Data);
				if (!(h == chunks1[cur].Digest))
				{
					assert(false, format!"wrong digest, want %(%02x%), got %(%02x%)"
					(
						chunks1[cur].Digest, h));
				}
			}

			chunks++;
			cur++;
		}
	}

	stderr.writefln!"%d chunks, average chunk size: %d bytes"(chunks, size/chunks);
}

version (benchmark) @(`BenchmarkChunkerWithSHA256`) unittest
{
	benchmarkChunker(true);
}

version (benchmark) @(`BenchmarkChunker`) unittest
{
	benchmarkChunker(false);
}

version (benchmark) @(`BenchmarkNewChunker`) unittest
{
	auto p = RandomPolynomial();

	// b.ResetTimer();

	for (auto i = 0; i < N; i++)
	{
		New(bufFile(null), p);
	}
}

