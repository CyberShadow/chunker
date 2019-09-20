/// Helpers used in tests and examples.
module chunker.internal.helpers;

package(chunker):

private import chunker.internal.gorng;

ubyte[] getRandom(int seed, int count)
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

// Temporary D shim
import std.stdio : File;
File bufFile(ubyte[] buf)
{
	File("temp.bin", "wb").rawWrite(buf);
	return File("temp.bin", "rb");
}
