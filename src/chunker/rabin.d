module chunker.rabin;

import chunker.polynomials;

/// Calculates a streaming Rabin Fingerprint.
struct RabinHash
{
	// cache precomputed tables, these are read-only anyway
	private struct Cache
	{
	static:
		Tables[Pol] entries;
		Object mutex;
	}
	static Cache cache;

	static this()
	{
		cache.mutex = new Object;
	}

	private struct Tables
	{
		private Pol[256] out_;
		private Pol[256] mod;
	}

	/// Precomputed tables used for hashing.
	private Tables tables;
	/// Whether `tables` have already been computed.
	private bool tablesInitialized;

	/// `fillTables` calculates `out_table` and `mod_table` for
	/// optimization. This implementation uses a cache in the global
	/// variable `cache`.
	private void fillTables(Pol pol)
	{
		tablesInitialized = true;

		// test if the tables are cached for this polynomial
		synchronized(cache.mutex)
		{
			if (auto t = pol in cache.entries)
			{
				tables = *t;
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

				h = appendByte(h, cast(ubyte)b, pol);
				foreach (i; 0 .. RabinHash.windowSize-1)
					h = appendByte(h, 0, pol);
				tables.out_[b] = h;
			}

			// calculate table for reduction mod Polynomial
			auto k = pol.deg();
			foreach (b; 0 .. 256)
			{
				// mod_table[b] = A | B, where A = (b(x) * x^k mod pol) and  B = b(x) * x^k
				//
				// The 8 bits above deg(Polynomial) determine what happens next and so
				// these bits are used as a lookup to this table. The value is split in
				// two parts: Part A contains the result of the modulus operation, part
				// B is used to cancel out the 8 top bits so that one XOR operation is
				// enough to reduce modulo Polynomial
				tables.mod[b] = Pol(
					(Pol(Pol.Base(b) << uint(k)) % pol).value |
					    (Pol.Base(b) << uint(k)));
			}

			cache.entries[pol] = tables;
		}
	}

	/// Bits to shift the digest when updating the hash.
	/// Calculated from the polynomial's degree.
	private uint polShift;

	/// Initializes this `RabinHash` with the given polynomial.
	this(Pol pol)
	{
		polShift = uint(pol.deg() - 8);
		fillTables(pol);
	}

	// ---------------------------------------------------------------------

	/// Base type for the calculated digest.
	alias Digest = ulong;

	/// Size of the sliding window.
	public enum size_t windowSize = 64;

	/// Current contents of the sliding window.
	private ubyte[windowSize] window;

	/// Contains frequently-accessed scalar values.
	/// See Writer for details.
	private struct Scalars
	{
		/// Sliding window cursor.
		private size_t wpos;

		/// Current hash digest value.
		private Digest digest;
	}
	/// ditto
	private Scalars scalars;

	/// Reset internal state.
	void start()
	{
		foreach (i; 0 .. windowSize)
			window[i] = 0;
		scalars.digest = 0;
		scalars.wpos = 0;
	}

	/// Return current digest.
	@property Digest peek() const { return scalars.digest; }

	/// Return current digest, then reset.
	Digest finish()
	{
		auto result = scalars.digest;
		start();
		return result;
	}

	// ---------------------------------------------------------------------

	/// Type for fast writing of successive bytes.
	/// Must contain only scalar values, so that optimizing compilers
	/// may break it up and place the fields in registers.
	struct Writer
	{
		/// Contained scalar values.
		private Scalars scalars;
		/// Pointer to the rest of the object.
		private RabinHash* hash;

		/// Return current digest.
		public @property Digest peek() const { return scalars.digest; }

		/// Update the hash with one byte.
		public void slide(ubyte b)
		{
			slideImpl(hash.window, scalars.wpos, scalars.digest, hash.tables.out_, hash.tables.mod, hash.polShift, b);
		}
	}

	/// Return a fast writer object.
	/// Must be committed at the end using `commit`.
	Writer getWriter()
	{
		assert(tablesInitialized, "tables for polynomial computation not initialized");
		return Writer(scalars, &this);
	}

	/// Commit a `Writer`.
	void commit(ref Writer writer)
	{
		assert(writer.hash == &this);
		this.scalars = writer.scalars;
	}

	/// Update the hash with one byte.
	public void slide(ubyte b)
	{
		slideImpl(window, scalars.wpos, scalars.digest, tables.out_, tables.mod, polShift, b);
	}

	/// Implementation for `slide`.
	private static void slideImpl(ref ubyte[windowSize] window, ref size_t wpos, ref Digest digest, ref Pol[256] tabout, in ref Pol[256] tabmod, uint polShift, ubyte b)
	{
		auto out_ = window[wpos];
		window[wpos] = b;
		digest ^= Digest(tabout[out_].value);
		wpos++;
		if (wpos >= windowSize)
			wpos = 0;

		digest = updateDigest(digest, polShift, tabmod, b);
	}

	/// ditto
	private static Digest /*newDigest*/ updateDigest(Digest digest, uint polShift, in ref Pol[256] tabmod, ubyte b)
	{
		auto index = cast(ubyte)(digest >> polShift);
		digest <<= 8;
		digest |= Digest(b);

		digest ^= Digest(tabmod[index].value);
		return digest;
	}
}

private Pol appendByte(Pol hash, ubyte b, Pol pol)
{
	hash.value <<= 8;
	hash.value |= Pol.Base(b);

	return hash % pol;
}
