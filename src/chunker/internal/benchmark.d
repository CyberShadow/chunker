/// A minimal benchmarking harness.
module chunker.internal.benchmark;

/**
 * To make a module benchmarkable, add:
 * ---
 * version (benchmarkYourModule) { import chunker.internal.benchmark; mixin BenchmarkThisModule; }
 * ---
 *
 * Then, declare benchmarks by declaring private functions with a name
 * starting with "benchmark", e.g. `void benchmarkFoo() { ... }`.
 *
 * To benchmark a benchmarkable module, run:
 * ---
 * dmd -version=benchmarkYourModule -run module.d [--N=<iterations>] [BenchmarkName ...]
 * ---
 */
mixin template BenchmarkThisModule()
{
	void main()
	{
		Benchmark.run();
	}

	struct Benchmark
	{
	static:
		immutable ulong iterations;
		private immutable string[] benchmarks;

		shared static this()
		{
			import core.runtime : Runtime;
			import std.getopt : getopt;
			ulong n = 1;
			auto args = Runtime.args;
			getopt(args,
				"N", &n,
			);
			Benchmark.iterations = n;

			auto benchmarks = args[1..$].idup;
			if (!benchmarks.length)
				benchmarks = ["*"];
			Benchmark.benchmarks = benchmarks;
		}

		private void run()
		{
			import std.stdio : stderr;
			import std.path : globMatch;

			alias mod = __traits(parent, main);
			foreach (name; __traits(allMembers, mod))
				static if (is(typeof(__traits(getMember, mod, name))))
				{
					alias member = __traits(getMember, mod, name);
					static if (__traits(isStaticFunction, member))
					{
						static if (name.length > 9 && name[0..9] == "benchmark")
						{
							bool found = false;
							foreach (mask; benchmarks)
								if (globMatch(name[9..$], mask))
								{
									found = true;
									break;
								}
							if (!found)
								continue;

							stderr.writefln("Running benchmark %s.%s (%d iterations):",
								__traits(identifier, mod), name[9..$], iterations);
							best = Duration.max;
							benchmarked = false;
							member();
							assert(benchmarked, "No block was benchmarked during this benchmark.");
							stderr.writefln("  -> %s", best);
						}
					}
				}
		}

		import std.datetime.stopwatch : StopWatch;
		import core.time : Duration;

		private Duration best;
		private bool benchmarked;

		void benchmark(void delegate() dg)
		{
			assert(!benchmarked, "A block was already benchmarked during this benchmark.");
			benchmarked = true;

			StopWatch sw;
			foreach (iteration; 0 .. iterations)
			{
				sw.start();
				dg();
				sw.stop();
				if (best > sw.peek())
					best = sw.peek();
				sw.reset();
			}
		}
	}
}
