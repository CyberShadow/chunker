/// A minimal benchmarking harness.
module chunker.benchmark;

/**
 * To make a module benchmarkable, add:
 * version (benchmarkYourModule) { import chunker.benchmark; mixin BenchmarkThisModule; }
 *
 * Then, declare benchmarks by declaring private functions with a name
 * starting with "benchmark", e.g. `void benchmarkFoo() { ... }`.
 *
 * To benchmark a benchmarkable module, run:
 * dmd -version=benchmarkYourModule -run module.d [--N=<iterations>] [BenchmarkName ...]
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
		immutable ulong N;
		private immutable string[] benchmarks;

		shared static this()
		{
			import core.runtime : Runtime;
			import std.getopt : getopt;
			ulong N = 1;
			auto args = Runtime.args;
			getopt(args,
				"N", &N,
			);
			Benchmark.N = N;

			auto benchmarks = args[1..$].idup;
			if (!benchmarks.length)
				benchmarks = ["*"];
			Benchmark.benchmarks = benchmarks;
		}

		import std.datetime.stopwatch : StopWatch;
		private StopWatch sw;

		void resetTimer()
		{
			sw.reset();
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
								__traits(identifier, mod), name[9..$], N);
							sw.reset();
							sw.start();
							member();
							auto time = sw.peek();
							stderr.writefln("  -> %s", time);
						}
					}
				}
		}
	}
}
