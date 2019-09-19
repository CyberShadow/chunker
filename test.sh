#!/bin/bash
set -eEuo pipefail

dmd -Isrc -i -g -unittest -main               -run src/chunker/package
dmd -Isrc -i -g                               -run src/chunker/example
dmd -Isrc -i -g -version=benchmarkChunker     -run src/chunker/package
dmd -Isrc -i -g -version=benchmarkPolynomials -run src/chunker/polynomials

dub build
dub test
dub run chunker:example
