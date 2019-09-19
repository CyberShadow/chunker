[![DUB Package](https://img.shields.io/dub/v/chunker.svg)](https://code.dlang.org/packages/chunker)
[![Build Status](https://travis-ci.org/CyberShadow/chunker.svg?branch=master)](https://travis-ci.org/CyberShadow/chunker)

The package `chunker` implements content-defined-chunking (CDC) based on a
rolling Rabin Hash. This library is a D port of [the Go
package](https://github.com/restic/chunker), which is part of the [restic
backup program](https://github.com/restic/restic).

An introduction to Content Defined Chunking can be found in the restic blog
post [Foundation - Introducing Content Defined Chunking (CDC)](https://restic.github.io/blog/2015-09-12/restic-foundation1-cdc).

You can find the API documentation at
https://chunker.dpldocs.info/
