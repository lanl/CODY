CODY: Continuum Dynamics Evaluation and Test Suite
===========

Description
-----------

CODY is a development framework and suite of small applications, or
"mini-apps", characteristic of continuum dynamics applications that
will be used for research in new programming models, software
environments, and the evaluation of new computer architectures.

On what platforms does it run?

CODY is initially a suite of applications written in C and C++ with
extensions -- including CUDA and OpenCL -- that enable the use of
emerging computer architectures.  Our goals are operating system,
computer architecture, and programming language agnostic, however, so
we expect evolution over time.

For whom is it designed?

CODY is designed for computer science researchers to assess and
explore new directions in computer architecture and programming
models.

Why was it developed?

CODY was developed at Los Alamos National Laboratory to provide a set
of small, manageable applications that are representative of some
scientific workloads at LANL, in a form that can be readily shared
with research partners at other organizations.

Installation
------------

If necessary, set the `GOROOT` environment variable to the directory
containing `src/Make.inc` and `src/Make.pkg` and the `PAPI_INCDIR`
environment variable to the directory containing `papi.h`.  Also,
ensure that the directory containing `libpapi.so` is listed in your
`LD_LIBRARY_PATH`.

Afterwards, you can follow the usual Go package installation
procedure:

<pre>
    git clone http://github.com/losalamos/go-papi $GOROOT/src/pkg/github.com/losalamos/go-papi
    cd $GOROOT/src/pkg/github.com/losalamos/go-papi
    gomake
    gotest
    gomake install
</pre>

It is then safe to do a `gomake clean`.

At the time of this writing,
[`goinstall`](http://golang.org/cmd/goinstall/) is unable to install
packages such as go-papi that require
[`cgo`](http://golang.org/cmd/cgo/).  If this is ever fixed, the
preceding steps can be simplified into

<pre>
    goinstall github.com/losalamos/go-papi
</pre>


Documentation
-------------

Pre-built documentation for the core part of the go-papi API is
available online at
<http://gopkgdoc.appspot.com/pkg/github.com/losalamos/go-papi>,
courtesy of [GoPkgDoc](http://gopkgdoc.appspot.com/).  Unfortunately,
the online documentation omits descriptions of all constants,
variables, etc. that are generated during the build process,
specifically the list of PAPI events (`papi-event.go`), event
modifiers (`papi-emod.go`), and error values (`papi-errno.go`).

Once you install go-papi, you can view the complete go-papi API with
[`godoc`](http://golang.org/cmd/godoc/), for example by running

<pre>
    godoc -http=:6060
</pre>

to start a local Web server then viewing the documentation at
<http://localhost:6060/pkg/github.com/losalamos/go-papi/> in your
favorite browser.

For code examples, take a look at the `*_test.go` files in the go-papi
source distribution.  `papi_hl_test.go` utilizes PAPI's high-level
API; `papi_ll_test.go` utilizes PAPI's low-level API; and
`papi_test.go` utilizes a few miscellaneous functions.


License
-------

BSD-ish with a "modifications must be indicated" clause.  See
<http://github.com/losalamos/go-papi/blob/master/LICENSE> for the full
text.


Author
------

Scott Pakin, <pakin@lanl.gov>


