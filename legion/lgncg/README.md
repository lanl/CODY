lgncg
=====

### NOTE TO USERS
At this moment, this code is under active development and is **NOT** stable nor
optimized.

### Description
A [Legion](http://legion.stanford.edu/) implementation of
[HPCG](https://software.sandia.gov/hpcg/).

### Building
```bash
export LG_RT_DIR=/path/to/legion/runtime
make
```

### Running
```bash
./drivers/hpcg/lgn-hpcg
```
or
```bash
./drivers/hpcg/lgn-hpcg -no-precond
```

### Help
./drivers/hpcg/lgn-hpcg -h

### TODO
- Optimize SPMV.
- Optimize SYGS.

**LA-CC 10-123**
