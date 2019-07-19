## Variable sized arrays


## Memory Allocation

I first used `memset()` for memory allocation. `memset()` can be used easily for initiallizing arrays with -1, 0  or 1. It can also be modified for initialising with some other value. However, with floating point values, using `memset()` gives incorrect values. NOTE: `<string.h>` must be included to use `memset()`.

```c
// for 2D array
memset(array, 0, sizeof(array[0][0]) * Nx * Ny);
```

`memset()` in this case results in such an output. (Refer notes/memset.c for details)

```bash
conc[100][97] = 0.000000
conc[100][98] = -346351012320094344724873216.000000
conc[100][99] = 0.000000
conc[100][100] = -346352082231250619878866944.000000
```

## Using a function to initialise an array (using simple for loop)

I wrote a function that initialises an array to zero. My first guess was that since I was passing an empty array to the function and not returning it, the value of the array was getting updated within the function, but not getting returned and hence the error. [This answer](https://stackoverflow.com/a/27535098/9523246) may help understand such a situation.

However, this does not seem to be the case. As can be in the different versions of the `initialise_func_error.c` programs, the size of the array is returned incorrectly after function call. In `initialise_func_error_v3.c` the error seems to have been rectified but the reason why this works is not clear at the time of writing this note.

### Using gdb for debugging

To figure out where the segmentation fault lied in the erroneous versions, I tried using gdb. Frequently, it got stuck with the following output:

```bash
Reading symbols from ./a.out...
Reading symbols from /Users/jitinnair/OneDrive - Indian Institute of Science/TheGrind/Code/phasefieldC/spinodal/notes/a.out.dSYM/Contents/Resources/DWARF/a.out...
(gdb) r
Starting program: /Users/jitinnair/OneDrive - Indian Institute of Science/TheGrind/Code/phasefieldC/spinodal/notes/a.out
[New Thread 0x1003 of process 74135]
```

## Check Size of array

```c
//check array size
long arraySize = sizeof(array)/sizeof(array[0][0]);
printf("%ld\n", arraySize);
```

## Using malloc() and calloc()

`calloc()` has the advantage that it allocates a memory block and also initialises it to zero. An example of `malloc()` would be:

```c
// Syntax of malloc()
// ptr = (cast-type*) malloc(byte-size)
ptr = (int*) malloc(100 * sizeof(int));

// Syntax of calloc()
// ptr = (cast-type*)calloc(n, element-size);
ptr = (float*) calloc(25, sizeof(float));
```