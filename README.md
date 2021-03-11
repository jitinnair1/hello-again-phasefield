# About
This repo has simple implementations of classical diffusion equation and spinodal decomposition in C solving which are usually the first steps of getting started with Phasefield modelling as well.

For examples implemented in MATLAB, [check out this repo](https://github.com/jitinnair1/hello-phasefield)

I have learned a lot looking at other people's code. I am posting my code here so that somebody starting out could use this as a template to build upon or to tinker around to understand what is happening.

# Build and Use Instructions
Use the following commands in you terminal:
```bash
git clone git@github.com:jitinnair1/hello-again-phasefield.git
cd hello-again-phasefield && mkdir build && cd build
cmake ..
make
```

Once the project is built, you will find the excutables under the `bin` folder of your copy of the code. To run the executable use `./` followed by the name of the executable like:

```bash
./diffusion_eqn_explicit
```

This will generate an output folder with `.vtk` files which can be opened using [ParaView](https://www.paraview.org/download/)

# a word of caution
This code was written at a time when I had just started learning C. So this code is not optimised for performance. Please create an issue if you encountered errors or unexpected behavior.
