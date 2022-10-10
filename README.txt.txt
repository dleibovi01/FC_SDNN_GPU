GPU Accelerated FC solver for hyperbolic conservation laws : GPU-Accelerated Fourier Continuation differentiation

CPU and GPU code
Author: Daniel Leibovici.

This is an implementation of a GPU-accelerated 3D Fourier Continuation-based spatial differentiation.
While the overall project (mentionned in the initial proposal and the CPU demo) remains to accelerate the whole
PDE solver, and thus the time stepping scheme in which this differentiation routine is called repeatedly, for the
final submission I only focused on this central and delicate routine. While not exactly trivial, the parallelization
of the rest of the solver on a GPU should not be harder in difficulty, but will require more time.

-----------------------------------------------------------------------------------------------------------------
Usage instructions:
I installed the Intel MKL libraries necessary to compile and run the code in my directory.
One should just run ./shockSolver

-----------------------------------------------------------------------------------------------------------------
Project description:
The program creates a 3D mesh of points on the cube [0, 1] x [0, 1] x [0, 1] and sets an initial condition defined
on these points. It then computes the derivative of this function with respect to the "y" coordinate, using the
3D Fourier-continuation method. In detail, a smooth-periodic extension is constructed in the y- direction by
the FcontGramBlend3D_y routine (defined in FC_3D.cpp), and then the derivative is computed by producing the Fourier
series of this extended function and differentiating in frequency space. Two routines are compared: diff3D_y and diff3D_y_GPU, 
both of which are defined in FC_3D.cpp. The Fourier series are produced by computing 1D FFTs by batch, then multiplying these
resulting FFTs by the spectral differentiation coefficients (done in the cudaImCmpProdScaleKernel in "VectorOperations.cu" 
for the GPU version), then recovering the derivative in physical space by computing inverse FFTs by batch again, and
restricting these extended solutions to the initial domain using the cudaRestrictionKernel.
-----------------------------------------------------------------------------------------------------------------
Results:
The test is performed on a 101 x 101 x 101 -point mesh. Both CPU and GPU versions give the same result, as shown
by the variable "discrepancy", indicating a discrepancy of order 1e-14.
The performance is measured in two ways: by including or excluding the time necessary to copy the data from host to device and back.
Excluding it makes sense, in the context of a larger solver, because the data should not be copied back and forth between
host and device at every differentiation step; only once at the first time step, and once at the final time step.
As can be seen in "screenshot results.png", a 17x speedup can be obtained if we exclude the copying of the data.

--------------------------------------------------------------------------------------------------------------------
Performance analysis:
A 17x speedup is satisfactory. However, it remains to be seen if the same speedups can be obtained with the other
routines, such as FcontGramBlend_y, as the objective is to speed-up the whole solver.







