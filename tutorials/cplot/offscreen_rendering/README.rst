This tutorial demonstrates how to use offscreen rendering and dump an image or a mpeg movie.

      Rotating cube by C. Benoit.

CPlot can render in an image file or in a mpeg movie file.
If you are using CPlot on a cluster (with no graphic card), it must be installed with mesa. Then, when using CPlot.display, offscreen is triggered with offscreen=1.
If you are using CPlot on a computer with a graphic card, you must use offscreen=2.
In this tutorial, we consider a Cartesian grid, make it rotate with T.rotate, and dump image to a file at each step. The movie must be finalized at the end. A wait time is added before finalizing to be sure that render is finished.

[Download python script].