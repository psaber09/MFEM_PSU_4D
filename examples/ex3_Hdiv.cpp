//                                MFEM Example 3
//
// Compile with: make ex3
//
// Sample runs:  ex3 -m ../data/star.mesh
//               ex3 -m ../data/beam-tri.mesh -o 2
//               ex3 -m ../data/beam-tet.mesh
//               ex3 -m ../data/beam-hex.mesh
//               ex3 -m ../data/beam-hex.mesh -o 2 -pa
//               ex3 -m ../data/escher.mesh
//               ex3 -m ../data/escher.mesh -o 2
//               ex3 -m ../data/fichera.mesh
//               ex3 -m ../data/fichera-q2.vtk
//               ex3 -m ../data/fichera-q3.mesh
//               ex3 -m ../data/square-disc-nurbs.mesh
//               ex3 -m ../data/beam-hex-nurbs.mesh
//               ex3 -m ../data/amr-hex.mesh
//               ex3 -m ../data/fichera-amr.mesh
//               ex3 -m ../data/star-surf.mesh -o 1
//               ex3 -m ../data/mobius-strip.mesh -f 0.1
//               ex3 -m ../data/klein-bottle.mesh -f 0.1
//
// Device sample runs:
//               ex3 -m ../data/star.mesh -pa -d cuda
//               ex3 -m ../data/star.mesh -pa -d raja-cuda
//               ex3 -m ../data/star.mesh -pa -d raja-omp
//               ex3 -m ../data/beam-hex.mesh -pa -d cuda
//
// Description:  This example code solves a simple electromagnetic diffusion
//               problem corresponding to the second order definite Maxwell
//               equation curl curl E + E = f with boundary condition
//               E x n = <given tangential field>. Here, we use a given exact
//               solution E and compute the corresponding r.h.s. f.
//               We discretize with Nedelec finite elements in 2D or 3D.
//
//               The example demonstrates the use of H(curl) finite element
//               spaces with the curl-curl and the (vector finite element) mass
//               bilinear form, as well as the computation of discretization
//               error when the exact solution is known. Static condensation is
//               also illustrated.
//
//               We recommend viewing examples 1-2 before viewing this example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
void E_exact(const Vector &, Vector &);
void Q_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
double freq = 1.0, kappa;
int dim;

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/beam-tet.mesh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = 1;
   int ref_levels = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                  " solution.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ref_levels, "-ref", "--ref-levels", "");

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);
   kappa = freq * M_PI;

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
      //int ref_levels = (int)floor(log(50000./mesh->GetNE())/log(2.)/dim);
       std::cout << "ref_levels = " << ref_levels << std::endl;
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }
   //mesh->ReorientTetMesh();

   // 5. Define a finite element space on the mesh. Here we use the Nedelec
   //    finite elements of the specified order.
   FiniteElementCollection *fec = new RT_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

   

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x by projecting the exact
   //    solution. Note that only values from the boundary edges will be used
   //    when eliminating the non-homogeneous boundary condition to modify the
   //    r.h.s. vector b.
   GridFunction x(fespace);
   //VectorFunctionCoefficient E(sdim, E_exact);
   //x.ProjectCoefficient(E);
    
   VectorFunctionCoefficient Q(sdim, Q_exact);
   x.ProjectCoefficient(Q);
    
    const IntegrationRule* irs[Geometry::NumGeom];
    for (int i = 0; i < Geometry::NumGeom; i++)
       irs[i] = &(IntRules.Get(i, 12));

     // 8a. Compute and print the L^2 norm of the error.
    cout << "\n Initial || E_h - E ||_{L^2} = " << x.ComputeL2Error(Q, irs) << '\n' << endl;
    cout << "Made it" << endl;


   return 0;
}


void E_exact(const Vector &x, Vector &E)
{
   if (dim == 3)
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(2));
      E(2) = sin(kappa * x(0));
   }
   else
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(0));
      if (x.Size() == 3) { E(2) = 0.0; }
   }
}

void Q_exact(const Vector &x, Vector &Q)
{
   if (dim == 4)
   {
      Q(0) = sin(kappa * x(1));
      Q(1) = sin(kappa * x(2));
      Q(2) = sin(kappa * x(3));
      Q(3) = sin(kappa * x(0));
//        Q(0) = x(0) * x(1) * x(3) * x(4);
//        Q(1) = x(0) * x(1) * x(3) * x(4);
//        Q(2) = x(0) * x(1) * x(3) * x(4);
//        Q(3) = x(0) * x(1) * x(3) * x(4);


   }
  if (dim == 3)
  {
      Q(0) = sin(kappa * x(1));
      Q(1) = sin(kappa * x(2));
      Q(2) = sin(kappa * x(0));
  }else
   {
      Q(0) = sin(kappa * x(1));
      Q(1) = sin(kappa * x(0));
      if (x.Size() == 3) { Q(2) = 0.0; }
   }
}

void f_exact(const Vector &x, Vector &f)
{
   if (dim == 3)
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(2));
      f(2) = (1. + kappa * kappa) * sin(kappa * x(0));
   }
   else
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(0));
      if (x.Size() == 3) { f(2) = 0.0; }
   }
}
