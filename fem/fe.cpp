// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Finite Element classes

#include "fe.hpp"
#include "fe_coll.hpp"
#include "../mesh/nurbs.hpp"
#include "bilininteg.hpp"
#include <cmath>

namespace mfem
{

using namespace std;

FiniteElement::FiniteElement(int D, Geometry::Type G, int Do, int O, int F)
   : Nodes(Do)
{
   dim = D ; geom_type = G ; dof = Do ; order = O ; func_space = F;
   range_type = SCALAR;
   map_type = VALUE;
   deriv_type = NONE;
   deriv_range_type = SCALAR;
   deriv_map_type = VALUE;
   for (int i = 0; i < Geometry::MaxDim; i++) { orders[i] = -1; }
#ifndef MFEM_THREAD_SAFE
   vshape.SetSize(dof, dim);
#endif
}

void FiniteElement::CalcVShape (
   const IntegrationPoint &ip, DenseMatrix &shape) const
{
   MFEM_ABORT("method is not implemented for this class");
}

void FiniteElement::CalcVShape (
   ElementTransformation &Trans, DenseMatrix &shape) const
{
   MFEM_ABORT("method is not implemented for this class");
}

void FiniteElement::CalcDivShape (
   const IntegrationPoint &ip, Vector &divshape) const
{
   MFEM_ABORT("method is not implemented for this class");
}

void FiniteElement::CalcPhysDivShape(
   ElementTransformation &Trans, Vector &div_shape) const
{
   CalcDivShape(Trans.GetIntPoint(), div_shape);
   div_shape *= (1.0 / Trans.Weight());
}

void FiniteElement::CalcCurlShape(const IntegrationPoint &ip,
                                  DenseMatrix &curl_shape) const
{
   MFEM_ABORT("method is not implemented for this class");
}

void FiniteElement::CalcPhysCurlShape(ElementTransformation &Trans,
                                      DenseMatrix &curl_shape) const
{
   switch (dim)
   {
      case 3:
      {
#ifdef MFEM_THREAD_SAFE
         DenseMatrix vshape(dof, dim);
#endif
         CalcCurlShape(Trans.GetIntPoint(), vshape);
         MultABt(vshape, Trans.Jacobian(), curl_shape);
         curl_shape *= (1.0 / Trans.Weight());
         break;
      }
      case 2:
         // This is valid for both 2x2 and 3x2 Jacobians
         CalcCurlShape(Trans.GetIntPoint(), curl_shape);
         curl_shape *= (1.0 / Trans.Weight());
         break;
      default:
         MFEM_ABORT("Invalid dimension, Dim = " << dim);
   }
}

void FiniteElement::CalcDivSkewShape(const IntegrationPoint &ip,
                                     DenseMatrix &curl_shape) const
{
   mfem_error ("FiniteElement::CalcDivSkewShape (ip, ...)\n"
               "   is not implemented for this class!");
}


void FiniteElement::GetFaceDofs(int face, int **dofs, int *ndofs) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::CalcHessian (const IntegrationPoint &ip,
                                 DenseMatrix &h) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::GetLocalInterpolation (ElementTransformation &Trans,
                                           DenseMatrix &I) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::GetLocalRestriction(ElementTransformation &,
                                        DenseMatrix &) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::GetTransferMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &I) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::Project (
   Coefficient &coeff, ElementTransformation &Trans, Vector &dofs) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::ProjectFromNodes(Vector &vc, ElementTransformation &Trans,
                                     Vector &dofs) const
{
   mfem_error ("FiniteElement::ProjectFromNodes() (vector) is not overloaded!");
}

void FiniteElement::ProjectMatrixCoefficient(
   MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const
{
   MFEM_ABORT("method is not overloaded");
}

void FiniteElement::ProjectDelta(int vertex, Vector &dofs) const
{
   MFEM_ABORT("method is not implemented for this element");
}

void FiniteElement::Project(
   const FiniteElement &fe, ElementTransformation &Trans, DenseMatrix &I) const
{
   MFEM_ABORT("method is not implemented for this element");
}

void FiniteElement::ProjectGrad(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &grad) const
{
   MFEM_ABORT("method is not implemented for this element");
}

void FiniteElement::ProjectCurl(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &curl) const
{
   MFEM_ABORT("method is not implemented for this element");
}

void FiniteElement::ProjectDivSkew(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &DivSkew) const
{
   mfem_error("FiniteElement::ProjectDivSkew(...) is not implemented for "
              "this element!");
}

void FiniteElement::ProjectDiv(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &div) const
{
   MFEM_ABORT("method is not implemented for this element");
}

void FiniteElement::CalcPhysShape(ElementTransformation &Trans,
                                  Vector &shape) const
{
   CalcShape(Trans.GetIntPoint(), shape);
   if (map_type == INTEGRAL)
   {
      shape /= Trans.Weight();
   }
}

void FiniteElement::CalcPhysDShape(ElementTransformation &Trans,
                                   DenseMatrix &dshape) const
{
   MFEM_ASSERT(map_type == VALUE, "");
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif
   CalcDShape(Trans.GetIntPoint(), vshape);
   Mult(vshape, Trans.InverseJacobian(), dshape);
}

void FiniteElement::CalcPhysLaplacian(ElementTransformation &Trans,
                                      Vector &Laplacian) const
{
   MFEM_ASSERT(map_type == VALUE, "");

   // Simpler routine if mapping is affine
   if (Trans.Hessian().FNorm2() < 1e-20)
   {
      CalcPhysLinLaplacian(Trans, Laplacian);
      return;
   }

   // Compute full Hessian first if non-affine
   int size = (dim*(dim+1))/2;
   DenseMatrix hess(dof, size);
   CalcPhysHessian(Trans,hess);

   if (dim == 3)
   {
      for (int nd = 0; nd < dof; nd++)
      {
         Laplacian[nd] = hess(nd,0) + hess(nd,4) + hess(nd,5);
      }
   }
   else if (dim == 2)
   {
      for (int nd = 0; nd < dof; nd++)
      {
         Laplacian[nd] = hess(nd,0) + hess(nd,2);
      }
   }
   else
   {
      for (int nd = 0; nd < dof; nd++)
      {
         Laplacian[nd] = hess(nd,0);
      }
   }
}


// Assume a linear mapping
void FiniteElement::CalcPhysLinLaplacian(ElementTransformation &Trans,
                                         Vector &Laplacian) const
{
   MFEM_ASSERT(map_type == VALUE, "");
   int size = (dim*(dim+1))/2;
   DenseMatrix hess(dof, size);
   DenseMatrix Gij(dim,dim);
   Vector scale(size);

   CalcHessian (Trans.GetIntPoint(), hess);
   MultAAt(Trans.InverseJacobian(), Gij);

   if (dim == 3)
   {
      scale[0] =   Gij(0,0);
      scale[1] = 2*Gij(0,1);
      scale[2] = 2*Gij(0,2);

      scale[3] = 2*Gij(1,2);
      scale[4] =   Gij(2,2);

      scale[5] =   Gij(1,1);
   }
   else if (dim == 2)
   {
      scale[0] =   Gij(0,0);
      scale[1] = 2*Gij(0,1);
      scale[2] =   Gij(1,1);
   }
   else
   {
      scale[0] =   Gij(0,0);
   }

   for (int nd = 0; nd < dof; nd++)
   {
      Laplacian[nd] = 0.0;
      for (int ii = 0; ii < size; ii++)
      {
         Laplacian[nd] += hess(nd,ii)*scale[ii];
      }
   }

}

void  FiniteElement::CalcPhysHessian(ElementTransformation &Trans,
                                     DenseMatrix& Hessian) const
{
   MFEM_ASSERT(map_type == VALUE, "");

   // Roll 2-Tensors in vectors and 4-Tensor in Matrix, exploiting symmetry
   Array<int> map(dim*dim);
   if (dim == 3)
   {
      map[0] = 0;
      map[1] = 1;
      map[2] = 2;

      map[3] = 1;
      map[4] = 5;
      map[5] = 3;

      map[6] = 2;
      map[7] = 3;
      map[8] = 4;
   }
   else if (dim == 2)
   {
      map[0] = 0;
      map[1] = 1;

      map[2] = 1;
      map[3] = 2;
   }
   else
   {
      map[0] = 0;
   }

   // Hessian in ref coords
   int size = (dim*(dim+1))/2;
   DenseMatrix hess(dof, size);
   CalcHessian(Trans.GetIntPoint(), hess);

   // Gradient in physical coords
   if (Trans.Hessian().FNorm2() > 1e-10)
   {
      DenseMatrix grad(dof, dim);
      CalcPhysDShape(Trans, grad);
      DenseMatrix gmap(dof, size);
      Mult(grad,Trans.Hessian(),gmap);
      hess -= gmap;
   }

   // LHM
   DenseMatrix lhm(size,size);
   DenseMatrix invJ = Trans.Jacobian();
   lhm = 0.0;
   for (int i = 0; i < dim; i++)
   {
      for (int j = 0; j < dim; j++)
      {
         for (int k = 0; k < dim; k++)
         {
            for (int l = 0; l < dim; l++)
            {
               lhm(map[i*dim+j],map[k*dim+l]) += invJ(i,k)*invJ(j,l);
            }
         }
      }
   }
   // Correct multiplicity
   Vector mult(size);
   mult = 0.0;
   for (int i = 0; i < dim*dim; i++) { mult[map[i]]++; }
   lhm.InvRightScaling(mult);

   // Hessian in physical coords
   lhm.Invert();
   Mult( hess, lhm, Hessian);
}

const DofToQuad &FiniteElement::GetDofToQuad(const IntegrationRule &,
                                             DofToQuad::Mode) const
{
   MFEM_ABORT("method is not implemented for this element");
   return *dof2quad_array[0]; // suppress a warning
}

FiniteElement::~FiniteElement()
{
   for (int i = 0; i < dof2quad_array.Size(); i++)
   {
      delete dof2quad_array[i];
   }
}


void ScalarFiniteElement::NodalLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I,
   const ScalarFiniteElement &fine_fe) const
{
   double v[Geometry::MaxDim];
   Vector vv (v, dim);
   IntegrationPoint f_ip;

#ifdef MFEM_THREAD_SAFE
   Vector c_shape(dof);
#endif

   MFEM_ASSERT(map_type == fine_fe.GetMapType(), "");

   I.SetSize(fine_fe.dof, dof);
   for (int i = 0; i < fine_fe.dof; i++)
   {
      Trans.Transform(fine_fe.Nodes.IntPoint(i), vv);
      f_ip.Set(v, dim);
      CalcShape(f_ip, c_shape);
      for (int j = 0; j < dof; j++)
      {
         if (fabs(I(i,j) = c_shape(j)) < 1.0e-12)
         {
            I(i,j) = 0.0;
         }
      }
   }
   if (map_type == INTEGRAL)
   {
      // assuming Trans is linear; this should be ok for all refinement types
      Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
      I *= Trans.Weight();
   }
}

void ScalarFiniteElement::ScalarLocalInterpolation(
   ElementTransformation &Trans, DenseMatrix &I,
   const ScalarFiniteElement &fine_fe) const
{
   // General "interpolation", defined by L2 projection

   double v[Geometry::MaxDim];
   Vector vv (v, dim);
   IntegrationPoint f_ip;

   const int fs = fine_fe.GetDof(), cs = this->GetDof();
   I.SetSize(fs, cs);
   Vector fine_shape(fs), coarse_shape(cs);
   DenseMatrix fine_mass(fs), fine_coarse_mass(fs, cs); // initialized with 0
   const int ir_order =
      std::max(GetOrder(), fine_fe.GetOrder()) + fine_fe.GetOrder();
   const IntegrationRule &ir = IntRules.Get(fine_fe.GetGeomType(), ir_order);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      fine_fe.CalcShape(ip, fine_shape);
      Trans.Transform(ip, vv);
      f_ip.Set(v, dim);
      this->CalcShape(f_ip, coarse_shape);

      AddMult_a_VVt(ip.weight, fine_shape, fine_mass);
      AddMult_a_VWt(ip.weight, fine_shape, coarse_shape, fine_coarse_mass);
   }

   DenseMatrixInverse fine_mass_inv(fine_mass);
   fine_mass_inv.Mult(fine_coarse_mass, I);

   if (map_type == INTEGRAL)
   {
      // assuming Trans is linear; this should be ok for all refinement types
      Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
      I *= Trans.Weight();
   }
}

void ScalarFiniteElement::ScalarLocalRestriction(
   ElementTransformation &Trans, DenseMatrix &R,
   const ScalarFiniteElement &coarse_fe) const
{
   // General "restriction", defined by L2 projection
   double v[Geometry::MaxDim];
   Vector vv (v, dim);
   IntegrationPoint f_ip;

   const int cs = coarse_fe.GetDof(), fs = this->GetDof();
   R.SetSize(cs, fs);
   Vector fine_shape(fs), coarse_shape(cs);
   DenseMatrix coarse_mass(cs), coarse_fine_mass(cs, fs); // initialized with 0
   const int ir_order = GetOrder() + coarse_fe.GetOrder();
   const IntegrationRule &ir = IntRules.Get(coarse_fe.GetGeomType(), ir_order);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      this->CalcShape(ip, fine_shape);
      Trans.Transform(ip, vv);
      f_ip.Set(v, dim);
      coarse_fe.CalcShape(f_ip, coarse_shape);

      AddMult_a_VVt(ip.weight, coarse_shape, coarse_mass);
      AddMult_a_VWt(ip.weight, coarse_shape, fine_shape, coarse_fine_mass);
   }

   DenseMatrixInverse coarse_mass_inv(coarse_mass);
   coarse_mass_inv.Mult(coarse_fine_mass, R);

   if (map_type == INTEGRAL)
   {
      // assuming Trans is linear; this should be ok for all refinement types
      Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
      R *= 1.0 / Trans.Weight();
   }
}
const DofToQuad &ScalarFiniteElement::GetDofToQuad(const IntegrationRule &ir,
                                                   DofToQuad::Mode mode) const
{
   MFEM_VERIFY(mode == DofToQuad::FULL, "invalid mode requested");

   for (int i = 0; i < dof2quad_array.Size(); i++)
   {
      const DofToQuad &d2q = *dof2quad_array[i];
      if (d2q.IntRule == &ir && d2q.mode == mode) { return d2q; }
   }

   DofToQuad *d2q = new DofToQuad;
   const int nqpt = ir.GetNPoints();
   d2q->FE = this;
   d2q->IntRule = &ir;
   d2q->mode = mode;
   d2q->ndof = dof;
   d2q->nqpt = nqpt;
   d2q->B.SetSize(nqpt*dof);
   d2q->Bt.SetSize(dof*nqpt);
   d2q->G.SetSize(nqpt*dim*dof);
   d2q->Gt.SetSize(dof*nqpt*dim);
#ifdef MFEM_THREAD_SAFE
   Vector c_shape(dof);
   DenseMatrix vshape(dof, dim);
#endif
   for (int i = 0; i < nqpt; i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      CalcShape(ip, c_shape);
      for (int j = 0; j < dof; j++)
      {
         d2q->B[i+nqpt*j] = d2q->Bt[j+dof*i] = c_shape(j);
      }
      CalcDShape(ip, vshape);
      for (int d = 0; d < dim; d++)
      {
         for (int j = 0; j < dof; j++)
         {
            d2q->G[i+nqpt*(d+dim*j)] = d2q->Gt[j+dof*(i+nqpt*d)] = vshape(j,d);
         }
      }
   }
   dof2quad_array.Append(d2q);
   return *d2q;
}

// protected method
const DofToQuad &ScalarFiniteElement::GetTensorDofToQuad(
   const TensorBasisElement &tb,
   const IntegrationRule &ir, DofToQuad::Mode mode) const
{
   MFEM_VERIFY(mode == DofToQuad::TENSOR, "invalid mode requested");

   for (int i = 0; i < dof2quad_array.Size(); i++)
   {
      const DofToQuad &d2q = *dof2quad_array[i];
      if (d2q.IntRule == &ir && d2q.mode == mode) { return d2q; }
   }

   DofToQuad *d2q = new DofToQuad;
   const Poly_1D::Basis &basis_1d = tb.GetBasis1D();
   const int ndof = order + 1;
   const int nqpt = (int)floor(pow(ir.GetNPoints(), 1.0/dim) + 0.5);
   d2q->FE = this;
   d2q->IntRule = &ir;
   d2q->mode = mode;
   d2q->ndof = ndof;
   d2q->nqpt = nqpt;
   d2q->B.SetSize(nqpt*ndof);
   d2q->Bt.SetSize(ndof*nqpt);
   d2q->G.SetSize(nqpt*ndof);
   d2q->Gt.SetSize(ndof*nqpt);
   Vector val(ndof), grad(ndof);
   for (int i = 0; i < nqpt; i++)
   {
      // The first 'nqpt' points in 'ir' have the same x-coordinates as those
      // of the 1D rule.
      basis_1d.Eval(ir.IntPoint(i).x, val, grad);
      for (int j = 0; j < ndof; j++)
      {
         d2q->B[i+nqpt*j] = d2q->Bt[j+ndof*i] = val(j);
         d2q->G[i+nqpt*j] = d2q->Gt[j+ndof*i] = grad(j);
      }
   }
   dof2quad_array.Append(d2q);
   return *d2q;
}


void NodalFiniteElement::ProjectCurl_2D(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &curl) const
{
   DenseMatrix curl_shape(fe.GetDof(), 1);

   curl.SetSize(dof, fe.GetDof());
   for (int i = 0; i < dof; i++)
   {
      fe.CalcCurlShape(Nodes.IntPoint(i), curl_shape);

      double w = 1.0;
      if (GetMapType() == FiniteElement::VALUE)
      {
         Trans.SetIntPoint(&Nodes.IntPoint(i));
         w /= Trans.Weight();
      }
      for (int j = 0; j < fe.GetDof(); j++)
      {
         curl(i,j) = w * curl_shape(j,0);
      }
   }
}

void InvertLinearTrans(ElementTransformation &trans,
                       const IntegrationPoint &pt, Vector &x)
{
   // invert a linear transform with one Newton step
   IntegrationPoint p0;
   p0.Set3(0, 0, 0);
   trans.Transform(p0, x);

   double store[3];
   Vector v(store, x.Size());
   pt.Get(v, x.Size());
   v -= x;

   trans.InverseJacobian().Mult(v, x);
}

void NodalFiniteElement::GetLocalRestriction(ElementTransformation &Trans,
                                             DenseMatrix &R) const
{
   IntegrationPoint ipt;
   Vector pt(&ipt.x, dim);

#ifdef MFEM_THREAD_SAFE
   Vector c_shape(dof);
#endif

   Trans.SetIntPoint(&Nodes[0]);

   for (int j = 0; j < dof; j++)
   {
      InvertLinearTrans(Trans, Nodes[j], pt);
      if (Geometries.CheckPoint(geom_type, ipt)) // do we need an epsilon here?
      {
         CalcShape(ipt, c_shape);
         R.SetRow(j, c_shape);
      }
      else
      {
         // Set the whole row to avoid valgrind warnings in R.Threshold().
         R.SetRow(j, infinity());
      }
   }
   R.Threshold(1e-12);
}

void NodalFiniteElement::Project (
   Coefficient &coeff, ElementTransformation &Trans, Vector &dofs) const
{
   for (int i = 0; i < dof; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      // some coefficients expect that Trans.IntPoint is the same
      // as the second argument of Eval
      Trans.SetIntPoint(&ip);
      dofs(i) = coeff.Eval (Trans, ip);
      if (map_type == INTEGRAL)
      {
         dofs(i) *= Trans.Weight();
      }
   }
}

void NodalFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   MFEM_ASSERT(dofs.Size() == vc.GetVDim()*dof, "");
   Vector x(vc.GetVDim());

   for (int i = 0; i < dof; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      Trans.SetIntPoint(&ip);
      vc.Eval (x, Trans, ip);
      if (map_type == INTEGRAL)
      {
         x *= Trans.Weight();
      }
      for (int j = 0; j < x.Size(); j++)
      {
         dofs(dof*j+i) = x(j);
      }
   }
}

void NodalFiniteElement::ProjectMatrixCoefficient(
   MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const
{
   // (mc.height x mc.width) @ DOFs -> (dof x mc.width x mc.height) in dofs
   MFEM_ASSERT(dofs.Size() == mc.GetHeight()*mc.GetWidth()*dof, "");
   DenseMatrix MQ(mc.GetHeight(), mc.GetWidth());

   for (int k = 0; k < dof; k++)
   {
      T.SetIntPoint(&Nodes.IntPoint(k));
      mc.Eval(MQ, T, Nodes.IntPoint(k));
      if (map_type == INTEGRAL) { MQ *= T.Weight(); }
      for (int r = 0; r < MQ.Height(); r++)
      {
         for (int d = 0; d < MQ.Width(); d++)
         {
            dofs(k+dof*(d+MQ.Width()*r)) = MQ(r,d);
         }
      }
   }
}

void NodalFiniteElement::Project(
   const FiniteElement &fe, ElementTransformation &Trans, DenseMatrix &I) const
{
   if (fe.GetRangeType() == SCALAR)
   {
      Vector shape(fe.GetDof());

      I.SetSize(dof, fe.GetDof());
      if (map_type == fe.GetMapType())
      {
         for (int k = 0; k < dof; k++)
         {
            fe.CalcShape(Nodes.IntPoint(k), shape);
            for (int j = 0; j < shape.Size(); j++)
            {
               I(k,j) = (fabs(shape(j)) < 1e-12) ? 0.0 : shape(j);
            }
         }
      }
      else
      {
         for (int k = 0; k < dof; k++)
         {
            Trans.SetIntPoint(&Nodes.IntPoint(k));
            fe.CalcPhysShape(Trans, shape);
            if (map_type == INTEGRAL)
            {
               shape *= Trans.Weight();
            }
            for (int j = 0; j < shape.Size(); j++)
            {
               I(k,j) = (fabs(shape(j)) < 1e-12) ? 0.0 : shape(j);
            }
         }
      }
   }
   else
   {
      DenseMatrix vshape(fe.GetDof(), Trans.GetSpaceDim());

      I.SetSize(vshape.Width()*dof, fe.GetDof());
      for (int k = 0; k < dof; k++)
      {
         Trans.SetIntPoint(&Nodes.IntPoint(k));
         fe.CalcVShape(Trans, vshape);
         if (map_type == INTEGRAL)
         {
            vshape *= Trans.Weight();
         }
         for (int j = 0; j < vshape.Height(); j++)
            for (int d = 0; d < vshape.Width(); d++)
            {
               I(k+d*dof,j) = vshape(j,d);
            }
      }
   }
}

void NodalFiniteElement::ProjectGrad(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &grad) const
{
   MFEM_ASSERT(fe.GetMapType() == VALUE, "");
   MFEM_ASSERT(Trans.GetSpaceDim() == dim, "")

   DenseMatrix dshape(fe.GetDof(), dim), grad_k(fe.GetDof(), dim), Jinv(dim);

   grad.SetSize(dim*dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(k);
      fe.CalcDShape(ip, dshape);
      Trans.SetIntPoint(&ip);
      CalcInverse(Trans.Jacobian(), Jinv);
      Mult(dshape, Jinv, grad_k);
      if (map_type == INTEGRAL)
      {
         grad_k *= Trans.Weight();
      }
      for (int j = 0; j < grad_k.Height(); j++)
         for (int d = 0; d < dim; d++)
         {
            grad(k+d*dof,j) = grad_k(j,d);
         }
   }
}

void NodalFiniteElement::ProjectDiv(
   const FiniteElement &fe, ElementTransformation &Trans,
   DenseMatrix &div) const
{
   double detJ;
   Vector div_shape(fe.GetDof());

   div.SetSize(dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(k);
      fe.CalcDivShape(ip, div_shape);
      if (map_type == VALUE)
      {
         Trans.SetIntPoint(&ip);
         detJ = Trans.Weight();
         for (int j = 0; j < div_shape.Size(); j++)
         {
            div(k,j) = (fabs(div_shape(j)) < 1e-12) ? 0.0 : div_shape(j)/detJ;
         }
      }
      else
      {
         for (int j = 0; j < div_shape.Size(); j++)
         {
            div(k,j) = (fabs(div_shape(j)) < 1e-12) ? 0.0 : div_shape(j);
         }
      }
   }
}


void PositiveFiniteElement::Project(
   Coefficient &coeff, ElementTransformation &Trans, Vector &dofs) const
{
   for (int i = 0; i < dof; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      Trans.SetIntPoint(&ip);
      dofs(i) = coeff.Eval(Trans, ip);
   }
}

void PositiveFiniteElement::Project(
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   MFEM_ASSERT(dofs.Size() == vc.GetVDim()*dof, "");
   Vector x(vc.GetVDim());

   for (int i = 0; i < dof; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      Trans.SetIntPoint(&ip);
      vc.Eval (x, Trans, ip);
      for (int j = 0; j < x.Size(); j++)
      {
         dofs(dof*j+i) = x(j);
      }
   }
}

void PositiveFiniteElement::Project(
   const FiniteElement &fe, ElementTransformation &Trans, DenseMatrix &I) const
{
   const NodalFiniteElement *nfe =
      dynamic_cast<const NodalFiniteElement *>(&fe);

   if (nfe && dof == nfe->GetDof())
   {
      nfe->Project(*this, Trans, I);
      I.Invert();
   }
   else
   {
      // local L2 projection
      DenseMatrix pos_mass, mixed_mass;
      MassIntegrator mass_integ;

      mass_integ.AssembleElementMatrix(*this, Trans, pos_mass);
      mass_integ.AssembleElementMatrix2(fe, *this, Trans, mixed_mass);

      DenseMatrixInverse pos_mass_inv(pos_mass);
      I.SetSize(dof, fe.GetDof());
      pos_mass_inv.Mult(mixed_mass, I);
   }
}


void VectorFiniteElement::CalcShape (
   const IntegrationPoint &ip, Vector &shape ) const
{
   mfem_error ("Error: Cannot use scalar CalcShape(...) function with\n"
               "   VectorFiniteElements!");
}

void VectorFiniteElement::CalcDShape (
   const IntegrationPoint &ip, DenseMatrix &dshape ) const
{
   mfem_error ("Error: Cannot use scalar CalcDShape(...) function with\n"
               "   VectorFiniteElements!");
}

void VectorFiniteElement::SetDerivMembers()
{
   switch (map_type)
   {
      case H_DIV:
         deriv_type = DIV;
         deriv_range_type = SCALAR;
         deriv_map_type = INTEGRAL;
         break;
      case H_DIV_SKEW:
         deriv_type = DIV_SKEW;
         deriv_range_type = VECTOR;
         deriv_map_type = H_DIV;
         break;
      case H_CURL:
         switch (dim)
         {
            case 4: // curl: 4D H_CURL -> 4D H_DIV(skew)
               deriv_type = CURL;
               deriv_range_type = MAT_SKEW;
               deriv_map_type = H_DIV_SKEW;
               break;
            case 3: // curl: 3D H_CURL -> 3D H_DIV
               deriv_type = CURL;
               deriv_range_type = VECTOR;
               deriv_map_type = H_DIV;
               break;
            case 2:
               // curl: 2D H_CURL -> INTEGRAL
               deriv_type = CURL;
               deriv_range_type = SCALAR;
               deriv_map_type = INTEGRAL;
               break;
            case 1:
               deriv_type = NONE;
               deriv_range_type = SCALAR;
               deriv_map_type = INTEGRAL;
               break;
            default:
               MFEM_ABORT("Invalid dimension, Dim = " << dim);
         }
         break;
      default:
         MFEM_ABORT("Invalid MapType = " << map_type);
   }
}

void VectorFiniteElement::CalcVShape_RT (
   ElementTransformation &Trans, DenseMatrix &shape) const
{
   MFEM_ASSERT(map_type == H_DIV, "");
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif
   CalcVShape(Trans.GetIntPoint(), vshape);
   MultABt(vshape, Trans.Jacobian(), shape);
   shape *= (1.0 / Trans.Weight());
}

void VectorFiniteElement::CalcVShape_DivSkew (
   ElementTransformation &Trans, DenseMatrix &shape) const
{
   if (dim!=4) { return; }

   MFEM_ASSERT(map_type == H_DIV_SKEW, "");
   const DenseMatrix &J = Trans.Jacobian();

#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim*dim);
   DenseMatrix Jinv(J.Width(), J.Height());
#else
   Jinv.SetSize(J.Width(), J.Height());
#endif

   if (vshape.Width()!=dim*dim) { vshape.SetSize(dof,dim*dim); }

   CalcVShape(Trans.GetIntPoint(), vshape);

   CalcInverse(J, Jinv);
   DenseMatrix invJtr(Jinv); invJtr.Transpose();

   CalcVShape(Trans.GetIntPoint(), vshape);

   DenseMatrix mat(dim,dim); mat = 0.0;
   DenseMatrix tempMat(dim,dim);

   for (int o=0; o<dof; o++)
   {
      //    for(int ik=0; ik<dim; ik++)
      //       for(int jk=0; jk<dim; jk++)
      //       {
      //          mat(ik,jk) = vshape(o,dim*ik+jk);
      //       }
      //
      //    Mult(mat, Jinv, tempMat);
      //    Mult(invJtr, tempMat, mat);
      //
      //    for(int ik=0; ik<dim; ik++)
      //       for(int jk=0; jk<dim; jk++)
      //       {
      //          shape(o,dim*ik+jk) = mat(ik,jk);
      //       }


      mat(0,0) =  0.0;            mat(0,1) =  vshape(o,11);
      mat(0,2) =  vshape(o,13); mat(0,3) =  vshape(o,6);
      mat(1,0) =  vshape(o,14); mat(1,1) =  0.0;
      mat(1,2) =  vshape(o,3);  mat(1,3) =  vshape(o,8);
      mat(2,0) =  vshape(o,7);  mat(2,1) =  vshape(o,12); mat(2,2) =  0.0;
      mat(2,3) =  vshape(o,1);
      mat(3,0) =  vshape(o,9);  mat(3,1) =  vshape(o,2);
      mat(3,2) =  vshape(o,4);  mat(3,3) =  0.0;

      Mult(mat, Jinv, tempMat);
      Mult(invJtr, tempMat, mat);

      shape(o,0) =  0.0;      shape(o,1) =  mat(2,3); shape(o,2) =  mat(3,1);
      shape(o,3) =  mat(1,2);
      shape(o,4) =  mat(3,2); shape(o,5) =  0.0;      shape(o,6) =  mat(0,3);
      shape(o,7) =  mat(2,0);
      shape(o,8) =  mat(1,3); shape(o,9) =  mat(3,0); shape(o,10) = 0.0;
      shape(o,11) = mat(0,1);
      shape(o,12) = mat(2,1); shape(o,13) = mat(0,2); shape(o,14) = mat(1,0);
      shape(o,15) = 0.0;
   }
}

void VectorFiniteElement::CalcVShape_ND (
   ElementTransformation &Trans, DenseMatrix &shape) const
{
   MFEM_ASSERT(map_type == H_CURL, "");
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif
   CalcVShape(Trans.GetIntPoint(), vshape);
   Mult(vshape, Trans.InverseJacobian(), shape);
}

void VectorFiniteElement::Project_RT(
   const double *nk, const Array<int> &d2n,
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   double vk[Geometry::MaxDim];
   const int sdim = Trans.GetSpaceDim();
   MFEM_ASSERT(vc.GetVDim() == sdim, "");
   Vector xk(vk, sdim);
   const bool square_J = (dim == sdim);

   for (int k = 0; k < dof; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      vc.Eval(xk, Trans, Nodes.IntPoint(k));
      // dof_k = nk^t adj(J) xk
      dofs(k) = Trans.AdjugateJacobian().InnerProduct(vk, nk + d2n[k]*dim);
      if (!square_J) { dofs(k) /= Trans.Weight(); }
   }
}

void VectorFiniteElement::Project_RT(
   const double *nk, const Array<int> &d2n,
   Vector &vc, ElementTransformation &Trans, Vector &dofs) const
{
   const int sdim = Trans.GetSpaceDim();
   const bool square_J = (dim == sdim);

   for (int k = 0; k < dof; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      // dof_k = nk^t adj(J) xk
      Vector vk(vc.GetData()+k*sdim, sdim);
      dofs(k) = Trans.AdjugateJacobian().InnerProduct(vk, nk + d2n[k]*dim);
      if (!square_J) { dofs(k) /= Trans.Weight(); }
   }
}

void VectorFiniteElement::ProjectMatrixCoefficient_RT(
   const double *nk, const Array<int> &d2n,
   MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const
{
   // project the rows of the matrix coefficient in an RT space

   const int sdim = T.GetSpaceDim();
   MFEM_ASSERT(mc.GetWidth() == sdim, "");
   const bool square_J = (dim == sdim);
   DenseMatrix MQ(mc.GetHeight(), mc.GetWidth());
   Vector nk_phys(sdim), dofs_k(MQ.Height());
   MFEM_ASSERT(dofs.Size() == dof*MQ.Height(), "");

   for (int k = 0; k < dof; k++)
   {
      T.SetIntPoint(&Nodes.IntPoint(k));
      mc.Eval(MQ, T, Nodes.IntPoint(k));
      // nk_phys = adj(J)^t nk
      T.AdjugateJacobian().MultTranspose(nk + d2n[k]*dim, nk_phys);
      if (!square_J) { nk_phys /= T.Weight(); }
      MQ.Mult(nk_phys, dofs_k);
      for (int r = 0; r < MQ.Height(); r++)
      {
         dofs(k+dof*r) = dofs_k(r);
      }
   }
}

void VectorFiniteElement::Project_RT(
   const double *nk, const Array<int> &d2n, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &I) const
{
   if (fe.GetRangeType() == SCALAR)
   {
      double vk[Geometry::MaxDim];
      Vector shape(fe.GetDof());
      int sdim = Trans.GetSpaceDim();

      I.SetSize(dof, sdim*fe.GetDof());
      for (int k = 0; k < dof; k++)
      {
         const IntegrationPoint &ip = Nodes.IntPoint(k);

         fe.CalcShape(ip, shape);
         Trans.SetIntPoint(&ip);
         // Transform RT face normals from reference to physical space
         // vk = adj(J)^T nk
         Trans.AdjugateJacobian().MultTranspose(nk + d2n[k]*dim, vk);
         if (fe.GetMapType() == INTEGRAL)
         {
            double w = 1.0/Trans.Weight();
            for (int d = 0; d < dim; d++)
            {
               vk[d] *= w;
            }
         }

         for (int j = 0; j < shape.Size(); j++)
         {
            double s = shape(j);
            if (fabs(s) < 1e-12)
            {
               s = 0.0;
            }
            // Project scalar basis function multiplied by each coordinate
            // direction onto the transformed face normals
            for (int d = 0; d < sdim; d++)
            {
               I(k,j+d*shape.Size()) = s*vk[d];
            }
         }
      }
   }
   else
   {
      int sdim = Trans.GetSpaceDim();
      double vk[Geometry::MaxDim];
      DenseMatrix vshape(fe.GetDof(), sdim);
      Vector vshapenk(fe.GetDof());
      const bool square_J = (dim == sdim);

      I.SetSize(dof, fe.GetDof());
      for (int k = 0; k < dof; k++)
      {
         const IntegrationPoint &ip = Nodes.IntPoint(k);

         Trans.SetIntPoint(&ip);
         // Transform RT face normals from reference to physical space
         // vk = adj(J)^T nk
         Trans.AdjugateJacobian().MultTranspose(nk + d2n[k]*dim, vk);
         // Compute fe basis functions in physical space
         fe.CalcVShape(Trans, vshape);
         // Project fe basis functions onto transformed face normals
         vshape.Mult(vk, vshapenk);
         if (!square_J) { vshapenk /= Trans.Weight(); }
         for (int j=0; j<vshapenk.Size(); j++)
         {
            I(k,j) = vshapenk(j);
         }
      }
   }
}

void VectorFiniteElement::ProjectGrad_RT(
   const double *nk, const Array<int> &d2n, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &grad) const
{
   if (dim != 2)
   {
      mfem_error("VectorFiniteElement::ProjectGrad_RT works only in 2D!");
   }

   DenseMatrix dshape(fe.GetDof(), fe.GetDim());
   Vector grad_k(fe.GetDof());
   double tk[2];

   grad.SetSize(dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      fe.CalcDShape(Nodes.IntPoint(k), dshape);
      tk[0] = nk[d2n[k]*dim+1];
      tk[1] = -nk[d2n[k]*dim];
      dshape.Mult(tk, grad_k);
      for (int j = 0; j < grad_k.Size(); j++)
      {
         grad(k,j) = (fabs(grad_k(j)) < 1e-12) ? 0.0 : grad_k(j);
      }
   }
}

void VectorFiniteElement::ProjectCurl_ND(
   const double *tk, const Array<int> &d2t, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &curl) const
{
#ifdef MFEM_THREAD_SAFE
   DenseMatrix curlshape(fe.GetDof(), dim);
   DenseMatrix curlshape_J(fe.GetDof(), dim);
   DenseMatrix J(dim, dim);
#else
   curlshape.SetSize(fe.GetDof(), dim);
   curlshape_J.SetSize(fe.GetDof(), dim);
   J.SetSize(dim, dim);
#endif

   Vector curl_k(fe.GetDof());

   curl.SetSize(dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(k);

      // calculate J^t * J / |J|
      Trans.SetIntPoint(&ip);
      MultAtB(Trans.Jacobian(), Trans.Jacobian(), J);
      J *= 1.0 / Trans.Weight();

      // transform curl of shapes (rows) by J^t * J / |J|
      fe.CalcCurlShape(ip, curlshape);
      Mult(curlshape, J, curlshape_J);

      curlshape_J.Mult(tk + d2t[k]*dim, curl_k);
      for (int j = 0; j < curl_k.Size(); j++)
      {
         curl(k,j) = (fabs(curl_k(j)) < 1e-12) ? 0.0 : curl_k(j);
      }
   }
}

void VectorFiniteElement::ProjectCurl_RT(
   const double *nk, const Array<int> &d2n, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &curl) const
{
   DenseMatrix curl_shape(fe.GetDof(), dim);
   Vector curl_k(fe.GetDof());

   curl.SetSize(dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      fe.CalcCurlShape(Nodes.IntPoint(k), curl_shape);
      curl_shape.Mult(nk + d2n[k]*dim, curl_k);
      for (int j = 0; j < curl_k.Size(); j++)
      {
         curl(k,j) = (fabs(curl_k(j)) < 1e-12) ? 0.0 : curl_k(j);
      }
   }
}

void VectorFiniteElement::Project_ND(
   const double *tk, const Array<int> &d2t,
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   double vk[Geometry::MaxDim];
   Vector xk(vk, vc.GetVDim());

   for (int k = 0; k < dof; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));

      vc.Eval(xk, Trans, Nodes.IntPoint(k));
      // dof_k = xk^t J tk
      dofs(k) = Trans.Jacobian().InnerProduct(tk + d2t[k]*dim, vk);
   }
}

void VectorFiniteElement::Project_ND(
   const double *tk, const Array<int> &d2t,
   Vector &vc, ElementTransformation &Trans, Vector &dofs) const
{
   for (int k = 0; k < dof; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      Vector vk(vc.GetData()+k*dim, dim);
      // dof_k = xk^t J tk
      dofs(k) = Trans.Jacobian().InnerProduct(tk + d2t[k]*dim, vk);
   }
}

void VectorFiniteElement::ProjectMatrixCoefficient_ND(
   const double *tk, const Array<int> &d2t,
   MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const
{
   // project the rows of the matrix coefficient in an ND space

   const int sdim = T.GetSpaceDim();
   MFEM_ASSERT(mc.GetWidth() == sdim, "");
   DenseMatrix MQ(mc.GetHeight(), mc.GetWidth());
   Vector tk_phys(sdim), dofs_k(MQ.Height());
   MFEM_ASSERT(dofs.Size() == dof*MQ.Height(), "");

   for (int k = 0; k < dof; k++)
   {
      T.SetIntPoint(&Nodes.IntPoint(k));
      mc.Eval(MQ, T, Nodes.IntPoint(k));
      // tk_phys = J tk
      T.Jacobian().Mult(tk + d2t[k]*dim, tk_phys);
      MQ.Mult(tk_phys, dofs_k);
      for (int r = 0; r < MQ.Height(); r++)
      {
         dofs(k+dof*r) = dofs_k(r);
      }
   }
}

void VectorFiniteElement::Project_ND(
   const double *tk, const Array<int> &d2t, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &I) const
{
   if (fe.GetRangeType() == SCALAR)
   {
      int sdim = Trans.GetSpaceDim();
      double vk[Geometry::MaxDim];
      Vector shape(fe.GetDof());

      I.SetSize(dof, sdim*fe.GetDof());
      for (int k = 0; k < dof; k++)
      {
         const IntegrationPoint &ip = Nodes.IntPoint(k);

         fe.CalcShape(ip, shape);
         Trans.SetIntPoint(&ip);
         // Transform ND edge tengents from reference to physical space
         // vk = J tk
         Trans.Jacobian().Mult(tk + d2t[k]*dim, vk);
         if (fe.GetMapType() == INTEGRAL)
         {
            double w = 1.0/Trans.Weight();
            for (int d = 0; d < sdim; d++)
            {
               vk[d] *= w;
            }
         }

         for (int j = 0; j < shape.Size(); j++)
         {
            double s = shape(j);
            if (fabs(s) < 1e-12)
            {
               s = 0.0;
            }
            // Project scalar basis function multiplied by each coordinate
            // direction onto the transformed edge tangents
            for (int d = 0; d < sdim; d++)
            {
               I(k, j + d*shape.Size()) = s*vk[d];
            }
         }
      }
   }
   else
   {
      int sdim = Trans.GetSpaceDim();
      double vk[Geometry::MaxDim];
      DenseMatrix vshape(fe.GetDof(), sdim);
      Vector vshapetk(fe.GetDof());

      I.SetSize(dof, fe.GetDof());
      for (int k = 0; k < dof; k++)
      {
         const IntegrationPoint &ip = Nodes.IntPoint(k);

         Trans.SetIntPoint(&ip);
         // Transform ND edge tangents from reference to physical space
         // vk = J tk
         Trans.Jacobian().Mult(tk + d2t[k]*dim, vk);
         // Compute fe basis functions in physical space
         fe.CalcVShape(Trans, vshape);
         // Project fe basis functions onto transformed edge tangents
         vshape.Mult(vk, vshapetk);
         for (int j=0; j<vshapetk.Size(); j++)
         {
            I(k, j) = vshapetk(j);
         }
      }
   }
}

void VectorFiniteElement::ProjectGrad_ND(
   const double *tk, const Array<int> &d2t, const FiniteElement &fe,
   ElementTransformation &Trans, DenseMatrix &grad) const
{
   MFEM_ASSERT(fe.GetMapType() == VALUE, "");

   DenseMatrix dshape(fe.GetDof(), fe.GetDim());
   Vector grad_k(fe.GetDof());

   grad.SetSize(dof, fe.GetDof());
   for (int k = 0; k < dof; k++)
   {
      fe.CalcDShape(Nodes.IntPoint(k), dshape);
      dshape.Mult(tk + d2t[k]*dim, grad_k);
      for (int j = 0; j < grad_k.Size(); j++)
      {
         grad(k,j) = (fabs(grad_k(j)) < 1e-12) ? 0.0 : grad_k(j);
      }
   }
}

void VectorFiniteElement::LocalInterpolation_RT(
   const VectorFiniteElement &cfe, const double *nk, const Array<int> &d2n,
   ElementTransformation &Trans, DenseMatrix &I) const
{
   MFEM_ASSERT(map_type == cfe.GetMapType(), "");

   double vk[Geometry::MaxDim];
   Vector xk(vk, dim);
   IntegrationPoint ip;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(cfe.GetDof(), cfe.GetDim());
#else
   DenseMatrix vshape(cfe.vshape.Data(), cfe.GetDof(), cfe.GetDim());
#endif
   I.SetSize(dof, vshape.Height());

   // assuming Trans is linear; this should be ok for all refinement types
   Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
   const DenseMatrix &adjJ = Trans.AdjugateJacobian();
   for (int k = 0; k < dof; k++)
   {
      Trans.Transform(Nodes.IntPoint(k), xk);
      ip.Set4(vk);
      cfe.CalcVShape(ip, vshape);
      // xk = |J| J^{-t} n_k
      adjJ.MultTranspose(nk + d2n[k]*dim, vk);
      // I_k = vshape_k.adj(J)^t.n_k, k=1,...,dof
      for (int j = 0; j < vshape.Height(); j++)
      {
         double Ikj = 0.;
         for (int i = 0; i < dim; i++)
         {
            Ikj += vshape(j, i) * vk[i];
         }
         I(k, j) = (fabs(Ikj) < 1e-12) ? 0.0 : Ikj;
      }
   }
}

void VectorFiniteElement::LocalInterpolation_ND(
   const VectorFiniteElement &cfe, const double *tk, const Array<int> &d2t,
   ElementTransformation &Trans, DenseMatrix &I) const
{
   double vk[Geometry::MaxDim];
   Vector xk(vk, dim);
   IntegrationPoint ip;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(cfe.GetDof(), cfe.GetDim());
#else
   DenseMatrix vshape(cfe.vshape.Data(), cfe.GetDof(), cfe.GetDim());
#endif
   I.SetSize(dof, vshape.Height());

   // assuming Trans is linear; this should be ok for all refinement types
   Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
   const DenseMatrix &J = Trans.Jacobian();
   for (int k = 0; k < dof; k++)
   {
      Trans.Transform(Nodes.IntPoint(k), xk);
      ip.Set3(vk);
      cfe.CalcVShape(ip, vshape);
      // xk = J t_k
      J.Mult(tk + d2t[k]*dim, vk);
      // I_k = vshape_k.J.t_k, k=1,...,Dof
      for (int j = 0; j < vshape.Height(); j++)
      {
         double Ikj = 0.;
         for (int i = 0; i < dim; i++)
         {
            Ikj += vshape(j, i) * vk[i];
         }
         I(k, j) = (fabs(Ikj) < 1e-12) ? 0.0 : Ikj;
      }
   }
}

void VectorFiniteElement::LocalRestriction_RT(
   const double *nk, const Array<int> &d2n, ElementTransformation &Trans,
   DenseMatrix &R) const
{
   double pt_data[Geometry::MaxDim];
   IntegrationPoint ip;
   Vector pt(pt_data, dim);

#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif

   Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
   const DenseMatrix &J = Trans.Jacobian();
   const double weight = Trans.Weight();
   for (int j = 0; j < dof; j++)
   {
      InvertLinearTrans(Trans, Nodes.IntPoint(j), pt);
      ip.Set(pt_data, dim);
      if (Geometries.CheckPoint(geom_type, ip)) // do we need an epsilon here?
      {
         CalcVShape(ip, vshape);
         J.MultTranspose(nk+dim*d2n[j], pt_data);
         pt /= weight;
         for (int k = 0; k < dof; k++)
         {
            double R_jk = 0.0;
            for (int d = 0; d < dim; d++)
            {
               R_jk += vshape(k,d)*pt_data[d];
            }
            R(j,k) = R_jk;
         }
      }
      else
      {
         // Set the whole row to avoid valgrind warnings in R.Threshold().
         R.SetRow(j, infinity());
      }
   }
   R.Threshold(1e-12);
}

void VectorFiniteElement::LocalRestriction_ND(
   const double *tk, const Array<int> &d2t, ElementTransformation &Trans,
   DenseMatrix &R) const
{
   double pt_data[Geometry::MaxDim];
   IntegrationPoint ip;
   Vector pt(pt_data, dim);

#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif

   Trans.SetIntPoint(&Geometries.GetCenter(geom_type));
   const DenseMatrix &Jinv = Trans.InverseJacobian();
   for (int j = 0; j < dof; j++)
   {
      InvertLinearTrans(Trans, Nodes.IntPoint(j), pt);
      ip.Set(pt_data, dim);
      if (Geometries.CheckPoint(geom_type, ip)) // do we need an epsilon here?
      {
         CalcVShape(ip, vshape);
         Jinv.Mult(tk+dim*d2t[j], pt_data);
         for (int k = 0; k < dof; k++)
         {
            double R_jk = 0.0;
            for (int d = 0; d < dim; d++)
            {
               R_jk += vshape(k,d)*pt_data[d];
            }
            R(j,k) = R_jk;
         }
      }
      else
      {
         // Set the whole row to avoid valgrind warnings in R.Threshold().
         R.SetRow(j, infinity());
      }
   }
   R.Threshold(1e-12);
}


PointFiniteElement::PointFiniteElement()
   : NodalFiniteElement(0, Geometry::POINT, 1, 0)
{
   lex_ordering.SetSize(1);
   lex_ordering[0] = 0;
   Nodes.IntPoint(0).x = 0.0;
}

void PointFiniteElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
   shape(0) = 1.;
}

void PointFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
   // dshape is (1 x 0) - nothing to compute
}

Linear1DFiniteElement::Linear1DFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 2, 1)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
}

void Linear1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   shape(0) = 1. - ip.x;
   shape(1) = ip.x;
}

void Linear1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   dshape(0,0) = -1.;
   dshape(1,0) =  1.;
}

Linear2DFiniteElement::Linear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 3, 1)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
}

void Linear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   shape(0) = 1. - ip.x - ip.y;
   shape(1) = ip.x;
   shape(2) = ip.y;
}

void Linear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   dshape(0,0) = -1.; dshape(0,1) = -1.;
   dshape(1,0) =  1.; dshape(1,1) =  0.;
   dshape(2,0) =  0.; dshape(2,1) =  1.;
}

BiLinear2DFiniteElement::BiLinear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 4, 1, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
}

void BiLinear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   shape(0) = (1. - ip.x) * (1. - ip.y) ;
   shape(1) = ip.x * (1. - ip.y) ;
   shape(2) = ip.x * ip.y ;
   shape(3) = (1. - ip.x) * ip.y ;
}

void BiLinear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   dshape(0,0) = -1. + ip.y; dshape(0,1) = -1. + ip.x ;
   dshape(1,0) =  1. - ip.y; dshape(1,1) = -ip.x ;
   dshape(2,0) =  ip.y ;     dshape(2,1) = ip.x ;
   dshape(3,0) = -ip.y ;     dshape(3,1) = 1. - ip.x ;
}

void BiLinear2DFiniteElement::CalcHessian(
   const IntegrationPoint &ip, DenseMatrix &h) const
{
   h(0,0) = 0.;   h(0,1) =  1.;   h(0,2) = 0.;
   h(1,0) = 0.;   h(1,1) = -1.;   h(1,2) = 0.;
   h(2,0) = 0.;   h(2,1) =  1.;   h(2,2) = 0.;
   h(3,0) = 0.;   h(3,1) = -1.;   h(3,2) = 0.;
}


GaussLinear2DFiniteElement::GaussLinear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 3, 1, FunctionSpace::Pk)
{
   Nodes.IntPoint(0).x = 1./6.;
   Nodes.IntPoint(0).y = 1./6.;
   Nodes.IntPoint(1).x = 2./3.;
   Nodes.IntPoint(1).y = 1./6.;
   Nodes.IntPoint(2).x = 1./6.;
   Nodes.IntPoint(2).y = 2./3.;
}

void GaussLinear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                           Vector &shape) const
{
   const double x = ip.x, y = ip.y;

   shape(0) = 5./3. - 2. * (x + y);
   shape(1) = 2. * (x - 1./6.);
   shape(2) = 2. * (y - 1./6.);
}

void GaussLinear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                            DenseMatrix &dshape) const
{
   dshape(0,0) = -2.;  dshape(0,1) = -2.;
   dshape(1,0) =  2.;  dshape(1,1) =  0.;
   dshape(2,0) =  0.;  dshape(2,1) =  2.;
}

void GaussLinear2DFiniteElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs(vertex)       = 2./3.;
   dofs((vertex+1)%3) = 1./6.;
   dofs((vertex+2)%3) = 1./6.;
}


// 0.5-0.5/sqrt(3) and 0.5+0.5/sqrt(3)
const double GaussBiLinear2DFiniteElement::p[] =
{ 0.2113248654051871177454256, 0.7886751345948128822545744 };

GaussBiLinear2DFiniteElement::GaussBiLinear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 4, 1, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = p[0];
   Nodes.IntPoint(0).y = p[0];
   Nodes.IntPoint(1).x = p[1];
   Nodes.IntPoint(1).y = p[0];
   Nodes.IntPoint(2).x = p[1];
   Nodes.IntPoint(2).y = p[1];
   Nodes.IntPoint(3).x = p[0];
   Nodes.IntPoint(3).y = p[1];
}

void GaussBiLinear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   const double x = ip.x, y = ip.y;

   shape(0) = 3. * (p[1] - x) * (p[1] - y);
   shape(1) = 3. * (x - p[0]) * (p[1] - y);
   shape(2) = 3. * (x - p[0]) * (y - p[0]);
   shape(3) = 3. * (p[1] - x) * (y - p[0]);
}

void GaussBiLinear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   const double x = ip.x, y = ip.y;

   dshape(0,0) = 3. * (y - p[1]);  dshape(0,1) = 3. * (x - p[1]);
   dshape(1,0) = 3. * (p[1] - y);  dshape(1,1) = 3. * (p[0] - x);
   dshape(2,0) = 3. * (y - p[0]);  dshape(2,1) = 3. * (x - p[0]);
   dshape(3,0) = 3. * (p[0] - y);  dshape(3,1) = 3. * (p[1] - x);
}

void GaussBiLinear2DFiniteElement::ProjectDelta(int vertex, Vector &dofs) const
{
#if 1
   dofs(vertex)       = p[1]*p[1];
   dofs((vertex+1)%4) = p[0]*p[1];
   dofs((vertex+2)%4) = p[0]*p[0];
   dofs((vertex+3)%4) = p[0]*p[1];
#else
   dofs = 1.0;
#endif
}


P1OnQuadFiniteElement::P1OnQuadFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 3, 1, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
}

void P1OnQuadFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   shape(0) = 1. - ip.x - ip.y;
   shape(1) = ip.x;
   shape(2) = ip.y;
}

void P1OnQuadFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   dshape(0,0) = -1.; dshape(0,1) = -1.;
   dshape(1,0) =  1.; dshape(1,1) =  0.;
   dshape(2,0) =  0.; dshape(2,1) =  1.;
}


Quad1DFiniteElement::Quad1DFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 3, 2)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(2).x = 0.5;
}

void Quad1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                    Vector &shape) const
{
   double x = ip.x;
   double l1 = 1.0 - x, l2 = x, l3 = 2. * x - 1.;

   shape(0) = l1 * (-l3);
   shape(1) = l2 * l3;
   shape(2) = 4. * l1 * l2;
}

void Quad1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   double x = ip.x;

   dshape(0,0) = 4. * x - 3.;
   dshape(1,0) = 4. * x - 1.;
   dshape(2,0) = 4. - 8. * x;
}


QuadPos1DFiniteElement::QuadPos1DFiniteElement()
   : PositiveFiniteElement(1, Geometry::SEGMENT, 3, 2)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(2).x = 0.5;
}

void QuadPos1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                       Vector &shape) const
{
   const double x = ip.x, x1 = 1. - x;

   shape(0) = x1 * x1;
   shape(1) = x * x;
   shape(2) = 2. * x * x1;
}

void QuadPos1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                        DenseMatrix &dshape) const
{
   const double x = ip.x;

   dshape(0,0) = 2. * x - 2.;
   dshape(1,0) = 2. * x;
   dshape(2,0) = 2. - 4. * x;
}

Quad2DFiniteElement::Quad2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 6, 2)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.5;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.5;
}

void Quad2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                    Vector &shape) const
{
   double x = ip.x, y = ip.y;
   double l1 = 1.-x-y, l2 = x, l3 = y;

   shape(0) = l1 * (2. * l1 - 1.);
   shape(1) = l2 * (2. * l2 - 1.);
   shape(2) = l3 * (2. * l3 - 1.);
   shape(3) = 4. * l1 * l2;
   shape(4) = 4. * l2 * l3;
   shape(5) = 4. * l3 * l1;
}

void Quad2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y;

   dshape(0,0) =
      dshape(0,1) = 4. * (x + y) - 3.;

   dshape(1,0) = 4. * x - 1.;
   dshape(1,1) = 0.;

   dshape(2,0) = 0.;
   dshape(2,1) = 4. * y - 1.;

   dshape(3,0) = -4. * (2. * x + y - 1.);
   dshape(3,1) = -4. * x;

   dshape(4,0) = 4. * y;
   dshape(4,1) = 4. * x;

   dshape(5,0) = -4. * y;
   dshape(5,1) = -4. * (x + 2. * y - 1.);
}

void Quad2DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                       DenseMatrix &h) const
{
   h(0,0) = 4.;
   h(0,1) = 4.;
   h(0,2) = 4.;

   h(1,0) = 4.;
   h(1,1) = 0.;
   h(1,2) = 0.;

   h(2,0) = 0.;
   h(2,1) = 0.;
   h(2,2) = 4.;

   h(3,0) = -8.;
   h(3,1) = -4.;
   h(3,2) =  0.;

   h(4,0) = 0.;
   h(4,1) = 4.;
   h(4,2) = 0.;

   h(5,0) =  0.;
   h(5,1) = -4.;
   h(5,2) = -8.;
}

void Quad2DFiniteElement::ProjectDelta(int vertex, Vector &dofs) const
{
#if 0
   dofs = 1.;
#else
   dofs = 0.;
   dofs(vertex) = 1.;
   switch (vertex)
   {
      case 0: dofs(3) = 0.25; dofs(5) = 0.25; break;
      case 1: dofs(3) = 0.25; dofs(4) = 0.25; break;
      case 2: dofs(4) = 0.25; dofs(5) = 0.25; break;
   }
#endif
}


const double GaussQuad2DFiniteElement::p[] =
{ 0.0915762135097707434595714634022015, 0.445948490915964886318329253883051 };

GaussQuad2DFiniteElement::GaussQuad2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 6, 2), A(6), D(6,2), pol(6)
{
   Nodes.IntPoint(0).x = p[0];
   Nodes.IntPoint(0).y = p[0];
   Nodes.IntPoint(1).x = 1. - 2. * p[0];
   Nodes.IntPoint(1).y = p[0];
   Nodes.IntPoint(2).x = p[0];
   Nodes.IntPoint(2).y = 1. - 2. * p[0];
   Nodes.IntPoint(3).x = p[1];
   Nodes.IntPoint(3).y = p[1];
   Nodes.IntPoint(4).x = 1. - 2. * p[1];
   Nodes.IntPoint(4).y = p[1];
   Nodes.IntPoint(5).x = p[1];
   Nodes.IntPoint(5).y = 1. - 2. * p[1];

   for (int i = 0; i < 6; i++)
   {
      const double x = Nodes.IntPoint(i).x, y = Nodes.IntPoint(i).y;
      A(0,i) = 1.;
      A(1,i) = x;
      A(2,i) = y;
      A(3,i) = x * x;
      A(4,i) = x * y;
      A(5,i) = y * y;
   }

   A.Invert();
}

void GaussQuad2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   const double x = ip.x, y = ip.y;
   pol(0) = 1.;
   pol(1) = x;
   pol(2) = y;
   pol(3) = x * x;
   pol(4) = x * y;
   pol(5) = y * y;

   A.Mult(pol, shape);
}

void GaussQuad2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   const double x = ip.x, y = ip.y;
   D(0,0) = 0.;      D(0,1) = 0.;
   D(1,0) = 1.;      D(1,1) = 0.;
   D(2,0) = 0.;      D(2,1) = 1.;
   D(3,0) = 2. *  x; D(3,1) = 0.;
   D(4,0) = y;       D(4,1) = x;
   D(5,0) = 0.;      D(5,1) = 2. * y;

   Mult(A, D, dshape);
}


BiQuad2DFiniteElement::BiQuad2DFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 9, 2, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.5;
}

void BiQuad2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   double x = ip.x, y = ip.y;
   double l1x, l2x, l3x, l1y, l2y, l3y;

   l1x = (x - 1.) * (2. * x - 1);
   l2x = 4. * x * (1. - x);
   l3x = x * (2. * x - 1.);
   l1y = (y - 1.) * (2. * y - 1);
   l2y = 4. * y * (1. - y);
   l3y = y * (2. * y - 1.);

   shape(0) = l1x * l1y;
   shape(4) = l2x * l1y;
   shape(1) = l3x * l1y;
   shape(7) = l1x * l2y;
   shape(8) = l2x * l2y;
   shape(5) = l3x * l2y;
   shape(3) = l1x * l3y;
   shape(6) = l2x * l3y;
   shape(2) = l3x * l3y;
}

void BiQuad2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y;
   double l1x, l2x, l3x, l1y, l2y, l3y;
   double d1x, d2x, d3x, d1y, d2y, d3y;

   l1x = (x - 1.) * (2. * x - 1);
   l2x = 4. * x * (1. - x);
   l3x = x * (2. * x - 1.);
   l1y = (y - 1.) * (2. * y - 1);
   l2y = 4. * y * (1. - y);
   l3y = y * (2. * y - 1.);

   d1x = 4. * x - 3.;
   d2x = 4. - 8. * x;
   d3x = 4. * x - 1.;
   d1y = 4. * y - 3.;
   d2y = 4. - 8. * y;
   d3y = 4. * y - 1.;

   dshape(0,0) = d1x * l1y;
   dshape(0,1) = l1x * d1y;

   dshape(4,0) = d2x * l1y;
   dshape(4,1) = l2x * d1y;

   dshape(1,0) = d3x * l1y;
   dshape(1,1) = l3x * d1y;

   dshape(7,0) = d1x * l2y;
   dshape(7,1) = l1x * d2y;

   dshape(8,0) = d2x * l2y;
   dshape(8,1) = l2x * d2y;

   dshape(5,0) = d3x * l2y;
   dshape(5,1) = l3x * d2y;

   dshape(3,0) = d1x * l3y;
   dshape(3,1) = l1x * d3y;

   dshape(6,0) = d2x * l3y;
   dshape(6,1) = l2x * d3y;

   dshape(2,0) = d3x * l3y;
   dshape(2,1) = l3x * d3y;
}

void BiQuad2DFiniteElement::ProjectDelta(int vertex, Vector &dofs) const
{
#if 0
   dofs = 1.;
#else
   dofs = 0.;
   dofs(vertex) = 1.;
   switch (vertex)
   {
      case 0: dofs(4) = 0.25; dofs(7) = 0.25; break;
      case 1: dofs(4) = 0.25; dofs(5) = 0.25; break;
      case 2: dofs(5) = 0.25; dofs(6) = 0.25; break;
      case 3: dofs(6) = 0.25; dofs(7) = 0.25; break;
   }
   dofs(8) = 1./16.;
#endif
}


H1Ser_QuadrilateralElement::H1Ser_QuadrilateralElement(const int p)
   : ScalarFiniteElement(2, Geometry::SQUARE, (p*p + 3*p +6) / 2, p,
                         FunctionSpace::Qk)
{
   // Store the dof_map of the associated TensorBasisElement, which will be used
   // to create the serendipity dof map.  Its size is larger than the size of
   // the serendipity element.
   TensorBasisElement tbeTemp =
      TensorBasisElement(2, p, BasisType::GaussLobatto,
                         TensorBasisElement::DofMapType::Sr_DOF_MAP);
   const Array<int> tp_dof_map = tbeTemp.GetDofMap();

   const double *cp = poly1d.ClosedPoints(p, BasisType::GaussLobatto);

   // Fixing the Nodes is exactly the same as the H1_QuadrilateralElement
   // constructor except we only use those values of the associated tensor
   // product dof_map that are <= the number of serendipity Dofs e.g. only DoFs
   // 0-7 out of the 9 tensor product dofs (at quadratic order)
   int o = 0;

   for (int j = 0; j <= p; j++)
   {
      for (int i = 0; i <= p; i++)
      {
         if (tp_dof_map[o] < Nodes.Size())
         {
            Nodes.IntPoint(tp_dof_map[o]).x = cp[i];
            Nodes.IntPoint(tp_dof_map[o]).y = cp[j];
         }
         o++;
      }
   }
}

void H1Ser_QuadrilateralElement::CalcShape(const IntegrationPoint &ip,
                                           Vector &shape) const
{
   int p = (this)->GetOrder();
   double x = ip.x, y = ip.y;

   Poly_1D::Basis edgeNodalBasis(poly1d.GetBasis(p, BasisType::GaussLobatto));
   Vector nodalX(p+1);
   Vector nodalY(p+1);

   edgeNodalBasis.Eval(x, nodalX);
   edgeNodalBasis.Eval(y, nodalY);

   // First, fix edge-based shape functions. Use a nodal interpolant for edge
   // points, weighted by the linear function that vanishes on opposite edge.
   for (int i = 0; i < p-1; i++)
   {
      shape(4 + 0*(p-1) + i) = (nodalX(i+1))*(1.-y);         // south edge 0->1
      shape(4 + 1*(p-1) + i) = (nodalY(i+1))*x;              // east edge  1->2
      shape(4 + 3*(p-1) - i - 1) = (nodalX(i+1)) * y;        // north edge 3->2
      shape(4 + 4*(p-1) - i - 1) = (nodalY(i+1)) * (1. - x); // west edge  0->3
   }

   BiLinear2DFiniteElement bilinear = BiLinear2DFiniteElement();
   Vector bilinearsAtIP(4);
   bilinear.CalcShape(ip, bilinearsAtIP);

   const double *edgePts(poly1d.ClosedPoints(p, BasisType::GaussLobatto));

   // Next, set the shape function associated with vertex V, evaluated at (x,y)
   // to be: bilinear function associated to V, evaluated at (x,y) - sum (shape
   // function at edge point P, weighted by bilinear function for V evaluated at
   // P) where the sum is taken only for points P on edges incident to V.

   double vtx0fix =0;
   double vtx1fix =0;
   double vtx2fix =0;
   double vtx3fix =0;
   for (int i = 0; i<p-1; i++)
   {
      vtx0fix += (1-edgePts[i+1])*(shape(4 + i) +
                                   shape(4 + 4*(p-1) - i - 1)); // bot+left edge
      vtx1fix += (1-edgePts[i+1])*(shape(4 + 1*(p-1) + i) +
                                   shape(4 + (p-2)-i));        // right+bot edge
      vtx2fix += (1-edgePts[i+1])*(shape(4 + 2*(p-1) + i) +
                                   shape(1 + 2*p-i));          // top+right edge
      vtx3fix += (1-edgePts[i+1])*(shape(4 + 3*(p-1) + i) +
                                   shape(3*p - i));            // left+top edge
   }
   shape(0) = bilinearsAtIP(0) - vtx0fix;
   shape(1) = bilinearsAtIP(1) - vtx1fix;
   shape(2) = bilinearsAtIP(2) - vtx2fix;
   shape(3) = bilinearsAtIP(3) - vtx3fix;

   // Interior basis functions appear starting at order p=4. These are non-nodal
   // bubble functions.
   if (p > 3)
   {
      double *legX = new double[p-1];
      double *legY = new double[p-1];
      Poly_1D *storeLegendre = new Poly_1D();

      storeLegendre->CalcLegendre(p-2, x, legX);
      storeLegendre->CalcLegendre(p-2, y, legY);

      int interior_total = 0;
      for (int j = 4; j < p + 1; j++)
      {
         for (int k = 0; k < j-3; k++)
         {
            shape(4 + 4*(p-1) + interior_total)
               = legX[k] * legY[j-4-k] * x * (1. - x) * y * (1. - y);
            interior_total++;
         }
      }

      delete[] legX;
      delete[] legY;
      delete storeLegendre;
   }
}

void H1Ser_QuadrilateralElement::CalcDShape(const IntegrationPoint &ip,
                                            DenseMatrix &dshape) const
{
   int p = (this)->GetOrder();
   double x = ip.x, y = ip.y;

   Poly_1D::Basis edgeNodalBasis(poly1d.GetBasis(p, BasisType::GaussLobatto));
   Vector nodalX(p+1);
   Vector DnodalX(p+1);
   Vector nodalY(p+1);
   Vector DnodalY(p+1);

   edgeNodalBasis.Eval(x, nodalX, DnodalX);
   edgeNodalBasis.Eval(y, nodalY, DnodalY);

   for (int i = 0; i < p-1; i++)
   {
      dshape(4 + 0*(p-1) + i,0) =  DnodalX(i+1) * (1.-y);
      dshape(4 + 0*(p-1) + i,1) = -nodalX(i+1);
      dshape(4 + 1*(p-1) + i,0) =  nodalY(i+1);
      dshape(4 + 1*(p-1) + i,1) =  DnodalY(i+1)*x;
      dshape(4 + 3*(p-1) - i - 1,0) =  DnodalX(i+1)*y;
      dshape(4 + 3*(p-1) - i - 1,1) =  nodalX(i+1);
      dshape(4 + 4*(p-1) - i - 1,0) = -nodalY(i+1);
      dshape(4 + 4*(p-1) - i - 1,1) =  DnodalY(i+1) * (1.-x);
   }

   BiLinear2DFiniteElement bilinear = BiLinear2DFiniteElement();
   DenseMatrix DbilinearsAtIP(4);
   bilinear.CalcDShape(ip, DbilinearsAtIP);

   const double *edgePts(poly1d.ClosedPoints(p, BasisType::GaussLobatto));

   dshape(0,0) = DbilinearsAtIP(0,0);
   dshape(0,1) = DbilinearsAtIP(0,1);
   dshape(1,0) = DbilinearsAtIP(1,0);
   dshape(1,1) = DbilinearsAtIP(1,1);
   dshape(2,0) = DbilinearsAtIP(2,0);
   dshape(2,1) = DbilinearsAtIP(2,1);
   dshape(3,0) = DbilinearsAtIP(3,0);
   dshape(3,1) = DbilinearsAtIP(3,1);

   for (int i = 0; i<p-1; i++)
   {
      dshape(0,0) -= (1-edgePts[i+1])*(dshape(4 + 0*(p-1) + i, 0) +
                                       dshape(4 + 4*(p-1) - i - 1,0));
      dshape(0,1) -= (1-edgePts[i+1])*(dshape(4 + 0*(p-1) + i, 1) +
                                       dshape(4 + 4*(p-1) - i - 1,1));
      dshape(1,0) -= (1-edgePts[i+1])*(dshape(4 + 1*(p-1) + i, 0) +
                                       dshape(4 + (p-2)-i, 0));
      dshape(1,1) -= (1-edgePts[i+1])*(dshape(4 + 1*(p-1) + i, 1) +
                                       dshape(4 + (p-2)-i, 1));
      dshape(2,0) -= (1-edgePts[i+1])*(dshape(4 + 2*(p-1) + i, 0) +
                                       dshape(1 + 2*p-i, 0));
      dshape(2,1) -= (1-edgePts[i+1])*(dshape(4 + 2*(p-1) + i, 1) +
                                       dshape(1 + 2*p-i, 1));
      dshape(3,0) -= (1-edgePts[i+1])*(dshape(4 + 3*(p-1) + i, 0) +
                                       dshape(3*p - i, 0));
      dshape(3,1) -= (1-edgePts[i+1])*(dshape(4 + 3*(p-1) + i, 1) +
                                       dshape(3*p - i, 1));
   }

   if (p > 3)
   {
      double *legX = new double[p-1];
      double *legY = new double[p-1];
      double *DlegX = new double[p-1];
      double *DlegY = new double[p-1];
      Poly_1D *storeLegendre = new Poly_1D();

      storeLegendre->CalcLegendre(p-2, x, legX, DlegX);
      storeLegendre->CalcLegendre(p-2, y, legY, DlegY);

      int interior_total = 0;
      for (int j = 4; j < p + 1; j++)
      {
         for (int k = 0; k < j-3; k++)
         {
            dshape(4 + 4*(p-1) + interior_total, 0) =
               legY[j-4-k]*y*(1-y) * (DlegX[k]*x*(1-x) + legX[k]*(1-2*x));
            dshape(4 + 4*(p-1) + interior_total, 1) =
               legX[k]*x*(1-x) * (DlegY[j-4-k]*y*(1-y) + legY[j-4-k]*(1-2*y));
            interior_total++;
         }
      }
      delete[] legX;
      delete[] legY;
      delete[] DlegX;
      delete[] DlegY;
      delete storeLegendre;
   }
}

void H1Ser_QuadrilateralElement::GetLocalInterpolation(ElementTransformation
                                                       &Trans,
                                                       DenseMatrix &I) const
{
   // For p<=4, the basis is nodal; for p>4, the quad-interior functions are
   // non-nodal.
   if (order <= 4)
   {
      NodalLocalInterpolation(Trans, I, *this);
   }
   else
   {
      ScalarLocalInterpolation(Trans, I, *this);
   }
}


BiQuadPos2DFiniteElement::BiQuadPos2DFiniteElement()
   : PositiveFiniteElement(2, Geometry::SQUARE, 9, 2, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.5;
}

void BiQuadPos2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   double x = ip.x, y = ip.y;
   double l1x, l2x, l3x, l1y, l2y, l3y;

   l1x = (1. - x) * (1. - x);
   l2x = 2. * x * (1. - x);
   l3x = x * x;
   l1y = (1. - y) * (1. - y);
   l2y = 2. * y * (1. - y);
   l3y = y * y;

   shape(0) = l1x * l1y;
   shape(4) = l2x * l1y;
   shape(1) = l3x * l1y;
   shape(7) = l1x * l2y;
   shape(8) = l2x * l2y;
   shape(5) = l3x * l2y;
   shape(3) = l1x * l3y;
   shape(6) = l2x * l3y;
   shape(2) = l3x * l3y;
}

void BiQuadPos2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y;
   double l1x, l2x, l3x, l1y, l2y, l3y;
   double d1x, d2x, d3x, d1y, d2y, d3y;

   l1x = (1. - x) * (1. - x);
   l2x = 2. * x * (1. - x);
   l3x = x * x;
   l1y = (1. - y) * (1. - y);
   l2y = 2. * y * (1. - y);
   l3y = y * y;

   d1x = 2. * x - 2.;
   d2x = 2. - 4. * x;
   d3x = 2. * x;
   d1y = 2. * y - 2.;
   d2y = 2. - 4. * y;
   d3y = 2. * y;

   dshape(0,0) = d1x * l1y;
   dshape(0,1) = l1x * d1y;

   dshape(4,0) = d2x * l1y;
   dshape(4,1) = l2x * d1y;

   dshape(1,0) = d3x * l1y;
   dshape(1,1) = l3x * d1y;

   dshape(7,0) = d1x * l2y;
   dshape(7,1) = l1x * d2y;

   dshape(8,0) = d2x * l2y;
   dshape(8,1) = l2x * d2y;

   dshape(5,0) = d3x * l2y;
   dshape(5,1) = l3x * d2y;

   dshape(3,0) = d1x * l3y;
   dshape(3,1) = l1x * d3y;

   dshape(6,0) = d2x * l3y;
   dshape(6,1) = l2x * d3y;

   dshape(2,0) = d3x * l3y;
   dshape(2,1) = l3x * d3y;
}

void BiQuadPos2DFiniteElement::GetLocalInterpolation(
   ElementTransformation &Trans, DenseMatrix &I) const
{
   double s[9];
   IntegrationPoint tr_ip;
   Vector xx(&tr_ip.x, 2), shape(s, 9);

   for (int i = 0; i < 9; i++)
   {
      Trans.Transform(Nodes.IntPoint(i), xx);
      CalcShape(tr_ip, shape);
      for (int j = 0; j < 9; j++)
         if (fabs(I(i,j) = s[j]) < 1.0e-12)
         {
            I(i,j) = 0.0;
         }
   }
   for (int i = 0; i < 9; i++)
   {
      double *d = &I(0,i);
      d[4] = 2. * d[4] - 0.5 * (d[0] + d[1]);
      d[5] = 2. * d[5] - 0.5 * (d[1] + d[2]);
      d[6] = 2. * d[6] - 0.5 * (d[2] + d[3]);
      d[7] = 2. * d[7] - 0.5 * (d[3] + d[0]);
      d[8] = 4. * d[8] - 0.5 * (d[4] + d[5] + d[6] + d[7]) -
             0.25 * (d[0] + d[1] + d[2] + d[3]);
   }
}

void BiQuadPos2DFiniteElement::Project(
   Coefficient &coeff, ElementTransformation &Trans, Vector &dofs) const
{
   double *d = dofs;

   for (int i = 0; i < 9; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      Trans.SetIntPoint(&ip);
      d[i] = coeff.Eval(Trans, ip);
   }
   d[4] = 2. * d[4] - 0.5 * (d[0] + d[1]);
   d[5] = 2. * d[5] - 0.5 * (d[1] + d[2]);
   d[6] = 2. * d[6] - 0.5 * (d[2] + d[3]);
   d[7] = 2. * d[7] - 0.5 * (d[3] + d[0]);
   d[8] = 4. * d[8] - 0.5 * (d[4] + d[5] + d[6] + d[7]) -
          0.25 * (d[0] + d[1] + d[2] + d[3]);
}

void BiQuadPos2DFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double v[3];
   Vector x (v, vc.GetVDim());

   for (int i = 0; i < 9; i++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(i);
      Trans.SetIntPoint(&ip);
      vc.Eval (x, Trans, ip);
      for (int j = 0; j < x.Size(); j++)
      {
         dofs(9*j+i) = v[j];
      }
   }
   for (int j = 0; j < x.Size(); j++)
   {
      double *d = &dofs(9*j);

      d[4] = 2. * d[4] - 0.5 * (d[0] + d[1]);
      d[5] = 2. * d[5] - 0.5 * (d[1] + d[2]);
      d[6] = 2. * d[6] - 0.5 * (d[2] + d[3]);
      d[7] = 2. * d[7] - 0.5 * (d[3] + d[0]);
      d[8] = 4. * d[8] - 0.5 * (d[4] + d[5] + d[6] + d[7]) -
             0.25 * (d[0] + d[1] + d[2] + d[3]);
   }
}


GaussBiQuad2DFiniteElement::GaussBiQuad2DFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 9, 2, FunctionSpace::Qk)
{
   const double p1 = 0.5*(1.-sqrt(3./5.));

   Nodes.IntPoint(0).x = p1;
   Nodes.IntPoint(0).y = p1;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = p1;
   Nodes.IntPoint(1).x = 1.-p1;
   Nodes.IntPoint(1).y = p1;
   Nodes.IntPoint(7).x = p1;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.5;
   Nodes.IntPoint(5).x = 1.-p1;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(3).x = p1;
   Nodes.IntPoint(3).y = 1.-p1;
   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 1.-p1;
   Nodes.IntPoint(2).x = 1.-p1;
   Nodes.IntPoint(2).y = 1.-p1;
}

void GaussBiQuad2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                           Vector &shape) const
{
   const double a = sqrt(5./3.);
   const double p1 = 0.5*(1.-sqrt(3./5.));

   double x = a*(ip.x-p1), y = a*(ip.y-p1);
   double l1x, l2x, l3x, l1y, l2y, l3y;

   l1x = (x - 1.) * (2. * x - 1);
   l2x = 4. * x * (1. - x);
   l3x = x * (2. * x - 1.);
   l1y = (y - 1.) * (2. * y - 1);
   l2y = 4. * y * (1. - y);
   l3y = y * (2. * y - 1.);

   shape(0) = l1x * l1y;
   shape(4) = l2x * l1y;
   shape(1) = l3x * l1y;
   shape(7) = l1x * l2y;
   shape(8) = l2x * l2y;
   shape(5) = l3x * l2y;
   shape(3) = l1x * l3y;
   shape(6) = l2x * l3y;
   shape(2) = l3x * l3y;
}

void GaussBiQuad2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                            DenseMatrix &dshape) const
{
   const double a = sqrt(5./3.);
   const double p1 = 0.5*(1.-sqrt(3./5.));

   double x = a*(ip.x-p1), y = a*(ip.y-p1);
   double l1x, l2x, l3x, l1y, l2y, l3y;
   double d1x, d2x, d3x, d1y, d2y, d3y;

   l1x = (x - 1.) * (2. * x - 1);
   l2x = 4. * x * (1. - x);
   l3x = x * (2. * x - 1.);
   l1y = (y - 1.) * (2. * y - 1);
   l2y = 4. * y * (1. - y);
   l3y = y * (2. * y - 1.);

   d1x = a * (4. * x - 3.);
   d2x = a * (4. - 8. * x);
   d3x = a * (4. * x - 1.);
   d1y = a * (4. * y - 3.);
   d2y = a * (4. - 8. * y);
   d3y = a * (4. * y - 1.);

   dshape(0,0) = d1x * l1y;
   dshape(0,1) = l1x * d1y;

   dshape(4,0) = d2x * l1y;
   dshape(4,1) = l2x * d1y;

   dshape(1,0) = d3x * l1y;
   dshape(1,1) = l3x * d1y;

   dshape(7,0) = d1x * l2y;
   dshape(7,1) = l1x * d2y;

   dshape(8,0) = d2x * l2y;
   dshape(8,1) = l2x * d2y;

   dshape(5,0) = d3x * l2y;
   dshape(5,1) = l3x * d2y;

   dshape(3,0) = d1x * l3y;
   dshape(3,1) = l1x * d3y;

   dshape(6,0) = d2x * l3y;
   dshape(6,1) = l2x * d3y;

   dshape(2,0) = d3x * l3y;
   dshape(2,1) = l3x * d3y;
}

BiCubic2DFiniteElement::BiCubic2DFiniteElement()
   : NodalFiniteElement (2, Geometry::SQUARE, 16, 3, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.;
   Nodes.IntPoint(0).y = 0.;
   Nodes.IntPoint(1).x = 1.;
   Nodes.IntPoint(1).y = 0.;
   Nodes.IntPoint(2).x = 1.;
   Nodes.IntPoint(2).y = 1.;
   Nodes.IntPoint(3).x = 0.;
   Nodes.IntPoint(3).y = 1.;
   Nodes.IntPoint(4).x = 1./3.;
   Nodes.IntPoint(4).y = 0.;
   Nodes.IntPoint(5).x = 2./3.;
   Nodes.IntPoint(5).y = 0.;
   Nodes.IntPoint(6).x = 1.;
   Nodes.IntPoint(6).y = 1./3.;
   Nodes.IntPoint(7).x = 1.;
   Nodes.IntPoint(7).y = 2./3.;
   Nodes.IntPoint(8).x = 2./3.;
   Nodes.IntPoint(8).y = 1.;
   Nodes.IntPoint(9).x = 1./3.;
   Nodes.IntPoint(9).y = 1.;
   Nodes.IntPoint(10).x = 0.;
   Nodes.IntPoint(10).y = 2./3.;
   Nodes.IntPoint(11).x = 0.;
   Nodes.IntPoint(11).y = 1./3.;
   Nodes.IntPoint(12).x = 1./3.;
   Nodes.IntPoint(12).y = 1./3.;
   Nodes.IntPoint(13).x = 2./3.;
   Nodes.IntPoint(13).y = 1./3.;
   Nodes.IntPoint(14).x = 1./3.;
   Nodes.IntPoint(14).y = 2./3.;
   Nodes.IntPoint(15).x = 2./3.;
   Nodes.IntPoint(15).y = 2./3.;
}

void BiCubic2DFiniteElement::CalcShape(
   const IntegrationPoint &ip, Vector &shape) const
{
   double x = ip.x, y = ip.y;

   double w1x, w2x, w3x, w1y, w2y, w3y;
   double l0x, l1x, l2x, l3x, l0y, l1y, l2y, l3y;

   w1x = x - 1./3.; w2x = x - 2./3.; w3x = x - 1.;
   w1y = y - 1./3.; w2y = y - 2./3.; w3y = y - 1.;

   l0x = (- 4.5) * w1x * w2x * w3x;
   l1x = ( 13.5) *   x * w2x * w3x;
   l2x = (-13.5) *   x * w1x * w3x;
   l3x = (  4.5) *   x * w1x * w2x;

   l0y = (- 4.5) * w1y * w2y * w3y;
   l1y = ( 13.5) *   y * w2y * w3y;
   l2y = (-13.5) *   y * w1y * w3y;
   l3y = (  4.5) *   y * w1y * w2y;

   shape(0)  = l0x * l0y;
   shape(1)  = l3x * l0y;
   shape(2)  = l3x * l3y;
   shape(3)  = l0x * l3y;
   shape(4)  = l1x * l0y;
   shape(5)  = l2x * l0y;
   shape(6)  = l3x * l1y;
   shape(7)  = l3x * l2y;
   shape(8)  = l2x * l3y;
   shape(9)  = l1x * l3y;
   shape(10) = l0x * l2y;
   shape(11) = l0x * l1y;
   shape(12) = l1x * l1y;
   shape(13) = l2x * l1y;
   shape(14) = l1x * l2y;
   shape(15) = l2x * l2y;
}

void BiCubic2DFiniteElement::CalcDShape(
   const IntegrationPoint &ip, DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y;

   double w1x, w2x, w3x, w1y, w2y, w3y;
   double l0x, l1x, l2x, l3x, l0y, l1y, l2y, l3y;
   double d0x, d1x, d2x, d3x, d0y, d1y, d2y, d3y;

   w1x = x - 1./3.; w2x = x - 2./3.; w3x = x - 1.;
   w1y = y - 1./3.; w2y = y - 2./3.; w3y = y - 1.;

   l0x = (- 4.5) * w1x * w2x * w3x;
   l1x = ( 13.5) *   x * w2x * w3x;
   l2x = (-13.5) *   x * w1x * w3x;
   l3x = (  4.5) *   x * w1x * w2x;

   l0y = (- 4.5) * w1y * w2y * w3y;
   l1y = ( 13.5) *   y * w2y * w3y;
   l2y = (-13.5) *   y * w1y * w3y;
   l3y = (  4.5) *   y * w1y * w2y;

   d0x = -5.5 + ( 18. - 13.5 * x) * x;
   d1x =  9.  + (-45. + 40.5 * x) * x;
   d2x = -4.5 + ( 36. - 40.5 * x) * x;
   d3x =  1.  + (- 9. + 13.5 * x) * x;

   d0y = -5.5 + ( 18. - 13.5 * y) * y;
   d1y =  9.  + (-45. + 40.5 * y) * y;
   d2y = -4.5 + ( 36. - 40.5 * y) * y;
   d3y =  1.  + (- 9. + 13.5 * y) * y;

   dshape( 0,0) = d0x * l0y;   dshape( 0,1) = l0x * d0y;
   dshape( 1,0) = d3x * l0y;   dshape( 1,1) = l3x * d0y;
   dshape( 2,0) = d3x * l3y;   dshape( 2,1) = l3x * d3y;
   dshape( 3,0) = d0x * l3y;   dshape( 3,1) = l0x * d3y;
   dshape( 4,0) = d1x * l0y;   dshape( 4,1) = l1x * d0y;
   dshape( 5,0) = d2x * l0y;   dshape( 5,1) = l2x * d0y;
   dshape( 6,0) = d3x * l1y;   dshape( 6,1) = l3x * d1y;
   dshape( 7,0) = d3x * l2y;   dshape( 7,1) = l3x * d2y;
   dshape( 8,0) = d2x * l3y;   dshape( 8,1) = l2x * d3y;
   dshape( 9,0) = d1x * l3y;   dshape( 9,1) = l1x * d3y;
   dshape(10,0) = d0x * l2y;   dshape(10,1) = l0x * d2y;
   dshape(11,0) = d0x * l1y;   dshape(11,1) = l0x * d1y;
   dshape(12,0) = d1x * l1y;   dshape(12,1) = l1x * d1y;
   dshape(13,0) = d2x * l1y;   dshape(13,1) = l2x * d1y;
   dshape(14,0) = d1x * l2y;   dshape(14,1) = l1x * d2y;
   dshape(15,0) = d2x * l2y;   dshape(15,1) = l2x * d2y;
}

void BiCubic2DFiniteElement::CalcHessian(
   const IntegrationPoint &ip, DenseMatrix &h) const
{
   double x = ip.x, y = ip.y;

   double w1x, w2x, w3x, w1y, w2y, w3y;
   double l0x, l1x, l2x, l3x, l0y, l1y, l2y, l3y;
   double d0x, d1x, d2x, d3x, d0y, d1y, d2y, d3y;
   double h0x, h1x, h2x, h3x, h0y, h1y, h2y, h3y;

   w1x = x - 1./3.; w2x = x - 2./3.; w3x = x - 1.;
   w1y = y - 1./3.; w2y = y - 2./3.; w3y = y - 1.;

   l0x = (- 4.5) * w1x * w2x * w3x;
   l1x = ( 13.5) *   x * w2x * w3x;
   l2x = (-13.5) *   x * w1x * w3x;
   l3x = (  4.5) *   x * w1x * w2x;

   l0y = (- 4.5) * w1y * w2y * w3y;
   l1y = ( 13.5) *   y * w2y * w3y;
   l2y = (-13.5) *   y * w1y * w3y;
   l3y = (  4.5) *   y * w1y * w2y;

   d0x = -5.5 + ( 18. - 13.5 * x) * x;
   d1x =  9.  + (-45. + 40.5 * x) * x;
   d2x = -4.5 + ( 36. - 40.5 * x) * x;
   d3x =  1.  + (- 9. + 13.5 * x) * x;

   d0y = -5.5 + ( 18. - 13.5 * y) * y;
   d1y =  9.  + (-45. + 40.5 * y) * y;
   d2y = -4.5 + ( 36. - 40.5 * y) * y;
   d3y =  1.  + (- 9. + 13.5 * y) * y;

   h0x = -27. * x + 18.;
   h1x =  81. * x - 45.;
   h2x = -81. * x + 36.;
   h3x =  27. * x -  9.;

   h0y = -27. * y + 18.;
   h1y =  81. * y - 45.;
   h2y = -81. * y + 36.;
   h3y =  27. * y -  9.;

   h( 0,0) = h0x * l0y;   h( 0,1) = d0x * d0y;   h( 0,2) = l0x * h0y;
   h( 1,0) = h3x * l0y;   h( 1,1) = d3x * d0y;   h( 1,2) = l3x * h0y;
   h( 2,0) = h3x * l3y;   h( 2,1) = d3x * d3y;   h( 2,2) = l3x * h3y;
   h( 3,0) = h0x * l3y;   h( 3,1) = d0x * d3y;   h( 3,2) = l0x * h3y;
   h( 4,0) = h1x * l0y;   h( 4,1) = d1x * d0y;   h( 4,2) = l1x * h0y;
   h( 5,0) = h2x * l0y;   h( 5,1) = d2x * d0y;   h( 5,2) = l2x * h0y;
   h( 6,0) = h3x * l1y;   h( 6,1) = d3x * d1y;   h( 6,2) = l3x * h1y;
   h( 7,0) = h3x * l2y;   h( 7,1) = d3x * d2y;   h( 7,2) = l3x * h2y;
   h( 8,0) = h2x * l3y;   h( 8,1) = d2x * d3y;   h( 8,2) = l2x * h3y;
   h( 9,0) = h1x * l3y;   h( 9,1) = d1x * d3y;   h( 9,2) = l1x * h3y;
   h(10,0) = h0x * l2y;   h(10,1) = d0x * d2y;   h(10,2) = l0x * h2y;
   h(11,0) = h0x * l1y;   h(11,1) = d0x * d1y;   h(11,2) = l0x * h1y;
   h(12,0) = h1x * l1y;   h(12,1) = d1x * d1y;   h(12,2) = l1x * h1y;
   h(13,0) = h2x * l1y;   h(13,1) = d2x * d1y;   h(13,2) = l2x * h1y;
   h(14,0) = h1x * l2y;   h(14,1) = d1x * d2y;   h(14,2) = l1x * h2y;
   h(15,0) = h2x * l2y;   h(15,1) = d2x * d2y;   h(15,2) = l2x * h2y;
}


Cubic1DFiniteElement::Cubic1DFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 4, 3)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(2).x = 0.33333333333333333333;
   Nodes.IntPoint(3).x = 0.66666666666666666667;
}

void Cubic1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   double x = ip.x;
   double l1 = x,
          l2 = (1.0-x),
          l3 = (0.33333333333333333333-x),
          l4 = (0.66666666666666666667-x);

   shape(0) =   4.5 * l2 * l3 * l4;
   shape(1) =   4.5 * l1 * l3 * l4;
   shape(2) =  13.5 * l1 * l2 * l4;
   shape(3) = -13.5 * l1 * l2 * l3;
}

void Cubic1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   double x = ip.x;

   dshape(0,0) = -5.5 + x * (18. - 13.5 * x);
   dshape(1,0) = 1. - x * (9. - 13.5 * x);
   dshape(2,0) = 9. - x * (45. - 40.5 * x);
   dshape(3,0) = -4.5 + x * (36. - 40.5 * x);
}


Cubic2DFiniteElement::Cubic2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 10, 3)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.33333333333333333333;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(4).x = 0.66666666666666666667;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(5).x = 0.66666666666666666667;
   Nodes.IntPoint(5).y = 0.33333333333333333333;
   Nodes.IntPoint(6).x = 0.33333333333333333333;
   Nodes.IntPoint(6).y = 0.66666666666666666667;
   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.66666666666666666667;
   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = 0.33333333333333333333;
   Nodes.IntPoint(9).x = 0.33333333333333333333;
   Nodes.IntPoint(9).y = 0.33333333333333333333;
}

void Cubic2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   double x = ip.x, y = ip.y;
   double l1 = (-1. + x + y),
          lx = (-1. + 3.*x),
          ly = (-1. + 3.*y);

   shape(0) = -0.5*l1*(3.*l1 + 1.)*(3.*l1 + 2.);
   shape(1) =  0.5*x*(lx - 1.)*lx;
   shape(2) =  0.5*y*(-1. + ly)*ly;
   shape(3) =  4.5*x*l1*(3.*l1 + 1.);
   shape(4) = -4.5*x*lx*l1;
   shape(5) =  4.5*x*lx*y;
   shape(6) =  4.5*x*y*ly;
   shape(7) = -4.5*y*l1*ly;
   shape(8) =  4.5*y*l1*(1. + 3.*l1);
   shape(9) = -27.*x*y*l1;
}

void Cubic2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y;

   dshape(0,0) =  0.5*(-11. + 36.*y - 9.*(x*(-4. + 3.*x) + 6.*x*y + 3.*y*y));
   dshape(1,0) =  1. + 4.5*x*(-2. + 3.*x);
   dshape(2,0) =  0.;
   dshape(3,0) =  4.5*(2. + 9.*x*x - 5.*y + 3.*y*y + 2.*x*(-5. + 6.*y));
   dshape(4,0) = -4.5*(1. - 1.*y + x*(-8. + 9.*x + 6.*y));
   dshape(5,0) =  4.5*(-1. + 6.*x)*y;
   dshape(6,0) =  4.5*y*(-1. + 3.*y);
   dshape(7,0) =  4.5*(1. - 3.*y)*y;
   dshape(8,0) =  4.5*y*(-5. + 6.*x + 6.*y);
   dshape(9,0) =  -27.*y*(-1. + 2.*x + y);

   dshape(0,1) =  0.5*(-11. + 36.*y - 9.*(x*(-4. + 3.*x) + 6.*x*y + 3.*y*y));
   dshape(1,1) =  0.;
   dshape(2,1) =  1. + 4.5*y*(-2. + 3.*y);
   dshape(3,1) =  4.5*x*(-5. + 6.*x + 6.*y);
   dshape(4,1) =  4.5*(1. - 3.*x)*x;
   dshape(5,1) =  4.5*x*(-1. + 3.*x);
   dshape(6,1) =  4.5*x*(-1. + 6.*y);
   dshape(7,1) = -4.5*(1. + x*(-1. + 6.*y) + y*(-8. + 9.*y));
   dshape(8,1) =  4.5*(2. + 3.*x*x + y*(-10. + 9.*y) + x*(-5. + 12.*y));
   dshape(9,1) = -27.*x*(-1. + x + 2.*y);
}

void Cubic2DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                        DenseMatrix &h) const
{
   double x = ip.x, y = ip.y;

   h(0,0) = 18.-27.*(x+y);
   h(0,1) = 18.-27.*(x+y);
   h(0,2) = 18.-27.*(x+y);

   h(1,0) = -9.+27.*x;
   h(1,1) = 0.;
   h(1,2) = 0.;

   h(2,0) = 0.;
   h(2,1) = 0.;
   h(2,2) = -9.+27.*y;

   h(3,0) = -45.+81.*x+54.*y;
   h(3,1) = -22.5+54.*x+27.*y;
   h(3,2) = 27.*x;

   h(4,0) = 36.-81.*x-27.*y;
   h(4,1) = 4.5-27.*x;
   h(4,2) = 0.;

   h(5,0) = 27.*y;
   h(5,1) = -4.5+27.*x;
   h(5,2) = 0.;

   h(6,0) = 0.;
   h(6,1) = -4.5+27.*y;
   h(6,2) = 27.*x;

   h(7,0) = 0.;
   h(7,1) = 4.5-27.*y;
   h(7,2) = 36.-27.*x-81.*y;

   h(8,0) = 27.*y;
   h(8,1) = -22.5+27.*x+54.*y;
   h(8,2) = -45.+54.*x+81.*y;

   h(9,0) = -54.*y;
   h(9,1) = 27.-54.*(x+y);
   h(9,2) = -54.*x;
}


Cubic3DFiniteElement::Cubic3DFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 20, 3)
{
   Nodes.IntPoint(0).x = 0;
   Nodes.IntPoint(0).y = 0;
   Nodes.IntPoint(0).z = 0;
   Nodes.IntPoint(1).x = 1.;
   Nodes.IntPoint(1).y = 0;
   Nodes.IntPoint(1).z = 0;
   Nodes.IntPoint(2).x = 0;
   Nodes.IntPoint(2).y = 1.;
   Nodes.IntPoint(2).z = 0;
   Nodes.IntPoint(3).x = 0;
   Nodes.IntPoint(3).y = 0;
   Nodes.IntPoint(3).z = 1.;
   Nodes.IntPoint(4).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(4).y = 0;
   Nodes.IntPoint(4).z = 0;
   Nodes.IntPoint(5).x = 0.6666666666666666666666666667;
   Nodes.IntPoint(5).y = 0;
   Nodes.IntPoint(5).z = 0;
   Nodes.IntPoint(6).x = 0;
   Nodes.IntPoint(6).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(6).z = 0;
   Nodes.IntPoint(7).x = 0;
   Nodes.IntPoint(7).y = 0.6666666666666666666666666667;
   Nodes.IntPoint(7).z = 0;
   Nodes.IntPoint(8).x = 0;
   Nodes.IntPoint(8).y = 0;
   Nodes.IntPoint(8).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(9).x = 0;
   Nodes.IntPoint(9).y = 0;
   Nodes.IntPoint(9).z = 0.6666666666666666666666666667;
   Nodes.IntPoint(10).x = 0.6666666666666666666666666667;
   Nodes.IntPoint(10).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(10).z = 0;
   Nodes.IntPoint(11).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(11).y = 0.6666666666666666666666666667;
   Nodes.IntPoint(11).z = 0;
   Nodes.IntPoint(12).x = 0.6666666666666666666666666667;
   Nodes.IntPoint(12).y = 0;
   Nodes.IntPoint(12).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(13).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(13).y = 0;
   Nodes.IntPoint(13).z = 0.6666666666666666666666666667;
   Nodes.IntPoint(14).x = 0;
   Nodes.IntPoint(14).y = 0.6666666666666666666666666667;
   Nodes.IntPoint(14).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(15).x = 0;
   Nodes.IntPoint(15).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(15).z = 0.6666666666666666666666666667;
   Nodes.IntPoint(16).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(16).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(16).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(17).x = 0;
   Nodes.IntPoint(17).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(17).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(18).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(18).y = 0;
   Nodes.IntPoint(18).z = 0.3333333333333333333333333333;
   Nodes.IntPoint(19).x = 0.3333333333333333333333333333;
   Nodes.IntPoint(19).y = 0.3333333333333333333333333333;
   Nodes.IntPoint(19).z = 0;
}

void Cubic3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;

   shape(0) = -((-1 + x + y + z)*(-2 + 3*x + 3*y + 3*z)*
                (-1 + 3*x + 3*y + 3*z))/2.;
   shape(4) = (9*x*(-1 + x + y + z)*(-2 + 3*x + 3*y + 3*z))/2.;
   shape(5) = (-9*x*(-1 + 3*x)*(-1 + x + y + z))/2.;
   shape(1) = (x*(2 + 9*(-1 + x)*x))/2.;
   shape(6) = (9*y*(-1 + x + y + z)*(-2 + 3*x + 3*y + 3*z))/2.;
   shape(19) = -27*x*y*(-1 + x + y + z);
   shape(10) = (9*x*(-1 + 3*x)*y)/2.;
   shape(7) = (-9*y*(-1 + 3*y)*(-1 + x + y + z))/2.;
   shape(11) = (9*x*y*(-1 + 3*y))/2.;
   shape(2) = (y*(2 + 9*(-1 + y)*y))/2.;
   shape(8) = (9*z*(-1 + x + y + z)*(-2 + 3*x + 3*y + 3*z))/2.;
   shape(18) = -27*x*z*(-1 + x + y + z);
   shape(12) = (9*x*(-1 + 3*x)*z)/2.;
   shape(17) = -27*y*z*(-1 + x + y + z);
   shape(16) = 27*x*y*z;
   shape(14) = (9*y*(-1 + 3*y)*z)/2.;
   shape(9) = (-9*z*(-1 + x + y + z)*(-1 + 3*z))/2.;
   shape(13) = (9*x*z*(-1 + 3*z))/2.;
   shape(15) = (9*y*z*(-1 + 3*z))/2.;
   shape(3) = (z*(2 + 9*(-1 + z)*z))/2.;
}

void Cubic3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y, z = ip.z;

   dshape(0,0) = (-11 + 36*y + 36*z - 9*(3*pow(x,2) + 3*pow(y + z,2) +
                                         x*(-4 + 6*y + 6*z)))/2.;
   dshape(0,1) = (-11 + 36*y + 36*z - 9*(3*pow(x,2) + 3*pow(y + z,2) +
                                         x*(-4 + 6*y + 6*z)))/2.;
   dshape(0,2) = (-11 + 36*y + 36*z - 9*(3*pow(x,2) + 3*pow(y + z,2) +
                                         x*(-4 + 6*y + 6*z)))/2.;
   dshape(4,0) = (9*(9*pow(x,2) + (-1 + y + z)*(-2 + 3*y + 3*z) +
                     2*x*(-5 + 6*y + 6*z)))/2.;
   dshape(4,1) = (9*x*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(4,2) = (9*x*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(5,0) = (-9*(1 - y - z + x*(-8 + 9*x + 6*y + 6*z)))/2.;
   dshape(5,1) = (9*(1 - 3*x)*x)/2.;
   dshape(5,2) = (9*(1 - 3*x)*x)/2.;
   dshape(1,0) = 1 + (9*x*(-2 + 3*x))/2.;
   dshape(1,1) = 0;
   dshape(1,2) = 0;
   dshape(6,0) = (9*y*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(6,1) = (9*(2 + 3*pow(x,2) - 10*y - 5*z + 3*(y + z)*(3*y + z) +
                     x*(-5 + 12*y + 6*z)))/2.;
   dshape(6,2) = (9*y*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(19,0) = -27*y*(-1 + 2*x + y + z);
   dshape(19,1) = -27*x*(-1 + x + 2*y + z);
   dshape(19,2) = -27*x*y;
   dshape(10,0) = (9*(-1 + 6*x)*y)/2.;
   dshape(10,1) = (9*x*(-1 + 3*x))/2.;
   dshape(10,2) = 0;
   dshape(7,0) = (9*(1 - 3*y)*y)/2.;
   dshape(7,1) = (-9*(1 + x*(-1 + 6*y) - z + y*(-8 + 9*y + 6*z)))/2.;
   dshape(7,2) = (9*(1 - 3*y)*y)/2.;
   dshape(11,0) = (9*y*(-1 + 3*y))/2.;
   dshape(11,1) = (9*x*(-1 + 6*y))/2.;
   dshape(11,2) = 0;
   dshape(2,0) = 0;
   dshape(2,1) = 1 + (9*y*(-2 + 3*y))/2.;
   dshape(2,2) = 0;
   dshape(8,0) = (9*z*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(8,1) = (9*z*(-5 + 6*x + 6*y + 6*z))/2.;
   dshape(8,2) = (9*(2 + 3*pow(x,2) - 5*y - 10*z + 3*(y + z)*(y + 3*z) +
                     x*(-5 + 6*y + 12*z)))/2.;
   dshape(18,0) = -27*z*(-1 + 2*x + y + z);
   dshape(18,1) = -27*x*z;
   dshape(18,2) = -27*x*(-1 + x + y + 2*z);
   dshape(12,0) = (9*(-1 + 6*x)*z)/2.;
   dshape(12,1) = 0;
   dshape(12,2) = (9*x*(-1 + 3*x))/2.;
   dshape(17,0) = -27*y*z;
   dshape(17,1) = -27*z*(-1 + x + 2*y + z);
   dshape(17,2) = -27*y*(-1 + x + y + 2*z);
   dshape(16,0) = 27*y*z;
   dshape(16,1) = 27*x*z;
   dshape(16,2) = 27*x*y;
   dshape(14,0) = 0;
   dshape(14,1) = (9*(-1 + 6*y)*z)/2.;
   dshape(14,2) = (9*y*(-1 + 3*y))/2.;
   dshape(9,0) = (9*(1 - 3*z)*z)/2.;
   dshape(9,1) = (9*(1 - 3*z)*z)/2.;
   dshape(9,2) = (9*(-1 + x + y + 8*z - 6*(x + y)*z - 9*pow(z,2)))/2.;
   dshape(13,0) = (9*z*(-1 + 3*z))/2.;
   dshape(13,1) = 0;
   dshape(13,2) = (9*x*(-1 + 6*z))/2.;
   dshape(15,0) = 0;
   dshape(15,1) = (9*z*(-1 + 3*z))/2.;
   dshape(15,2) = (9*y*(-1 + 6*z))/2.;
   dshape(3,0) = 0;
   dshape(3,1) = 0;
   dshape(3,2) = 1 + (9*z*(-2 + 3*z))/2.;
}

void Cubic3DFiniteElement::CalcHessian(const IntegrationPoint &ip,
                                       DenseMatrix &h) const
{
   double x = ip.x, y = ip.y, z = ip.z;


   // always in the order u_xx, u_xy, u_xz, u_yy, u_yz, u_zz

   h(0,0) = -9.*(-2. + 3.*x + 3.*y + 3.*z);
   h(0,1) = -9.*(-2. + 3.*x + 3.*y + 3.*z);
   h(0,2) = -9.*(-2. + 3.*x + 3.*y + 3.*z);
   h(0,3) = -9.*(-2. + 3.*x + 3.*y + 3.*z);
   h(0,4) = -9.*(-2. + 3.*x + 3.*y + 3.*z);
   h(0,5) = -9.*(-2. + 3.*x + 3.*y + 3.*z);

   h(1,0) = -9. + 27.*x;
   h(1,1) = 0.;
   h(1,2) = 0.;
   h(1,3) = 0.;
   h(1,4) = 0.;
   h(1,5) = 0.;

   h(2,0) = 0.;
   h(2,1) = 0.;
   h(2,2) = 0.;
   h(2,3) = -9. + 27.*y;
   h(2,4) = 0.;
   h(2,5) = 0.;

   h(3,0) = 0.;
   h(3,1) = 0.;
   h(3,2) = 0.;
   h(3,3) = 0.;
   h(3,4) = 0.;
   h(3,5) =-9. + 27.*z;

   h(4,0) = 9.*(-5. + 9.*x + 6.*y + 6.*z);
   h(4,1) = 4.5*(-5. + 12.*x + 6.*y + 6.*z);
   h(4,2) = 4.5*(-5. + 12.*x + 6.*y + 6.*z);
   h(4,3) = 27.*x;
   h(4,4) = 27.*x;
   h(4,5) = 27.*x;

   h(5,0) = -9.*(-4. + 9.*x + 3.*y + 3.*z);
   h(5,1) = 4.5 - 27.*x;
   h(5,2) = 4.5 - 27.*x;
   h(5,3) = 0.;
   h(5,4) = 0.;
   h(5,5) = 0.;

   h(6,0) = 27.*y;
   h(6,1) = 4.5*(-5. + 6.*x + 12.*y + 6.*z);
   h(6,2) = 27.*y;
   h(6,3) = 9.*(-5. + 6.*x + 9.*y + 6.*z);
   h(6,4) = 4.5*(-5. + 6.*x + 12.*y + 6.*z);
   h(6,5) = 27.*y;

   h(7,0)=0.;
   h(7,1)=4.5 - 27.*y;
   h(7,2)=0.;
   h(7,3)=-9.*(-4. + 3.*x + 9.*y + 3.*z);
   h(7,4)=4.5 - 27.*y;
   h(7,5)=0.;

   h(8,0)=27.*z;
   h(8,1)=27.*z;
   h(8,2)=4.5*(-5. + 6.*x + 6.*y + 12.*z);
   h(8,3)=27.*z;
   h(8,4)=4.5*(-5. + 6.*x + 6.*y + 12.*z);
   h(8,5)=9.*(-5. + 6.*x + 6.*y + 9.*z);

   h(9,0) = 0.;
   h(9,1) = 0.;
   h(9,2) = 4.5 - 27.*z;
   h(9,3) = 0.;
   h(9,4)=4.5 - 27.*z;
   h(9,5)=-9.*(-4. + 3.*x + 3.*y + 9.*z);

   h(10,0)=27.*y;
   h(10,1)=-4.5 + 27.*x;
   h(10,2)=0.;
   h(10,3)=0.;
   h(10,4)=0.;
   h(10,5)=0.;

   h(11,0)=0.;
   h(11,1)=-4.5 + 27.*y;
   h(11,2)=0.;
   h(11,3)=27.*x;
   h(11,4)=0.;
   h(11,5)=0.;

   h(12,0)=27.*z;
   h(12,1)=0.;
   h(12,2)=-4.5 + 27.*x;
   h(12,3)=0.;
   h(12,4)=0.;
   h(12,5)=0.;

   h(13,0)=0.;
   h(13,1)=0.;
   h(13,2)=-4.5 + 27.*z;
   h(13,3)=0.;
   h(13,4)=0.;
   h(13,5)=27.*x;

   h(14,0)=0.;
   h(14,1)=0.;
   h(14,2)=0.;
   h(14,3)=27.*z;
   h(14,4)=-4.5 + 27.*y;
   h(14,5)=0.;

   h(15,0)=0.;
   h(15,1)=0.;
   h(15,2)=0.;
   h(15,3)=0.;
   h(15,4)=-4.5 + 27.*z;
   h(15,5)=27.*y;

   h(16,0)=0.;
   h(16,1)=27.*z;
   h(16,2)=27.*y;
   h(16,3)=0.;
   h(16,4)=27.*x;
   h(16,5)=0.;

   h(17,0)=0.;
   h(17,1)=-27.*z;
   h(17,2)=-27.*y;
   h(17,3)=-54.*z;
   h(17,4)=-27.*(-1. + x + 2.*y + 2.*z);
   h(17,5)=-54.*y;

   h(18,0)=-54.*z;
   h(18,1)=-27.*z;
   h(18,2)=-27.*(-1. + 2.*x + y + 2.*z);
   h(18,3)=0.;
   h(18,4)=-27.*x;
   h(18,5)=-54.*x;

   h(19,0)=-54.*y;
   h(19,1)=-27.*(-1. + 2.*x + 2.*y + z);
   h(19,2)=-27.*y;
   h(19,3)=-54.*x;
   h(19,4)=-27.*x;
   h(19,5)=0.;
}


P0TriangleFiniteElement::P0TriangleFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 1, 0)
{
   Nodes.IntPoint(0).x = 0.333333333333333333;
   Nodes.IntPoint(0).y = 0.333333333333333333;
}

void P0TriangleFiniteElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   shape(0) = 1.0;
}

void P0TriangleFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   dshape(0,0) = 0.0;
   dshape(0,1) = 0.0;
}


P0QuadFiniteElement::P0QuadFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 1, 0, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.5;
}

void P0QuadFiniteElement::CalcShape(const IntegrationPoint &ip,
                                    Vector &shape) const
{
   shape(0) = 1.0;
}

void P0QuadFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   dshape(0,0) = 0.0;
   dshape(0,1) = 0.0;
}


Linear3DFiniteElement::Linear3DFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 4, 1)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 1.0;
}

void Linear3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   shape(0) = 1. - ip.x - ip.y - ip.z;
   shape(1) = ip.x;
   shape(2) = ip.y;
   shape(3) = ip.z;
}

void Linear3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   if (dshape.Height() == 4)
   {
      double *A = &dshape(0,0);
      A[0] = -1.; A[4] = -1.; A[8]  = -1.;
      A[1] =  1.; A[5] =  0.; A[9]  =  0.;
      A[2] =  0.; A[6] =  1.; A[10] =  0.;
      A[3] =  0.; A[7] =  0.; A[11] =  1.;
   }
   else
   {
      dshape(0,0) = -1.; dshape(0,1) = -1.; dshape(0,2) = -1.;
      dshape(1,0) =  1.; dshape(1,1) =  0.; dshape(1,2) =  0.;
      dshape(2,0) =  0.; dshape(2,1) =  1.; dshape(2,2) =  0.;
      dshape(3,0) =  0.; dshape(3,1) =  0.; dshape(3,2) =  1.;
   }
}

void Linear3DFiniteElement::GetFaceDofs (int face, int **dofs, int *ndofs)
const
{
   static int face_dofs[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

   *ndofs = 3;
   *dofs  = face_dofs[face];
}


Quadratic3DFiniteElement::Quadratic3DFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 10, 2)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 1.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 0.0;
   Nodes.IntPoint(6).x = 0.0;
   Nodes.IntPoint(6).y = 0.0;
   Nodes.IntPoint(6).z = 0.5;
   Nodes.IntPoint(7).x = 0.5;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(7).z = 0.0;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.5;
   Nodes.IntPoint(9).x = 0.0;
   Nodes.IntPoint(9).y = 0.5;
   Nodes.IntPoint(9).z = 0.5;
}

void Quadratic3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   double L0, L1, L2, L3;

   L0 = 1. - ip.x - ip.y - ip.z;
   L1 = ip.x;
   L2 = ip.y;
   L3 = ip.z;

   shape(0) = L0 * ( 2.0 * L0 - 1.0 );
   shape(1) = L1 * ( 2.0 * L1 - 1.0 );
   shape(2) = L2 * ( 2.0 * L2 - 1.0 );
   shape(3) = L3 * ( 2.0 * L3 - 1.0 );
   shape(4) = 4.0 * L0 * L1;
   shape(5) = 4.0 * L0 * L2;
   shape(6) = 4.0 * L0 * L3;
   shape(7) = 4.0 * L1 * L2;
   shape(8) = 4.0 * L1 * L3;
   shape(9) = 4.0 * L2 * L3;
}

void Quadratic3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   double x, y, z, L0;

   x = ip.x;
   y = ip.y;
   z = ip.z;
   L0 = 1.0 - x - y - z;

   dshape(0,0) = dshape(0,1) = dshape(0,2) = 1.0 - 4.0 * L0;
   dshape(1,0) = -1.0 + 4.0 * x; dshape(1,1) = 0.0; dshape(1,2) = 0.0;
   dshape(2,0) = 0.0; dshape(2,1) = -1.0 + 4.0 * y; dshape(2,2) = 0.0;
   dshape(3,0) = dshape(3,1) = 0.0; dshape(3,2) = -1.0 + 4.0 * z;
   dshape(4,0) = 4.0 * (L0 - x); dshape(4,1) = dshape(4,2) = -4.0 * x;
   dshape(5,0) = dshape(5,2) = -4.0 * y; dshape(5,1) = 4.0 * (L0 - y);
   dshape(6,0) = dshape(6,1) = -4.0 * z; dshape(6,2) = 4.0 * (L0 - z);
   dshape(7,0) = 4.0 * y; dshape(7,1) = 4.0 * x; dshape(7,2) = 0.0;
   dshape(8,0) = 4.0 * z; dshape(8,1) = 0.0; dshape(8,2) = 4.0 * x;
   dshape(9,0) = 0.0; dshape(9,1) = 4.0 * z; dshape(9,2) = 4.0 * y;
}

void Quadratic3DFiniteElement::CalcHessian(const IntegrationPoint &ip,
                                           DenseMatrix &h) const
{
   //u_xx
   h(0,0) = 4.;
   h(1,0) = 4.;
   h(2,0) = 0.;
   h(3,0) = 0.;
   h(4,0) = -8.;
   h(5,0) = 0.;
   h(6,0) = 0.;
   h(7,0) = 0.;
   h(8,0) = 0.;
   h(9,0) = 0.;

   //u_xy
   h(0,1) = 4.;
   h(1,1) = 0.;
   h(2,1) = 0.;
   h(3,1) = 0.;
   h(4,1) = -4.;
   h(5,1) = -4.;
   h(6,1) = 0.;
   h(7,1) = 4.;
   h(8,1) = 0.;
   h(9,1) = 0.;

   //u_xz
   h(0,2) = 4.;
   h(1,2) = 0.;
   h(2,2) = 0.;
   h(3,2) = 0.;
   h(4,2) = -4.;
   h(5,2) = 0.;
   h(6,2) = -4.;
   h(7,2) = 0.;
   h(8,2) = 4.;
   h(9,2) = 0.;

   //u_yy
   h(0,3) = 4.;
   h(1,3) = 0.;
   h(2,3) = 4.;
   h(3,3) = 0.;
   h(4,3) = 0.;
   h(5,3) = -8.;
   h(6,3) = 0.;
   h(7,3) = 0.;
   h(8,3) = 0.;
   h(9,3) = 0.;

   //u_yz
   h(0,4) = 4.;
   h(1,4) = 0.;
   h(2,4) = 0.;
   h(3,4) = 0.;
   h(4,4) = 0.;
   h(5,4) = -4.;
   h(6,4) = -4.;
   h(7,4) = 0.;
   h(8,4) = 0.;
   h(9,4) = 4.;

   //u_zz
   h(0,5) = 4.;
   h(1,5) = 0.;
   h(2,5) = 0.;
   h(3,5) = 4.;
   h(4,5) = 0.;
   h(5,5) = 0.;
   h(6,5) = -8.;
   h(7,5) = 0.;
   h(8,5) = 0.;
   h(9,5) = 0.;

}

TriLinear3DFiniteElement::TriLinear3DFiniteElement()
   : NodalFiniteElement(3, Geometry::CUBE, 8, 1, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;

   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(3).z = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 1.0;

   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 1.0;

   Nodes.IntPoint(6).x = 1.0;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(6).z = 1.0;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 1.0;
   Nodes.IntPoint(7).z = 1.0;
}

void TriLinear3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;
   double ox = 1.-x, oy = 1.-y, oz = 1.-z;

   shape(0) = ox * oy * oz;
   shape(1) =  x * oy * oz;
   shape(2) =  x *  y * oz;
   shape(3) = ox *  y * oz;
   shape(4) = ox * oy *  z;
   shape(5) =  x * oy *  z;
   shape(6) =  x *  y *  z;
   shape(7) = ox *  y *  z;
}

void TriLinear3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y, z = ip.z;
   double ox = 1.-x, oy = 1.-y, oz = 1.-z;

   dshape(0,0) = - oy * oz;
   dshape(0,1) = - ox * oz;
   dshape(0,2) = - ox * oy;

   dshape(1,0) =   oy * oz;
   dshape(1,1) = -  x * oz;
   dshape(1,2) = -  x * oy;

   dshape(2,0) =    y * oz;
   dshape(2,1) =    x * oz;
   dshape(2,2) = -  x *  y;

   dshape(3,0) = -  y * oz;
   dshape(3,1) =   ox * oz;
   dshape(3,2) = - ox *  y;

   dshape(4,0) = - oy *  z;
   dshape(4,1) = - ox *  z;
   dshape(4,2) =   ox * oy;

   dshape(5,0) =   oy *  z;
   dshape(5,1) = -  x *  z;
   dshape(5,2) =    x * oy;

   dshape(6,0) =    y *  z;
   dshape(6,1) =    x *  z;
   dshape(6,2) =    x *  y;

   dshape(7,0) = -  y *  z;
   dshape(7,1) =   ox *  z;
   dshape(7,2) =   ox *  y;
}

QuadLinear4DFiniteElement::QuadLinear4DFiniteElement()
   : NodalFiniteElement(4, Geometry::TESSERACT, 16, 1, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(2).t = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(3).z = 0.0;
   Nodes.IntPoint(3).t = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 1.0;
   Nodes.IntPoint(4).t = 0.0;

   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 1.0;
   Nodes.IntPoint(5).t = 0.0;

   Nodes.IntPoint(6).x = 1.0;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(6).z = 1.0;
   Nodes.IntPoint(6).t = 0.0;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 1.0;
   Nodes.IntPoint(7).z = 1.0;
   Nodes.IntPoint(7).t = 0.0;


   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.0;
   Nodes.IntPoint(8).t = 1.0;

   Nodes.IntPoint(9).x = 1.0;
   Nodes.IntPoint(9).y = 0.0;
   Nodes.IntPoint(9).z = 0.0;
   Nodes.IntPoint(9).t = 1.0;

   Nodes.IntPoint(10).x = 1.0;
   Nodes.IntPoint(10).y = 1.0;
   Nodes.IntPoint(10).z = 0.0;
   Nodes.IntPoint(10).t = 1.0;

   Nodes.IntPoint(11).x = 0.0;
   Nodes.IntPoint(11).y = 1.0;
   Nodes.IntPoint(11).z = 0.0;
   Nodes.IntPoint(11).t = 1.0;

   Nodes.IntPoint(12).x = 0.0;
   Nodes.IntPoint(12).y = 0.0;
   Nodes.IntPoint(12).z = 1.0;
   Nodes.IntPoint(12).t = 1.0;

   Nodes.IntPoint(13).x = 1.0;
   Nodes.IntPoint(13).y = 0.0;
   Nodes.IntPoint(13).z = 1.0;
   Nodes.IntPoint(13).t = 1.0;

   Nodes.IntPoint(14).x = 1.0;
   Nodes.IntPoint(14).y = 1.0;
   Nodes.IntPoint(14).z = 1.0;
   Nodes.IntPoint(14).t = 1.0;

   Nodes.IntPoint(15).x = 0.0;
   Nodes.IntPoint(15).y = 1.0;
   Nodes.IntPoint(15).z = 1.0;
   Nodes.IntPoint(15).t = 1.0;
}

void QuadLinear4DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                          Vector &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z, t = ip.t;
   double ox = 1.-x, oy = 1.-y, oz = 1.-z, ot = 1.-t;

   shape(0) = ox * oy * oz * ot;
   shape(1) =  x * oy * oz * ot;
   shape(2) =  x *  y * oz * ot;
   shape(3) = ox *  y * oz * ot;
   shape(4) = ox * oy *  z * ot;
   shape(5) =  x * oy *  z * ot;
   shape(6) =  x *  y *  z * ot;
   shape(7) = ox *  y *  z * ot;

   shape(8) = ox * oy * oz * t;
   shape(9) =  x * oy * oz * t;
   shape(10) =  x *  y * oz * t;
   shape(11) = ox *  y * oz * t;
   shape(12) = ox * oy *  z * t;
   shape(13) =  x * oy *  z * t;
   shape(14) =  x *  y *  z * t;
   shape(15) = ox *  y *  z * t;
}

void QuadLinear4DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                           DenseMatrix &dshape) const
{
   double x = ip.x, y = ip.y, z = ip.z, t = ip.t;
   double ox = 1.-x, oy = 1.-y, oz = 1.-z, ot = 1.-t;

   dshape(0,0) = - oy * oz * ot;
   dshape(0,1) = - ox * oz * ot;
   dshape(0,2) = - ox * oy * ot;
   dshape(0,3) = - ox * oy * oz;

   dshape(1,0) =   oy * oz * ot;
   dshape(1,1) = -  x * oz * ot;
   dshape(1,2) = -  x * oy * ot;
   dshape(1,3) = -  x * oy * oz;

   dshape(2,0) =    y * oz * ot;
   dshape(2,1) =    x * oz * ot;
   dshape(2,2) = -  x *  y * ot;
   dshape(2,3) = -  x *  y * oz;

   dshape(3,0) = -  y * oz * ot;
   dshape(3,1) =   ox * oz * ot;
   dshape(3,2) = - ox *  y * ot;
   dshape(3,3) = - ox *  y * oz;

   dshape(4,0) = - oy *  z * ot;
   dshape(4,1) = - ox *  z * ot;
   dshape(4,2) =   ox * oy * ot;
   dshape(4,3) = - ox * oy *  z;

   dshape(5,0) =   oy *  z * ot;
   dshape(5,1) = -  x *  z * ot;
   dshape(5,2) =    x * oy * ot;
   dshape(5,3) = -  x * oy *  z;

   dshape(6,0) =    y *  z * ot;
   dshape(6,1) =    x *  z * ot;
   dshape(6,2) =    x *  y * ot;
   dshape(6,3) = -  x *  y *  z;

   dshape(7,0) = -  y *  z * ot;
   dshape(7,1) =   ox *  z * ot;
   dshape(7,2) =   ox *  y * ot;
   dshape(7,3) = - ox *  y *  z;


   dshape(8,0) = - oy * oz * t;
   dshape(8,1) = - ox * oz * t;
   dshape(8,2) = - ox * oy * t;
   dshape(8,3) =   ox * oy * oz;

   dshape(9,0) =   oy * oz * t;
   dshape(9,1) = -  x * oz * t;
   dshape(9,2) = -  x * oy * t;
   dshape(9,3) =    x * oy * oz;

   dshape(10,0) =    y * oz * t;
   dshape(10,1) =    x * oz * t;
   dshape(10,2) = -  x *  y * t;
   dshape(10,3) =    x *  y * oz;

   dshape(11,0) = -  y * oz * t;
   dshape(11,1) =   ox * oz * t;
   dshape(11,2) = - ox *  y * t;
   dshape(11,3) =   ox *  y * oz;

   dshape(12,0) = - oy *  z * t;
   dshape(12,1) = - ox *  z * t;
   dshape(12,2) =   ox * oy * t;
   dshape(12,3) =   ox * oy *  z;

   dshape(13,0) =   oy *  z * t;
   dshape(13,1) = -  x *  z * t;
   dshape(13,2) =    x * oy * t;
   dshape(13,3) =    x * oy *  z;

   dshape(14,0) =    y *  z * t;
   dshape(14,1) =    x *  z * t;
   dshape(14,2) =    x *  y * t;
   dshape(14,3) =    x *  y *  z;

   dshape(15,0) = -  y *  z * t;
   dshape(15,1) =   ox *  z * t;
   dshape(15,2) =   ox *  y * t;
   dshape(15,3) =   ox *  y *  z;
}

void QuadLinear4DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                             DenseMatrix &h) const
{
   double x = ip.x, y = ip.y, z = ip.z, t = ip.t;
   double ox = 1.-x, oy = 1.-y, oz = 1.-z, ot = 1.-t;


   //u_xx
   h(0,0)  = 0.0;
   h(1,0)  = 0.0;
   h(2,0)  = 0.0;
   h(3,0)  = 0.0;
   h(4,0)  = 0.0;
   h(5,0)  = 0.0;
   h(6,0)  = 0.0;
   h(7,0)  = 0.0;
   h(8,0)  = 0.0;
   h(9,0)  = 0.0;
   h(10,0) = 0.0;
   h(11,0) = 0.0;
   h(12,0) = 0.0;
   h(13,0) = 0.0;
   h(14,0) = 0.0;
   h(15,0) = 0.0;

   //u_xy
   h(0,1)  =  ot * oz;
   h(1,1)  = -ot * oz;
   h(2,1)  =  ot * oz;
   h(3,1)  = -ot * oz;
   h(4,1)  =  ot * z;
   h(5,1)  = -ot * z;
   h(6,1)  =  ot * z;
   h(7,1)  = -ot * z;
   h(8,1)  =  t * oz;
   h(9,1)  = -t * oz;
   h(10,1) =  t * oz;
   h(11,1) = -t * oz;
   h(12,1) =  t * z;
   h(13,1) = -t * z;
   h(14,1) =  t * z;
   h(15,1) = -t * z;

   //u_xz
   h(0,2)  =  ot * oy;
   h(1,2)  = -ot * oy;
   h(2,2)  = -ot *  y;
   h(3,2)  =  ot *  y;
   h(4,2)  = -ot * oy;
   h(5,2)  =  ot * oy;
   h(6,2)  =  ot *  y;
   h(7,2)  = -ot *  y;
   h(8,2)  =   t * oy;
   h(9,2)  = - t * oy;
   h(10,2) = - t *  y;
   h(11,2) =   t *  y;
   h(12,2) = - t * oy;
   h(13,2) =   t * oy;
   h(14,2) =   t *  y;
   h(15,2) = - t *  y;

   //u_xt
   h(0,3)  =  oy * oz;
   h(1,3)  = -oy * oz;
   h(2,3)  = - y * oz;
   h(3,3)  =   y * oz;
   h(4,3)  =  oy *  z;
   h(5,3)  = -oy *  z;
   h(6,3)  = - y *  z;
   h(7,3)  =   y *  z;
   h(8,3)  = -oy * oz;
   h(9,3)  =  oy * oz;
   h(10,3) =   y * oz;
   h(11,3) = - y * oz;
   h(12,3) = -oy *  z;
   h(13,3) =  oy *  z;
   h(14,3) =   y *  z;
   h(15,3) = - y *  z;

   //u_yy
   h(0,4)  = 0.0;
   h(1,4)  = 0.0;
   h(2,4)  = 0.0;
   h(3,4)  = 0.0;
   h(4,4)  = 0.0;
   h(5,4)  = 0.0;
   h(6,4)  = 0.0;
   h(7,4)  = 0.0;
   h(8,4)  = 0.0;
   h(9,4)  = 0.0;
   h(10,4) = 0.0;
   h(11,4) = 0.0;
   h(12,4) = 0.0;
   h(13,4) = 0.0;
   h(14,4) = 0.0;
   h(15,4) = 0.0;

   //u_yz
   h(0,5)  =  ot * ox;
   h(1,5)  =  ot *  x;
   h(2,5)  = -ot *  x;
   h(3,5)  = -ot * ox;
   h(4,5)  = -ot * ox;
   h(5,5)  = -ot *  x;
   h(6,5)  =  ot *  x;
   h(7,5)  =  ot * ox;
   h(8,5)  =   t * ox;
   h(9,5)  =   t *  x;
   h(10,5) = - t *  x;
   h(11,5) = - t * ox;
   h(12,5) = - t * ox;
   h(13,5) = - t *  x;
   h(14,5) =   t *  x;
   h(15,5) =   t * ox;

   //u_yt
   h(0,6)  =  ox * oz;
   h(1,6)  =   x * oz;
   h(2,6)  = - x * oz;
   h(3,6)  = -ox * oz;
   h(4,6)  =  ox *  z;
   h(5,6)  =   x *  z;
   h(6,6)  = - x *  z;
   h(7,6)  = -ox *  z;
   h(8,6)  = -ox * oz;
   h(9,6)  = - x * oz;
   h(10,6) =   x * oz;
   h(11,6) =  ox * oz;
   h(12,6) = -ox *  z;
   h(13,6) = - x *  z;
   h(14,6) =   x *  z;
   h(15,6) =  ox *  z;

   //u_zz
   h(0,7)  = 0.0;
   h(1,7)  = 0.0;
   h(2,7)  = 0.0;
   h(3,7)  = 0.0;
   h(4,7)  = 0.0;
   h(5,7)  = 0.0;
   h(6,7)  = 0.0;
   h(7,7)  = 0.0;
   h(8,7)  = 0.0;
   h(9,7)  = 0.0;
   h(10,7) = 0.0;
   h(11,7) = 0.0;
   h(12,7) = 0.0;
   h(13,7) = 0.0;
   h(14,7) = 0.0;
   h(15,7) = 0.0;

   //u_zt
   h(0,8)  =  ox * oy;
   h(1,8)  =   x * oy;
   h(2,8)  =   x *  y;
   h(3,8)  =  ox *  y;
   h(4,8)  = -ox * oy;
   h(5,8)  = - x * oy;
   h(6,8)  = - x *  y;
   h(7,8)  = -ox *  y;
   h(8,8)  = -ox * oy;
   h(9,8)  = - x * oy;
   h(10,8) = - x *  y;
   h(11,8) = -ox *  y;
   h(12,8) =  ox * oy;
   h(13,8) =   x * oy;
   h(14,8) =   x *  y;
   h(15,8) =  ox *  y;

   //u_tt
   h(0,9)  = 0.0;
   h(1,9)  = 0.0;
   h(2,9)  = 0.0;
   h(3,9)  = 0.0;
   h(4,9)  = 0.0;
   h(5,9)  = 0.0;
   h(6,9)  = 0.0;
   h(7,9)  = 0.0;
   h(8,9)  = 0.0;
   h(9,9)  = 0.0;
   h(10,9) = 0.0;
   h(11,9) = 0.0;
   h(12,9) = 0.0;
   h(13,9) = 0.0;
   h(14,9) = 0.0;
   h(15,9) = 0.0;

}

Linear4DFiniteElement::Linear4DFiniteElement()
   : NodalFiniteElement(4, Geometry::PENTATOPE, 5, 1)
{

   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(2).t = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 1.0;
   Nodes.IntPoint(3).t = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(4).t = 1.0;
}

void Linear4DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   shape(0) = 1. - ip.x - ip.y - ip.z - ip.t;
   shape(1) = ip.x;
   shape(2) = ip.y;
   shape(3) = ip.z;
   shape(4) = ip.t;
}

void Linear4DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   if (dshape.Height() == 5)
   {
      double *A = &dshape(0,0);
      A[0] = -1.; A[5] = -1.; A[10] = -1.; A[15] = -1.;
      A[1] =  1.; A[6] =  0.; A[11] =  0.; A[16] =  0.;
      A[2] =  0.; A[7] =  1.; A[12] =  0.; A[17] =  0.;
      A[3] =  0.; A[8] =  0.; A[13] =  1.; A[18] =  0.;
      A[4] =  0.; A[9] =  0.; A[14] =  0.; A[19] =  1.;
   }
   else
   {
      dshape(0,0) = -1.; dshape(0,1) = -1.; dshape(0,2) = -1.; dshape(0,3) = -1.;
      dshape(1,0) =  1.; dshape(1,1) =  0.; dshape(1,2) =  0.; dshape(1,3) =  0.;
      dshape(2,0) =  0.; dshape(2,1) =  1.; dshape(2,2) =  0.; dshape(2,3) =  0.;
      dshape(3,0) =  0.; dshape(3,1) =  0.; dshape(3,2) =  1.; dshape(3,3) =  0.;
      dshape(4,0) =  0.; dshape(4,1) =  0.; dshape(4,2) =  0.; dshape(4,3) =  1.;
   }
}

void Linear4DFiniteElement::CalcHessian(const IntegrationPoint &ip,
                                        DenseMatrix &h) const
{
   h = 0.0;
}

void Linear4DFiniteElement::GetFaceDofs (int face, int **dofs, int *ndofs)
const
{
   static int face_dofs[5][4] = {{0, 1, 2, 3}, {1, 0, 2, 4}, {0, 1, 3, 4}, {2, 0, 3, 4}, {1, 2, 3, 4}};

   *ndofs = 4;
   *dofs  = face_dofs[face];
}



Quadratic4DFiniteElement::Quadratic4DFiniteElement()
   : NodalFiniteElement(4, Geometry::PENTATOPE, 15, 2)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(2).t = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 1.0;
   Nodes.IntPoint(3).t = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(4).t = 1.0;

   Nodes.IntPoint(5).x = 0.5;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 0.0;
   Nodes.IntPoint(5).t = 0.0;

   Nodes.IntPoint(6).x = 0.0;
   Nodes.IntPoint(6).y = 0.5;
   Nodes.IntPoint(6).z = 0.0;
   Nodes.IntPoint(6).t = 0.0;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.0;
   Nodes.IntPoint(7).z = 0.5;
   Nodes.IntPoint(7).t = 0.0;

   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.0;
   Nodes.IntPoint(8).t = 0.5;

   Nodes.IntPoint(9).x = 0.5;
   Nodes.IntPoint(9).y = 0.5;
   Nodes.IntPoint(9).z = 0.0;
   Nodes.IntPoint(9).t = 0.0;

   Nodes.IntPoint(10).x = 0.5;
   Nodes.IntPoint(10).y = 0.0;
   Nodes.IntPoint(10).z = 0.5;
   Nodes.IntPoint(10).t = 0.0;

   Nodes.IntPoint(11).x = 0.5;
   Nodes.IntPoint(11).y = 0.0;
   Nodes.IntPoint(11).z = 0.0;
   Nodes.IntPoint(11).t = 0.5;

   Nodes.IntPoint(12).x = 0.0;
   Nodes.IntPoint(12).y = 0.5;
   Nodes.IntPoint(12).z = 0.5;
   Nodes.IntPoint(12).t = 0.0;

   Nodes.IntPoint(13).x = 0.0;
   Nodes.IntPoint(13).y = 0.5;
   Nodes.IntPoint(13).z = 0.0;
   Nodes.IntPoint(13).t = 0.5;

   Nodes.IntPoint(14).x = 0.0;
   Nodes.IntPoint(14).y = 0.0;
   Nodes.IntPoint(14).z = 0.5;
   Nodes.IntPoint(14).t = 0.5;
}

void Quadratic4DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   double L0, L1, L2, L3, L4;

   L0 = 1. - ip.x - ip.y - ip.z - ip.t;
   L1 = ip.x;
   L2 = ip.y;
   L3 = ip.z;
   L4 = ip.t;

   shape(0) = L0 * ( 2.0 * L0 - 1.0 );
   shape(1) = L1 * ( 2.0 * L1 - 1.0 );
   shape(2) = L2 * ( 2.0 * L2 - 1.0 );
   shape(3) = L3 * ( 2.0 * L3 - 1.0 );
   shape(4) = L4 * ( 2.0 * L4 - 1.0 );
   shape(5) = 4.0 * L0 * L1;
   shape(6) = 4.0 * L0 * L2;
   shape(7) = 4.0 * L0 * L3;
   shape(8) = 4.0 * L0 * L4;
   shape(9) = 4.0 * L1 * L2;
   shape(10) = 4.0 * L1 * L3;
   shape(11) = 4.0 * L1 * L4;
   shape(12) = 4.0 * L2 * L3;
   shape(13) = 4.0 * L2 * L4;
   shape(14) = 4.0 * L3 * L4;
}

void Quadratic4DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   double x, y, z, t, L0;

   x = ip.x;
   y = ip.y;
   z = ip.z;
   t = ip.t;
   L0 = 1.0 - x - y - z - t;

   dshape(0,0) = dshape(0,1) = dshape(0,2) = dshape(0,3) = 1.0 - 4.0 * L0;

   dshape(1,0) = -1.0 + 4.0 * x; dshape(1,1) = 0.0; dshape(1,2) = 0.0;
   dshape(1,3) = 0.0;
   dshape(2,0) = 0.0; dshape(2,1) = -1.0 + 4.0 * y; dshape(2,2) = 0.0;
   dshape(2,3) = 0.0;
   dshape(3,0) = 0.0; dshape(3,1) = 0.0; dshape(3,2) = -1.0 + 4.0 * z;
   dshape(3,3) = 0.0;
   dshape(4,0) = 0.0; dshape(4,1) = 0.0; dshape(4,2) = 0.0;
   dshape(4,3) = -1.0 + 4.0 * t;

   dshape(5,0) = 4.0 * (L0 - x);
   dshape(5,1) = dshape(5,2) = dshape(5,3) = -4.0 * x;
   dshape(6,1) = 4.0 * (L0 - y);
   dshape(6,0) = dshape(6,2) = dshape(6,3) = -4.0 * y;
   dshape(7,2) = 4.0 * (L0 - z);
   dshape(7,0) = dshape(7,1) = dshape(7,3) = -4.0 * z;
   dshape(8,3) = 4.0 * (L0 - t);
   dshape(8,0) = dshape(8,1) = dshape(8,2) = -4.0 * t;

   dshape(9,0) = 4.0 * y; dshape(9,1) = 4.0 * x; dshape(9,2) = 0.0;
   dshape(9,3) = 0.0;
   dshape(10,0) = 4.0 * z; dshape(10,1) = 0.0; dshape(10,2) = 4.0 * x;
   dshape(10,3) = 0.0;
   dshape(11,0) = 4.0 * t; dshape(11,1) = 0.0; dshape(11,2) = 0.0;
   dshape(11,3) = 4.0 * x;

   dshape(12,0) = 0.0; dshape(12,1) = 4.0 * z; dshape(12,2) = 4.0 * y;
   dshape(12,3) = 0.0;
   dshape(13,0) = 0.0; dshape(13,1) = 4.0 * t; dshape(13,2) = 0.0;
   dshape(13,3) = 4.0 * y;
   dshape(14,0) = 0.0; dshape(14,1) = 0.0; dshape(14,2) = 4.0 * t;
   dshape(14,3) = 4.0 * z;
}

void Quadratic4DFiniteElement::CalcHessian(const IntegrationPoint &ip,
                                           DenseMatrix &h) const
{
   //u_xx
   h(0,0) = 4.;
   h(1,0) = 4.;
   h(2,0) = 0.;
   h(3,0) = 0.;
   h(4,0) = 0.;
   h(5,0) = -8.;
   h(6,0) = 0.;
   h(7,0) = 0.;
   h(8,0) = 0.;
   h(9,0) = 0.;
   h(10,0) = 0.;
   h(11,0) = 0.;
   h(12,0) = 0.;
   h(13,0) = 0.;
   h(14,0) = 0.;

   //u_xy
   h(0,1) = 4.;
   h(1,1) = 0.;
   h(2,1) = 0.;
   h(3,1) = 0.;
   h(4,1) = 0.;
   h(5,1) = -4.;
   h(6,1) = -4.;
   h(7,1) = 0.;
   h(8,1) = 0.;
   h(9,1) = 4.;
   h(10,1) = 0.;
   h(11,1) = 0.;
   h(12,1) = 0.;
   h(13,1) = 0.;
   h(14,1) = 0.;

   //u_xz
   h(0,2) = 4.;
   h(1,2) = 0.;
   h(2,2) = 0.;
   h(3,2) = 0.;
   h(4,2) = 0.;
   h(5,2) = -4.;
   h(6,2) = 0.;
   h(7,2) = -4.;
   h(8,2) = 0.;
   h(9,2) = 0.;
   h(10,2) = 4.;
   h(11,2) = 0.;
   h(12,2) = 0.;
   h(13,2) = 0.;
   h(14,2) = 0.;

   //u_xt
   h(0,3) = 4.;
   h(1,3) = 0.;
   h(2,3) = 0.;
   h(3,3) = 0.;
   h(4,3) = 0.;
   h(5,3) = -4.;
   h(6,3) = 0.;
   h(7,3) = 0.;
   h(8,3) = -4.;
   h(9,3) = 0.;
   h(10,3) = 0.;
   h(11,3) = 4.;
   h(12,3) = 0.;
   h(13,3) = 0.;
   h(14,3) = 0.;

   //u_yy
   h(0,4) = 4.;
   h(1,4) = 0.;
   h(2,4) = 4.;
   h(3,4) = 0.;
   h(4,4) = 0.;
   h(5,4) = 0.;
   h(6,4) = -8.;
   h(7,4) = 0.;
   h(8,4) = 0.;
   h(9,4) = 0.;
   h(10,4) = 0.;
   h(11,4) = 0.;
   h(12,4) = 0.;
   h(13,4) = 0.;
   h(14,4) = 0.;

   //u_yz
   h(0,5) = 4.;
   h(1,5) = 0.;
   h(2,5) = 0.;
   h(3,5) = 0.;
   h(4,5) = 0.;
   h(5,5) = 0.;
   h(6,5) = -4.;
   h(7,5) = -4.;
   h(8,5) = 0.;
   h(9,5) = 0.;
   h(10,5) = 0.;
   h(11,5) = 0.;
   h(12,5) = 4.;
   h(13,5) = 0.;
   h(14,5) = 0.;

   //u_yt
   h(0,6) = 4.;
   h(1,6) = 0.;
   h(2,6) = 0.;
   h(3,6) = 0.;
   h(4,6) = 0.;
   h(5,6) = 0.;
   h(6,6) = -4.;
   h(7,6) = 0.;
   h(8,6) = -4.;
   h(9,6) = 0.;
   h(10,6) = 0.;
   h(11,6) = 0.;
   h(12,6) = 0.;
   h(13,6) = 4.;
   h(14,6) = 0.;

   //u_zz
   h(0,7) = 4.;
   h(1,7) = 0.;
   h(2,7) = 0.;
   h(3,7) = 4.;
   h(4,7) = 0.;
   h(5,7) = 0.;
   h(6,7) = 0.;
   h(7,7) = -8.;
   h(8,7) = 0.;
   h(9,7) = 0.;
   h(10,7) = 0.;
   h(11,7) = 0.;
   h(12,7) = 0.;
   h(13,7) = 0.;
   h(14,7) = 0.;

   //u_zt
   h(0,8) = 4.;
   h(1,8) = 0.;
   h(2,8) = 0.;
   h(3,8) = 0.;
   h(4,8) = 0.;
   h(5,8) = 0.;
   h(6,8) = 0.;
   h(7,8) = -4.;
   h(8,8) = -4.;
   h(9,8) = 0.;
   h(10,8) = 0.;
   h(11,8) = 0.;
   h(12,8) = 0.;
   h(13,8) = 0.;
   h(14,8) = 4.;

   //u_tt
   h(0,9) = 4.;
   h(1,9) = 0.;
   h(2,9) = 0.;
   h(3,9) = 0.;
   h(4,9) = 4.;
   h(5,9) = 0.;
   h(6,9) = 0.;
   h(7,9) = 0.;
   h(8,9) = -8.;
   h(9,9) = 0.;
   h(10,9) = 0.;
   h(11,9) = 0.;
   h(12,9) = 0.;
   h(13,9) = 0.;
   h(14,9) = 0.;
}

P0SegmentFiniteElement::P0SegmentFiniteElement(int Ord)
   : NodalFiniteElement(1, Geometry::SEGMENT, 1, Ord)   // default Ord = 0
{
   Nodes.IntPoint(0).x = 0.5;
}

void P0SegmentFiniteElement::CalcShape(const IntegrationPoint &ip,
                                       Vector &shape) const
{
   shape(0) = 1.0;
}

void P0SegmentFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                        DenseMatrix &dshape) const
{
   dshape(0,0) = 0.0;
}

CrouzeixRaviartFiniteElement::CrouzeixRaviartFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 3, 1)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 0.5;
}

void CrouzeixRaviartFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   shape(0) =  1.0 - 2.0 * ip.y;
   shape(1) = -1.0 + 2.0 * ( ip.x + ip.y );
   shape(2) =  1.0 - 2.0 * ip.x;
}

void CrouzeixRaviartFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   dshape(0,0) =  0.0; dshape(0,1) = -2.0;
   dshape(1,0) =  2.0; dshape(1,1) =  2.0;
   dshape(2,0) = -2.0; dshape(2,1) =  0.0;
}

CrouzeixRaviartQuadFiniteElement::CrouzeixRaviartQuadFiniteElement()
// the FunctionSpace should be rotated (45 degrees) Q_1
// i.e. the span of { 1, x, y, x^2 - y^2 }
   : NodalFiniteElement(2, Geometry::SQUARE, 4, 2, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(2).x = 0.5;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.5;
}

void CrouzeixRaviartQuadFiniteElement::CalcShape(const IntegrationPoint &ip,
                                                 Vector &shape) const
{
   const double l1 = ip.x+ip.y-0.5, l2 = 1.-l1, l3 = ip.x-ip.y+0.5, l4 = 1.-l3;

   shape(0) = l2 * l3;
   shape(1) = l1 * l3;
   shape(2) = l1 * l4;
   shape(3) = l2 * l4;
}

void CrouzeixRaviartQuadFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                                  DenseMatrix &dshape) const
{
   const double x2 = 2.*ip.x, y2 = 2.*ip.y;

   dshape(0,0) =  1. - x2; dshape(0,1) = -2. + y2;
   dshape(1,0) =       x2; dshape(1,1) =  1. - y2;
   dshape(2,0) =  1. - x2; dshape(2,1) =       y2;
   dshape(3,0) = -2. + x2; dshape(3,1) =  1. - y2;
}


RT0TriangleFiniteElement::RT0TriangleFiniteElement()
   : VectorFiniteElement(2, Geometry::TRIANGLE, 3, 1, H_DIV)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 0.5;
}

void RT0TriangleFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                          DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   shape(0,0) = x;
   shape(0,1) = y - 1.;
   shape(1,0) = x;
   shape(1,1) = y;
   shape(2,0) = x - 1.;
   shape(2,1) = y;
}

void RT0TriangleFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                            Vector &divshape) const
{
   divshape(0) = 2.;
   divshape(1) = 2.;
   divshape(2) = 2.;
}

const double RT0TriangleFiniteElement::nk[3][2] =
{ {0, -1}, {1, 1}, {-1, 0} };

void RT0TriangleFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 3; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 3; j++)
      {
         double d = vshape(j,0)*nk[k][0]+vshape(j,1)*nk[k][1];
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT0TriangleFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[2];
   Vector xk (vk, 2);

   for (k = 0; k < 3; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1];
      for (j = 0; j < 3; j++)
         if (fabs (I(k,j) = vshape(j,0)*vk[0]+vshape(j,1)*vk[1]) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT0TriangleFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[2];
   Vector xk (vk, 2);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 3; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) = (vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1] ) +
                 vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1] ));
   }
}

RT0QuadFiniteElement::RT0QuadFiniteElement()
   : VectorFiniteElement(2, Geometry::SQUARE, 4, 1, H_DIV, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(2).x = 0.5;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.5;
}

void RT0QuadFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   shape(0,0) = 0;
   shape(0,1) = y - 1.;
   shape(1,0) = x;
   shape(1,1) = 0;
   shape(2,0) = 0;
   shape(2,1) = y;
   shape(3,0) = x - 1.;
   shape(3,1) = 0;
}

void RT0QuadFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   divshape(0) = 1.;
   divshape(1) = 1.;
   divshape(2) = 1.;
   divshape(3) = 1.;
}

const double RT0QuadFiniteElement::nk[4][2] =
{ {0, -1}, {1, 0}, {0, 1}, {-1, 0} };

void RT0QuadFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 4; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 4; j++)
      {
         double d = vshape(j,0)*nk[k][0]+vshape(j,1)*nk[k][1];
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT0QuadFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear (more to have embedding?)
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[2];
   Vector xk (vk, 2);

   for (k = 0; k < 4; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1];
      for (j = 0; j < 4; j++)
         if (fabs (I(k,j) = vshape(j,0)*vk[0]+vshape(j,1)*vk[1]) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT0QuadFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[2];
   Vector xk (vk, 2);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 4; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) = (vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1] ) +
                 vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1] ));
   }
}

RT1TriangleFiniteElement::RT1TriangleFiniteElement()
   : VectorFiniteElement(2, Geometry::TRIANGLE, 8, 2, H_DIV)
{
   Nodes.IntPoint(0).x = 0.33333333333333333333;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 0.66666666666666666667;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.66666666666666666667;
   Nodes.IntPoint(2).y = 0.33333333333333333333;
   Nodes.IntPoint(3).x = 0.33333333333333333333;
   Nodes.IntPoint(3).y = 0.66666666666666666667;
   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.66666666666666666667;
   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.33333333333333333333;
   Nodes.IntPoint(6).x = 0.33333333333333333333;
   Nodes.IntPoint(6).y = 0.33333333333333333333;
   Nodes.IntPoint(7).x = 0.33333333333333333333;
   Nodes.IntPoint(7).y = 0.33333333333333333333;
}

void RT1TriangleFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                          DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   shape(0,0) = -2 * x * (-1 + x + 2 * y);
   shape(0,1) = -2 * (-1 + y) * (-1 + x + 2 * y);
   shape(1,0) =  2 * x * (x - y);
   shape(1,1) =  2 * (x - y) * (-1 + y);
   shape(2,0) =  2 * x * (-1 + 2 * x + y);
   shape(2,1) =  2 * y * (-1 + 2 * x + y);
   shape(3,0) =  2 * x * (-1 + x + 2 * y);
   shape(3,1) =  2 * y * (-1 + x + 2 * y);
   shape(4,0) = -2 * (-1 + x) * (x - y);
   shape(4,1) =  2 * y * (-x + y);
   shape(5,0) = -2 * (-1 + x) * (-1 + 2 * x + y);
   shape(5,1) = -2 * y * (-1 + 2 * x + y);
   shape(6,0) = -3 * x * (-2 + 2 * x + y);
   shape(6,1) = -3 * y * (-1 + 2 * x + y);
   shape(7,0) = -3 * x * (-1 + x + 2 * y);
   shape(7,1) = -3 * y * (-2 + x + 2 * y);
}

void RT1TriangleFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                            Vector &divshape) const
{
   double x = ip.x, y = ip.y;

   divshape(0) = -2 * (-4 + 3 * x + 6 * y);
   divshape(1) =  2 + 6 * x - 6 * y;
   divshape(2) = -4 + 12 * x + 6 * y;
   divshape(3) = -4 + 6 * x + 12 * y;
   divshape(4) =  2 - 6 * x + 6 * y;
   divshape(5) = -2 * (-4 + 6 * x + 3 * y);
   divshape(6) = -9 * (-1 + 2 * x + y);
   divshape(7) = -9 * (-1 + x + 2 * y);
}

const double RT1TriangleFiniteElement::nk[8][2] =
{
   { 0,-1}, { 0,-1},
   { 1, 1}, { 1, 1},
   {-1, 0}, {-1, 0},
   { 1, 0}, { 0, 1}
};

void RT1TriangleFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 8; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 8; j++)
      {
         double d = vshape(j,0)*nk[k][0]+vshape(j,1)*nk[k][1];
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT1QuadFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear (more to have embedding?)
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[2];
   Vector xk (vk, 2);

   for (k = 0; k < 8; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1];
      for (j = 0; j < 8; j++)
         if (fabs (I(k,j) = vshape(j,0)*vk[0]+vshape(j,1)*vk[1]) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT1TriangleFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   double vk[2];
   Vector xk (vk, 2);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 8; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) = (vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1] ) +
                 vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1] ));
      dofs(k) *= 0.5;
   }
}

RT1QuadFiniteElement::RT1QuadFiniteElement()
   : VectorFiniteElement(2, Geometry::SQUARE, 12, 2, H_DIV, FunctionSpace::Qk)
{
   // y = 0
   Nodes.IntPoint(0).x  = 1./3.;
   Nodes.IntPoint(0).y  = 0.0;
   Nodes.IntPoint(1).x  = 2./3.;
   Nodes.IntPoint(1).y  = 0.0;
   // x = 1
   Nodes.IntPoint(2).x  = 1.0;
   Nodes.IntPoint(2).y  = 1./3.;
   Nodes.IntPoint(3).x  = 1.0;
   Nodes.IntPoint(3).y  = 2./3.;
   // y = 1
   Nodes.IntPoint(4).x  = 2./3.;
   Nodes.IntPoint(4).y  = 1.0;
   Nodes.IntPoint(5).x  = 1./3.;
   Nodes.IntPoint(5).y  = 1.0;
   // x = 0
   Nodes.IntPoint(6).x  = 0.0;
   Nodes.IntPoint(6).y  = 2./3.;
   Nodes.IntPoint(7).x  = 0.0;
   Nodes.IntPoint(7).y  = 1./3.;
   // x = 0.5 (interior)
   Nodes.IntPoint(8).x  = 0.5;
   Nodes.IntPoint(8).y  = 1./3.;
   Nodes.IntPoint(9).x  = 0.5;
   Nodes.IntPoint(9).y  = 2./3.;
   // y = 0.5 (interior)
   Nodes.IntPoint(10).x = 1./3.;
   Nodes.IntPoint(10).y = 0.5;
   Nodes.IntPoint(11).x = 2./3.;
   Nodes.IntPoint(11).y = 0.5;
}

void RT1QuadFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   // y = 0
   shape(0,0)  = 0;
   shape(0,1)  = -( 1. - 3.*y + 2.*y*y)*( 2. - 3.*x);
   shape(1,0)  = 0;
   shape(1,1)  = -( 1. - 3.*y + 2.*y*y)*(-1. + 3.*x);
   // x = 1
   shape(2,0)  = (-x + 2.*x*x)*( 2. - 3.*y);
   shape(2,1)  = 0;
   shape(3,0)  = (-x + 2.*x*x)*(-1. + 3.*y);
   shape(3,1)  = 0;
   // y = 1
   shape(4,0)  = 0;
   shape(4,1)  = (-y + 2.*y*y)*(-1. + 3.*x);
   shape(5,0)  = 0;
   shape(5,1)  = (-y + 2.*y*y)*( 2. - 3.*x);
   // x = 0
   shape(6,0)  = -(1. - 3.*x + 2.*x*x)*(-1. + 3.*y);
   shape(6,1)  = 0;
   shape(7,0)  = -(1. - 3.*x + 2.*x*x)*( 2. - 3.*y);
   shape(7,1)  = 0;
   // x = 0.5 (interior)
   shape(8,0)  = (4.*x - 4.*x*x)*( 2. - 3.*y);
   shape(8,1)  = 0;
   shape(9,0)  = (4.*x - 4.*x*x)*(-1. + 3.*y);
   shape(9,1)  = 0;
   // y = 0.5 (interior)
   shape(10,0) = 0;
   shape(10,1) = (4.*y - 4.*y*y)*( 2. - 3.*x);
   shape(11,0) = 0;
   shape(11,1) = (4.*y - 4.*y*y)*(-1. + 3.*x);
}

void RT1QuadFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   double x = ip.x, y = ip.y;

   divshape(0)  = -(-3. + 4.*y)*( 2. - 3.*x);
   divshape(1)  = -(-3. + 4.*y)*(-1. + 3.*x);
   divshape(2)  = (-1. + 4.*x)*( 2. - 3.*y);
   divshape(3)  = (-1. + 4.*x)*(-1. + 3.*y);
   divshape(4)  = (-1. + 4.*y)*(-1. + 3.*x);
   divshape(5)  = (-1. + 4.*y)*( 2. - 3.*x);
   divshape(6)  = -(-3. + 4.*x)*(-1. + 3.*y);
   divshape(7)  = -(-3. + 4.*x)*( 2. - 3.*y);
   divshape(8)  = ( 4. - 8.*x)*( 2. - 3.*y);
   divshape(9)  = ( 4. - 8.*x)*(-1. + 3.*y);
   divshape(10) = ( 4. - 8.*y)*( 2. - 3.*x);
   divshape(11) = ( 4. - 8.*y)*(-1. + 3.*x);
}

const double RT1QuadFiniteElement::nk[12][2] =
{
   // y = 0
   {0,-1}, {0,-1},
   // X = 1
   {1, 0}, {1, 0},
   // y = 1
   {0, 1}, {0, 1},
   // x = 0
   {-1,0}, {-1,0},
   // x = 0.5 (interior)
   {1, 0}, {1, 0},
   // y = 0.5 (interior)
   {0, 1}, {0, 1}
};

void RT1QuadFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 12; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 12; j++)
      {
         double d = vshape(j,0)*nk[k][0]+vshape(j,1)*nk[k][1];
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT1QuadFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear (more to have embedding?)
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[2];
   Vector xk (vk, 2);

   for (k = 0; k < 12; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1];
      for (j = 0; j < 12; j++)
         if (fabs (I(k,j) = vshape(j,0)*vk[0]+vshape(j,1)*vk[1]) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT1QuadFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   double vk[2];
   Vector xk (vk, 2);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 12; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) = (vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1] ) +
                 vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1] ));
   }
}

const double RT2TriangleFiniteElement::M[15][15] =
{
   {
      0, -5.3237900077244501311, 5.3237900077244501311, 16.647580015448900262,
      0, 24.442740046346700787, -16.647580015448900262, -12.,
      -19.118950038622250656, -47.237900077244501311, 0, -34.414110069520051180,
      12., 30.590320061795601049, 15.295160030897800524
   },
   {
      0, 1.5, -1.5, -15., 0, 2.625, 15., 15., -4.125, 30., 0, -14.625, -15.,
      -15., 10.5
   },
   {
      0, -0.67620999227554986889, 0.67620999227554986889, 7.3524199845510997378,
      0, -3.4427400463467007866, -7.3524199845510997378, -12.,
      4.1189500386222506555, -0.76209992275549868892, 0, 7.4141100695200511800,
      12., -6.5903200617956010489, -3.2951600308978005244
   },
   {
      0, 0, 1.5, 0, 0, 1.5, -11.471370023173350393, 0, 2.4713700231733503933,
      -11.471370023173350393, 0, 2.4713700231733503933, 15.295160030897800524,
      0, -3.2951600308978005244
   },
   {
      0, 0, 4.875, 0, 0, 4.875, -16.875, 0, -16.875, -16.875, 0, -16.875, 10.5,
      36., 10.5
   },
   {
      0, 0, 1.5, 0, 0, 1.5, 2.4713700231733503933, 0, -11.471370023173350393,
      2.4713700231733503933, 0, -11.471370023173350393, -3.2951600308978005244,
      0, 15.295160030897800524
   },
   {
      -0.67620999227554986889, 0, -3.4427400463467007866, 0,
         7.3524199845510997378, 0.67620999227554986889, 7.4141100695200511800, 0,
         -0.76209992275549868892, 4.1189500386222506555, -12.,
         -7.3524199845510997378, -3.2951600308978005244, -6.5903200617956010489,
         12.
      },
   {
      1.5, 0, 2.625, 0, -15., -1.5, -14.625, 0, 30., -4.125, 15., 15., 10.5,
      -15., -15.
   },
   {
      -5.3237900077244501311, 0, 24.442740046346700787, 0, 16.647580015448900262,
         5.3237900077244501311, -34.414110069520051180, 0, -47.237900077244501311,
         -19.118950038622250656, -12., -16.647580015448900262, 15.295160030897800524,
         30.590320061795601049, 12.
      },
   { 0, 0, 18., 0, 0, 6., -42., 0, -30., -26., 0, -14., 24., 32., 8.},
   { 0, 0, 6., 0, 0, 18., -14., 0, -26., -30., 0, -42., 8., 32., 24.},
   { 0, 0, -6., 0, 0, -4., 30., 0, 4., 22., 0, 4., -24., -16., 0},
   { 0, 0, -4., 0, 0, -8., 20., 0, 8., 36., 0, 8., -16., -32., 0},
   { 0, 0, -8., 0, 0, -4., 8., 0, 36., 8., 0, 20., 0, -32., -16.},
   { 0, 0, -4., 0, 0, -6., 4., 0, 22., 4., 0, 30., 0, -16., -24.}
};

RT2TriangleFiniteElement::RT2TriangleFiniteElement()
   : VectorFiniteElement(2, Geometry::TRIANGLE, 15, 3, H_DIV)
{
   const double p = 0.11270166537925831148;

   Nodes.IntPoint(0).x = p;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 1.-p;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(3).x = 1.-p;
   Nodes.IntPoint(3).y = p;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(5).x = p;
   Nodes.IntPoint(5).y = 1.-p;
   Nodes.IntPoint(6).x = 0.0;
   Nodes.IntPoint(6).y = 1.-p;
   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = p;
   Nodes.IntPoint(9).x  = 0.25;
   Nodes.IntPoint(9).y  = 0.25;
   Nodes.IntPoint(10).x = 0.25;
   Nodes.IntPoint(10).y = 0.25;
   Nodes.IntPoint(11).x = 0.5;
   Nodes.IntPoint(11).y = 0.25;
   Nodes.IntPoint(12).x = 0.5;
   Nodes.IntPoint(12).y = 0.25;
   Nodes.IntPoint(13).x = 0.25;
   Nodes.IntPoint(13).y = 0.5;
   Nodes.IntPoint(14).x = 0.25;
   Nodes.IntPoint(14).y = 0.5;
}

void RT2TriangleFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                          DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   double Bx[15] = {1., 0., x, 0., y, 0., x*x, 0., x*y, 0., y*y, 0., x*x*x,
                    x*x*y, x*y*y
                   };
   double By[15] = {0., 1., 0., x, 0., y, 0., x*x, 0., x*y, 0., y*y,
                    x*x*y, x*y*y, y*y*y
                   };

   for (int i = 0; i < 15; i++)
   {
      double cx = 0.0, cy = 0.0;
      for (int j = 0; j < 15; j++)
      {
         cx += M[i][j] * Bx[j];
         cy += M[i][j] * By[j];
      }
      shape(i,0) = cx;
      shape(i,1) = cy;
   }
}

void RT2TriangleFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                            Vector &divshape) const
{
   double x = ip.x, y = ip.y;

   double DivB[15] = {0., 0., 1., 0., 0., 1., 2.*x, 0., y, x, 0., 2.*y,
                      4.*x*x, 4.*x*y, 4.*y*y
                     };

   for (int i = 0; i < 15; i++)
   {
      double div = 0.0;
      for (int j = 0; j < 15; j++)
      {
         div += M[i][j] * DivB[j];
      }
      divshape(i) = div;
   }
}

const double RT2QuadFiniteElement::pt[4] = {0.,1./3.,2./3.,1.};

const double RT2QuadFiniteElement::dpt[3] = {0.25,0.5,0.75};

RT2QuadFiniteElement::RT2QuadFiniteElement()
   : VectorFiniteElement(2, Geometry::SQUARE, 24, 3, H_DIV, FunctionSpace::Qk)
{
   // y = 0 (pt[0])
   Nodes.IntPoint(0).x  = dpt[0];  Nodes.IntPoint(0).y  =  pt[0];
   Nodes.IntPoint(1).x  = dpt[1];  Nodes.IntPoint(1).y  =  pt[0];
   Nodes.IntPoint(2).x  = dpt[2];  Nodes.IntPoint(2).y  =  pt[0];
   // x = 1 (pt[3])
   Nodes.IntPoint(3).x  =  pt[3];  Nodes.IntPoint(3).y  = dpt[0];
   Nodes.IntPoint(4).x  =  pt[3];  Nodes.IntPoint(4).y  = dpt[1];
   Nodes.IntPoint(5).x  =  pt[3];  Nodes.IntPoint(5).y  = dpt[2];
   // y = 1 (pt[3])
   Nodes.IntPoint(6).x  = dpt[2];  Nodes.IntPoint(6).y  =  pt[3];
   Nodes.IntPoint(7).x  = dpt[1];  Nodes.IntPoint(7).y  =  pt[3];
   Nodes.IntPoint(8).x  = dpt[0];  Nodes.IntPoint(8).y  =  pt[3];
   // x = 0 (pt[0])
   Nodes.IntPoint(9).x  =  pt[0];  Nodes.IntPoint(9).y  = dpt[2];
   Nodes.IntPoint(10).x =  pt[0];  Nodes.IntPoint(10).y = dpt[1];
   Nodes.IntPoint(11).x =  pt[0];  Nodes.IntPoint(11).y = dpt[0];
   // x = pt[1] (interior)
   Nodes.IntPoint(12).x =  pt[1];  Nodes.IntPoint(12).y = dpt[0];
   Nodes.IntPoint(13).x =  pt[1];  Nodes.IntPoint(13).y = dpt[1];
   Nodes.IntPoint(14).x =  pt[1];  Nodes.IntPoint(14).y = dpt[2];
   // x = pt[2] (interior)
   Nodes.IntPoint(15).x =  pt[2];  Nodes.IntPoint(15).y = dpt[0];
   Nodes.IntPoint(16).x =  pt[2];  Nodes.IntPoint(16).y = dpt[1];
   Nodes.IntPoint(17).x =  pt[2];  Nodes.IntPoint(17).y = dpt[2];
   // y = pt[1] (interior)
   Nodes.IntPoint(18).x = dpt[0];  Nodes.IntPoint(18).y =  pt[1];
   Nodes.IntPoint(19).x = dpt[1];  Nodes.IntPoint(19).y =  pt[1];
   Nodes.IntPoint(20).x = dpt[2];  Nodes.IntPoint(20).y =  pt[1];
   // y = pt[2] (interior)
   Nodes.IntPoint(21).x = dpt[0];  Nodes.IntPoint(21).y =  pt[2];
   Nodes.IntPoint(22).x = dpt[1];  Nodes.IntPoint(22).y =  pt[2];
   Nodes.IntPoint(23).x = dpt[2];  Nodes.IntPoint(23).y =  pt[2];
}

void RT2QuadFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y;

   double ax0 =  pt[0] - x;
   double ax1 =  pt[1] - x;
   double ax2 =  pt[2] - x;
   double ax3 =  pt[3] - x;

   double by0 = dpt[0] - y;
   double by1 = dpt[1] - y;
   double by2 = dpt[2] - y;

   double ay0 =  pt[0] - y;
   double ay1 =  pt[1] - y;
   double ay2 =  pt[2] - y;
   double ay3 =  pt[3] - y;

   double bx0 = dpt[0] - x;
   double bx1 = dpt[1] - x;
   double bx2 = dpt[2] - x;

   double A01 =  pt[0] -  pt[1];
   double A02 =  pt[0] -  pt[2];
   double A12 =  pt[1] -  pt[2];
   double A03 =  pt[0] -  pt[3];
   double A13 =  pt[1] -  pt[3];
   double A23 =  pt[2] -  pt[3];

   double B01 = dpt[0] - dpt[1];
   double B02 = dpt[0] - dpt[2];
   double B12 = dpt[1] - dpt[2];

   double tx0 =  (bx1*bx2)/(B01*B02);
   double tx1 = -(bx0*bx2)/(B01*B12);
   double tx2 =  (bx0*bx1)/(B02*B12);

   double ty0 =  (by1*by2)/(B01*B02);
   double ty1 = -(by0*by2)/(B01*B12);
   double ty2 =  (by0*by1)/(B02*B12);

   // y = 0 (p[0])
   shape(0,  0) =  0;
   shape(0,  1) =  (ay1*ay2*ay3)/(A01*A02*A03)*tx0;
   shape(1,  0) =  0;
   shape(1,  1) =  (ay1*ay2*ay3)/(A01*A02*A03)*tx1;
   shape(2,  0) =  0;
   shape(2,  1) =  (ay1*ay2*ay3)/(A01*A02*A03)*tx2;
   // x = 1 (p[3])
   shape(3,  0) =  (ax0*ax1*ax2)/(A03*A13*A23)*ty0;
   shape(3,  1) =  0;
   shape(4,  0) =  (ax0*ax1*ax2)/(A03*A13*A23)*ty1;
   shape(4,  1) =  0;
   shape(5,  0) =  (ax0*ax1*ax2)/(A03*A13*A23)*ty2;
   shape(5,  1) =  0;
   // y = 1 (p[3])
   shape(6,  0) =  0;
   shape(6,  1) =  (ay0*ay1*ay2)/(A03*A13*A23)*tx2;
   shape(7,  0) =  0;
   shape(7,  1) =  (ay0*ay1*ay2)/(A03*A13*A23)*tx1;
   shape(8,  0) =  0;
   shape(8,  1) =  (ay0*ay1*ay2)/(A03*A13*A23)*tx0;
   // x = 0 (p[0])
   shape(9,  0) =  (ax1*ax2*ax3)/(A01*A02*A03)*ty2;
   shape(9,  1) =  0;
   shape(10, 0) =  (ax1*ax2*ax3)/(A01*A02*A03)*ty1;
   shape(10, 1) =  0;
   shape(11, 0) =  (ax1*ax2*ax3)/(A01*A02*A03)*ty0;
   shape(11, 1) =  0;
   // x = p[1] (interior)
   shape(12, 0) =  (ax0*ax2*ax3)/(A01*A12*A13)*ty0;
   shape(12, 1) =  0;
   shape(13, 0) =  (ax0*ax2*ax3)/(A01*A12*A13)*ty1;
   shape(13, 1) =  0;
   shape(14, 0) =  (ax0*ax2*ax3)/(A01*A12*A13)*ty2;
   shape(14, 1) =  0;
   // x = p[2] (interior)
   shape(15, 0) = -(ax0*ax1*ax3)/(A02*A12*A23)*ty0;
   shape(15, 1) =  0;
   shape(16, 0) = -(ax0*ax1*ax3)/(A02*A12*A23)*ty1;
   shape(16, 1) =  0;
   shape(17, 0) = -(ax0*ax1*ax3)/(A02*A12*A23)*ty2;
   shape(17, 1) =  0;
   // y = p[1] (interior)
   shape(18, 0) =  0;
   shape(18, 1) =  (ay0*ay2*ay3)/(A01*A12*A13)*tx0;
   shape(19, 0) =  0;
   shape(19, 1) =  (ay0*ay2*ay3)/(A01*A12*A13)*tx1;
   shape(20, 0) =  0;
   shape(20, 1) =  (ay0*ay2*ay3)/(A01*A12*A13)*tx2;
   // y = p[2] (interior)
   shape(21, 0) =  0;
   shape(21, 1) = -(ay0*ay1*ay3)/(A02*A12*A23)*tx0;
   shape(22, 0) =  0;
   shape(22, 1) = -(ay0*ay1*ay3)/(A02*A12*A23)*tx1;
   shape(23, 0) =  0;
   shape(23, 1) = -(ay0*ay1*ay3)/(A02*A12*A23)*tx2;
}

void RT2QuadFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   double x = ip.x, y = ip.y;

   double a01 =  pt[0]*pt[1];
   double a02 =  pt[0]*pt[2];
   double a12 =  pt[1]*pt[2];
   double a03 =  pt[0]*pt[3];
   double a13 =  pt[1]*pt[3];
   double a23 =  pt[2]*pt[3];

   double bx0 = dpt[0] - x;
   double bx1 = dpt[1] - x;
   double bx2 = dpt[2] - x;

   double by0 = dpt[0] - y;
   double by1 = dpt[1] - y;
   double by2 = dpt[2] - y;

   double A01 =  pt[0] -  pt[1];
   double A02 =  pt[0] -  pt[2];
   double A12 =  pt[1] -  pt[2];
   double A03 =  pt[0] -  pt[3];
   double A13 =  pt[1] -  pt[3];
   double A23 =  pt[2] -  pt[3];

   double A012 = pt[0] + pt[1] + pt[2];
   double A013 = pt[0] + pt[1] + pt[3];
   double A023 = pt[0] + pt[2] + pt[3];
   double A123 = pt[1] + pt[2] + pt[3];

   double B01 = dpt[0] - dpt[1];
   double B02 = dpt[0] - dpt[2];
   double B12 = dpt[1] - dpt[2];

   double tx0 =  (bx1*bx2)/(B01*B02);
   double tx1 = -(bx0*bx2)/(B01*B12);
   double tx2 =  (bx0*bx1)/(B02*B12);

   double ty0 =  (by1*by2)/(B01*B02);
   double ty1 = -(by0*by2)/(B01*B12);
   double ty2 =  (by0*by1)/(B02*B12);

   // y = 0 (p[0])
   divshape(0)  = -(a12 + a13 + a23 - 2.*A123*y + 3.*y*y)/(A01*A02*A03)*tx0;
   divshape(1)  = -(a12 + a13 + a23 - 2.*A123*y + 3.*y*y)/(A01*A02*A03)*tx1;
   divshape(2)  = -(a12 + a13 + a23 - 2.*A123*y + 3.*y*y)/(A01*A02*A03)*tx2;
   // x = 1 (p[3])
   divshape(3)  = -(a01 + a02 + a12 - 2.*A012*x + 3.*x*x)/(A03*A13*A23)*ty0;
   divshape(4)  = -(a01 + a02 + a12 - 2.*A012*x + 3.*x*x)/(A03*A13*A23)*ty1;
   divshape(5)  = -(a01 + a02 + a12 - 2.*A012*x + 3.*x*x)/(A03*A13*A23)*ty2;
   // y = 1 (p[3])
   divshape(6)  = -(a01 + a02 + a12 - 2.*A012*y + 3.*y*y)/(A03*A13*A23)*tx2;
   divshape(7)  = -(a01 + a02 + a12 - 2.*A012*y + 3.*y*y)/(A03*A13*A23)*tx1;
   divshape(8)  = -(a01 + a02 + a12 - 2.*A012*y + 3.*y*y)/(A03*A13*A23)*tx0;
   // x = 0 (p[0])
   divshape(9)  = -(a12 + a13 + a23 - 2.*A123*x + 3.*x*x)/(A01*A02*A03)*ty2;
   divshape(10) = -(a12 + a13 + a23 - 2.*A123*x + 3.*x*x)/(A01*A02*A03)*ty1;
   divshape(11) = -(a12 + a13 + a23 - 2.*A123*x + 3.*x*x)/(A01*A02*A03)*ty0;
   // x = p[1] (interior)
   divshape(12) = -(a02 + a03 + a23 - 2.*A023*x + 3.*x*x)/(A01*A12*A13)*ty0;
   divshape(13) = -(a02 + a03 + a23 - 2.*A023*x + 3.*x*x)/(A01*A12*A13)*ty1;
   divshape(14) = -(a02 + a03 + a23 - 2.*A023*x + 3.*x*x)/(A01*A12*A13)*ty2;
   // x = p[2] (interior)
   divshape(15) =  (a01 + a03 + a13 - 2.*A013*x + 3.*x*x)/(A02*A12*A23)*ty0;
   divshape(16) =  (a01 + a03 + a13 - 2.*A013*x + 3.*x*x)/(A02*A12*A23)*ty1;
   divshape(17) =  (a01 + a03 + a13 - 2.*A013*x + 3.*x*x)/(A02*A12*A23)*ty2;
   // y = p[1] (interior)
   divshape(18) = -(a02 + a03 + a23 - 2.*A023*y + 3.*y*y)/(A01*A12*A13)*tx0;
   divshape(19) = -(a02 + a03 + a23 - 2.*A023*y + 3.*y*y)/(A01*A12*A13)*tx1;
   divshape(20) = -(a02 + a03 + a23 - 2.*A023*y + 3.*y*y)/(A01*A12*A13)*tx2;
   // y = p[2] (interior)
   divshape(21) =  (a01 + a03 + a13 - 2.*A013*y + 3.*y*y)/(A02*A12*A23)*tx0;
   divshape(22) =  (a01 + a03 + a13 - 2.*A013*y + 3.*y*y)/(A02*A12*A23)*tx1;
   divshape(23) =  (a01 + a03 + a13 - 2.*A013*y + 3.*y*y)/(A02*A12*A23)*tx2;
}

const double RT2QuadFiniteElement::nk[24][2] =
{
   // y = 0
   {0,-1}, {0,-1}, {0,-1},
   // x = 1
   {1, 0}, {1, 0}, {1, 0},
   // y = 1
   {0, 1}, {0, 1}, {0, 1},
   // x = 0
   {-1,0}, {-1,0}, {-1,0},
   // x = p[1] (interior)
   {1, 0}, {1, 0}, {1, 0},
   // x = p[2] (interior)
   {1, 0}, {1, 0}, {1, 0},
   // y = p[1] (interior)
   {0, 1}, {0, 1}, {0, 1},
   // y = p[1] (interior)
   {0, 1}, {0, 1}, {0, 1}
};

void RT2QuadFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 24; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 24; j++)
      {
         double d = vshape(j,0)*nk[k][0]+vshape(j,1)*nk[k][1];
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT2QuadFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear (more to have embedding?)
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[2];
   Vector xk (vk, 2);

   for (k = 0; k < 24; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1];
      for (j = 0; j < 24; j++)
         if (fabs (I(k,j) = vshape(j,0)*vk[0]+vshape(j,1)*vk[1]) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT2QuadFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans, Vector &dofs) const
{
   double vk[2];
   Vector xk (vk, 2);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 24; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) = (vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1] ) +
                 vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1] ));
   }
}

P1SegmentFiniteElement::P1SegmentFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 2, 1)
{
   Nodes.IntPoint(0).x = 0.33333333333333333333;
   Nodes.IntPoint(1).x = 0.66666666666666666667;
}

void P1SegmentFiniteElement::CalcShape(const IntegrationPoint &ip,
                                       Vector &shape) const
{
   double x = ip.x;

   shape(0) = 2. - 3. * x;
   shape(1) = 3. * x - 1.;
}

void P1SegmentFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                        DenseMatrix &dshape) const
{
   dshape(0,0) = -3.;
   dshape(1,0) =  3.;
}


P2SegmentFiniteElement::P2SegmentFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 3, 2)
{
   const double p = 0.11270166537925831148;

   Nodes.IntPoint(0).x = p;
   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(2).x = 1.-p;
}

void P2SegmentFiniteElement::CalcShape(const IntegrationPoint &ip,
                                       Vector &shape) const
{
   const double p = 0.11270166537925831148;
   const double w = 1./((1-2*p)*(1-2*p));
   double x = ip.x;

   shape(0) = (2*x-1)*(x-1+p)*w;
   shape(1) = 4*(x-1+p)*(p-x)*w;
   shape(2) = (2*x-1)*(x-p)*w;
}

void P2SegmentFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                        DenseMatrix &dshape) const
{
   const double p = 0.11270166537925831148;
   const double w = 1./((1-2*p)*(1-2*p));
   double x = ip.x;

   dshape(0,0) = (-3+4*x+2*p)*w;
   dshape(1,0) = (4-8*x)*w;
   dshape(2,0) = (-1+4*x-2*p)*w;
}


Lagrange1DFiniteElement::Lagrange1DFiniteElement(int degree)
   : NodalFiniteElement(1, Geometry::SEGMENT, degree+1, degree)
{
   int i, m = degree;

   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   for (i = 1; i < m; i++)
   {
      Nodes.IntPoint(i+1).x = double(i) / m;
   }

   rwk.SetSize(degree+1);
#ifndef MFEM_THREAD_SAFE
   rxxk.SetSize(degree+1);
#endif

   rwk(0) = 1.0;
   for (i = 1; i <= m; i++)
   {
      rwk(i) = rwk(i-1) * ( (double)(m) / (double)(i) );
   }
   for (i = 0; i < m/2+1; i++)
   {
      rwk(m-i) = ( rwk(i) *= rwk(m-i) );
   }
   for (i = m-1; i >= 0; i -= 2)
   {
      rwk(i) = -rwk(i);
   }
}

void Lagrange1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   double w, wk, x = ip.x;
   int i, k, m = GetOrder();

#ifdef MFEM_THREAD_SAFE
   Vector rxxk(m+1);
#endif

   k = (int) floor ( m * x + 0.5 );
   k = k > m ? m : k < 0 ? 0 : k; // clamp k to [0,m]

   wk = 1.0;
   for (i = 0; i <= m; i++)
      if (i != k)
      {
         wk *= ( rxxk(i) = x - (double)(i) / m );
      }
   w = wk * ( rxxk(k) = x - (double)(k) / m );

   if (k != 0)
   {
      shape(0) = w * rwk(0) / rxxk(0);
   }
   else
   {
      shape(0) = wk * rwk(0);
   }
   if (k != m)
   {
      shape(1) = w * rwk(m) / rxxk(m);
   }
   else
   {
      shape(1) = wk * rwk(k);
   }
   for (i = 1; i < m; i++)
      if (i != k)
      {
         shape(i+1) = w * rwk(i) / rxxk(i);
      }
      else
      {
         shape(k+1) = wk * rwk(k);
      }
}

void Lagrange1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   double s, srx, w, wk, x = ip.x;
   int i, k, m = GetOrder();

#ifdef MFEM_THREAD_SAFE
   Vector rxxk(m+1);
#endif

   k = (int) floor ( m * x + 0.5 );
   k = k > m ? m : k < 0 ? 0 : k; // clamp k to [0,m]

   wk = 1.0;
   for (i = 0; i <= m; i++)
      if (i != k)
      {
         wk *= ( rxxk(i) = x - (double)(i) / m );
      }
   w = wk * ( rxxk(k) = x - (double)(k) / m );

   for (i = 0; i <= m; i++)
   {
      rxxk(i) = 1.0 / rxxk(i);
   }
   srx = 0.0;
   for (i = 0; i <= m; i++)
      if (i != k)
      {
         srx += rxxk(i);
      }
   s = w * srx + wk;

   if (k != 0)
   {
      dshape(0,0) = (s - w * rxxk(0)) * rwk(0) * rxxk(0);
   }
   else
   {
      dshape(0,0) = wk * srx * rwk(0);
   }
   if (k != m)
   {
      dshape(1,0) = (s - w * rxxk(m)) * rwk(m) * rxxk(m);
   }
   else
   {
      dshape(1,0) = wk * srx * rwk(k);
   }
   for (i = 1; i < m; i++)
      if (i != k)
      {
         dshape(i+1,0) = (s - w * rxxk(i)) * rwk(i) * rxxk(i);
      }
      else
      {
         dshape(k+1,0) = wk * srx * rwk(k);
      }
}


P1TetNonConfFiniteElement::P1TetNonConfFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 4, 1)
{
   Nodes.IntPoint(0).x = 0.33333333333333333333;
   Nodes.IntPoint(0).y = 0.33333333333333333333;
   Nodes.IntPoint(0).z = 0.33333333333333333333;

   Nodes.IntPoint(1).x = 0.0;
   Nodes.IntPoint(1).y = 0.33333333333333333333;
   Nodes.IntPoint(1).z = 0.33333333333333333333;

   Nodes.IntPoint(2).x = 0.33333333333333333333;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.33333333333333333333;

   Nodes.IntPoint(3).x = 0.33333333333333333333;
   Nodes.IntPoint(3).y = 0.33333333333333333333;
   Nodes.IntPoint(3).z = 0.0;

}

void P1TetNonConfFiniteElement::CalcShape(const IntegrationPoint &ip,
                                          Vector &shape) const
{
   double L0, L1, L2, L3;

   L1 = ip.x;  L2 = ip.y;  L3 = ip.z;  L0 = 1.0 - L1 - L2 - L3;
   shape(0) = 1.0 - 3.0 * L0;
   shape(1) = 1.0 - 3.0 * L1;
   shape(2) = 1.0 - 3.0 * L2;
   shape(3) = 1.0 - 3.0 * L3;
}

void P1TetNonConfFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                           DenseMatrix &dshape) const
{
   dshape(0,0) =  3.0; dshape(0,1) =  3.0; dshape(0,2) =  3.0;
   dshape(1,0) = -3.0; dshape(1,1) =  0.0; dshape(1,2) =  0.0;
   dshape(2,0) =  0.0; dshape(2,1) = -3.0; dshape(2,2) =  0.0;
   dshape(3,0) =  0.0; dshape(3,1) =  0.0; dshape(3,2) = -3.0;
}


P0TetFiniteElement::P0TetFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 1, 0)
{
   Nodes.IntPoint(0).x = 0.25;
   Nodes.IntPoint(0).y = 0.25;
   Nodes.IntPoint(0).z = 0.25;
}

void P0TetFiniteElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
   shape(0) = 1.0;
}

void P0TetFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
   dshape(0,0) =  0.0; dshape(0,1) =  0.0; dshape(0,2) = 0.0;
}


P0HexFiniteElement::P0HexFiniteElement()
   : NodalFiniteElement(3, Geometry::CUBE, 1, 0, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.5;
   Nodes.IntPoint(0).z = 0.5;
}

void P0HexFiniteElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
   shape(0) = 1.0;
}

void P0HexFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
   dshape(0,0) =  0.0; dshape(0,1) =  0.0; dshape(0,2) = 0.0;
}


LagrangeHexFiniteElement::LagrangeHexFiniteElement (int degree)
   : NodalFiniteElement(3, Geometry::CUBE, (degree+1)*(degree+1)*(degree+1),
                        degree, FunctionSpace::Qk)
{
   if (degree == 2)
   {
      I = new int[dof];
      J = new int[dof];
      K = new int[dof];
      // nodes
      I[ 0] = 0; J[ 0] = 0; K[ 0] = 0;
      I[ 1] = 1; J[ 1] = 0; K[ 1] = 0;
      I[ 2] = 1; J[ 2] = 1; K[ 2] = 0;
      I[ 3] = 0; J[ 3] = 1; K[ 3] = 0;
      I[ 4] = 0; J[ 4] = 0; K[ 4] = 1;
      I[ 5] = 1; J[ 5] = 0; K[ 5] = 1;
      I[ 6] = 1; J[ 6] = 1; K[ 6] = 1;
      I[ 7] = 0; J[ 7] = 1; K[ 7] = 1;
      // edges
      I[ 8] = 2; J[ 8] = 0; K[ 8] = 0;
      I[ 9] = 1; J[ 9] = 2; K[ 9] = 0;
      I[10] = 2; J[10] = 1; K[10] = 0;
      I[11] = 0; J[11] = 2; K[11] = 0;
      I[12] = 2; J[12] = 0; K[12] = 1;
      I[13] = 1; J[13] = 2; K[13] = 1;
      I[14] = 2; J[14] = 1; K[14] = 1;
      I[15] = 0; J[15] = 2; K[15] = 1;
      I[16] = 0; J[16] = 0; K[16] = 2;
      I[17] = 1; J[17] = 0; K[17] = 2;
      I[18] = 1; J[18] = 1; K[18] = 2;
      I[19] = 0; J[19] = 1; K[19] = 2;
      // faces
      I[20] = 2; J[20] = 2; K[20] = 0;
      I[21] = 2; J[21] = 0; K[21] = 2;
      I[22] = 1; J[22] = 2; K[22] = 2;
      I[23] = 2; J[23] = 1; K[23] = 2;
      I[24] = 0; J[24] = 2; K[24] = 2;
      I[25] = 2; J[25] = 2; K[25] = 1;
      // element
      I[26] = 2; J[26] = 2; K[26] = 2;
   }
   else if (degree == 3)
   {
      I = new int[dof];
      J = new int[dof];
      K = new int[dof];
      // nodes
      I[ 0] = 0; J[ 0] = 0; K[ 0] = 0;
      I[ 1] = 1; J[ 1] = 0; K[ 1] = 0;
      I[ 2] = 1; J[ 2] = 1; K[ 2] = 0;
      I[ 3] = 0; J[ 3] = 1; K[ 3] = 0;
      I[ 4] = 0; J[ 4] = 0; K[ 4] = 1;
      I[ 5] = 1; J[ 5] = 0; K[ 5] = 1;
      I[ 6] = 1; J[ 6] = 1; K[ 6] = 1;
      I[ 7] = 0; J[ 7] = 1; K[ 7] = 1;
      // edges
      I[ 8] = 2; J[ 8] = 0; K[ 8] = 0;
      I[ 9] = 3; J[ 9] = 0; K[ 9] = 0;
      I[10] = 1; J[10] = 2; K[10] = 0;
      I[11] = 1; J[11] = 3; K[11] = 0;
      I[12] = 2; J[12] = 1; K[12] = 0;
      I[13] = 3; J[13] = 1; K[13] = 0;
      I[14] = 0; J[14] = 2; K[14] = 0;
      I[15] = 0; J[15] = 3; K[15] = 0;
      I[16] = 2; J[16] = 0; K[16] = 1;
      I[17] = 3; J[17] = 0; K[17] = 1;
      I[18] = 1; J[18] = 2; K[18] = 1;
      I[19] = 1; J[19] = 3; K[19] = 1;
      I[20] = 2; J[20] = 1; K[20] = 1;
      I[21] = 3; J[21] = 1; K[21] = 1;
      I[22] = 0; J[22] = 2; K[22] = 1;
      I[23] = 0; J[23] = 3; K[23] = 1;
      I[24] = 0; J[24] = 0; K[24] = 2;
      I[25] = 0; J[25] = 0; K[25] = 3;
      I[26] = 1; J[26] = 0; K[26] = 2;
      I[27] = 1; J[27] = 0; K[27] = 3;
      I[28] = 1; J[28] = 1; K[28] = 2;
      I[29] = 1; J[29] = 1; K[29] = 3;
      I[30] = 0; J[30] = 1; K[30] = 2;
      I[31] = 0; J[31] = 1; K[31] = 3;
      // faces
      I[32] = 2; J[32] = 3; K[32] = 0;
      I[33] = 3; J[33] = 3; K[33] = 0;
      I[34] = 2; J[34] = 2; K[34] = 0;
      I[35] = 3; J[35] = 2; K[35] = 0;
      I[36] = 2; J[36] = 0; K[36] = 2;
      I[37] = 3; J[37] = 0; K[37] = 2;
      I[38] = 2; J[38] = 0; K[38] = 3;
      I[39] = 3; J[39] = 0; K[39] = 3;
      I[40] = 1; J[40] = 2; K[40] = 2;
      I[41] = 1; J[41] = 3; K[41] = 2;
      I[42] = 1; J[42] = 2; K[42] = 3;
      I[43] = 1; J[43] = 3; K[43] = 3;
      I[44] = 3; J[44] = 1; K[44] = 2;
      I[45] = 2; J[45] = 1; K[45] = 2;
      I[46] = 3; J[46] = 1; K[46] = 3;
      I[47] = 2; J[47] = 1; K[47] = 3;
      I[48] = 0; J[48] = 3; K[48] = 2;
      I[49] = 0; J[49] = 2; K[49] = 2;
      I[50] = 0; J[50] = 3; K[50] = 3;
      I[51] = 0; J[51] = 2; K[51] = 3;
      I[52] = 2; J[52] = 2; K[52] = 1;
      I[53] = 3; J[53] = 2; K[53] = 1;
      I[54] = 2; J[54] = 3; K[54] = 1;
      I[55] = 3; J[55] = 3; K[55] = 1;
      // element
      I[56] = 2; J[56] = 2; K[56] = 2;
      I[57] = 3; J[57] = 2; K[57] = 2;
      I[58] = 3; J[58] = 3; K[58] = 2;
      I[59] = 2; J[59] = 3; K[59] = 2;
      I[60] = 2; J[60] = 2; K[60] = 3;
      I[61] = 3; J[61] = 2; K[61] = 3;
      I[62] = 3; J[62] = 3; K[62] = 3;
      I[63] = 2; J[63] = 3; K[63] = 3;
   }
   else
   {
      mfem_error ("LagrangeHexFiniteElement::LagrangeHexFiniteElement");
   }

   fe1d = new Lagrange1DFiniteElement(degree);
   dof1d = fe1d -> GetDof();

#ifndef MFEM_THREAD_SAFE
   shape1dx.SetSize(dof1d);
   shape1dy.SetSize(dof1d);
   shape1dz.SetSize(dof1d);

   dshape1dx.SetSize(dof1d,1);
   dshape1dy.SetSize(dof1d,1);
   dshape1dz.SetSize(dof1d,1);
#endif

   for (int n = 0; n < dof; n++)
   {
      Nodes.IntPoint(n).x = fe1d -> GetNodes().IntPoint(I[n]).x;
      Nodes.IntPoint(n).y = fe1d -> GetNodes().IntPoint(J[n]).x;
      Nodes.IntPoint(n).z = fe1d -> GetNodes().IntPoint(K[n]).x;
   }
}

void LagrangeHexFiniteElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   IntegrationPoint ipy, ipz;
   ipy.x = ip.y;
   ipz.x = ip.z;

#ifdef MFEM_THREAD_SAFE
   Vector shape1dx(dof1d), shape1dy(dof1d), shape1dz(dof1d);
#endif

   fe1d -> CalcShape(ip,  shape1dx);
   fe1d -> CalcShape(ipy, shape1dy);
   fe1d -> CalcShape(ipz, shape1dz);

   for (int n = 0; n < dof; n++)
   {
      shape(n) = shape1dx(I[n]) *  shape1dy(J[n]) * shape1dz(K[n]);
   }
}

void LagrangeHexFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
   IntegrationPoint ipy, ipz;
   ipy.x = ip.y;
   ipz.x = ip.z;

#ifdef MFEM_THREAD_SAFE
   Vector shape1dx(dof1d), shape1dy(dof1d), shape1dz(dof1d);
   DenseMatrix dshape1dx(dof1d,1), dshape1dy(dof1d,1), dshape1dz(dof1d,1);
#endif

   fe1d -> CalcShape(ip,  shape1dx);
   fe1d -> CalcShape(ipy, shape1dy);
   fe1d -> CalcShape(ipz, shape1dz);

   fe1d -> CalcDShape(ip,  dshape1dx);
   fe1d -> CalcDShape(ipy, dshape1dy);
   fe1d -> CalcDShape(ipz, dshape1dz);

   for (int n = 0; n < dof; n++)
   {
      dshape(n,0) = dshape1dx(I[n],0) * shape1dy(J[n])    * shape1dz(K[n]);
      dshape(n,1) = shape1dx(I[n])    * dshape1dy(J[n],0) * shape1dz(K[n]);
      dshape(n,2) = shape1dx(I[n])    * shape1dy(J[n])    * dshape1dz(K[n],0);
   }
}

LagrangeHexFiniteElement::~LagrangeHexFiniteElement ()
{
   delete fe1d;

   delete [] I;
   delete [] J;
   delete [] K;
}


RefinedLinear1DFiniteElement::RefinedLinear1DFiniteElement()
   : NodalFiniteElement(1, Geometry::SEGMENT, 3, 4)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(2).x = 0.5;
}

void RefinedLinear1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   double x = ip.x;

   if (x <= 0.5)
   {
      shape(0) = 1.0 - 2.0 * x;
      shape(1) = 0.0;
      shape(2) = 2.0 * x;
   }
   else
   {
      shape(0) = 0.0;
      shape(1) = 2.0 * x - 1.0;
      shape(2) = 2.0 - 2.0 * x;
   }
}

void RefinedLinear1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   double x = ip.x;

   if (x <= 0.5)
   {
      dshape(0,0) = - 2.0;
      dshape(1,0) =   0.0;
      dshape(2,0) =   2.0;
   }
   else
   {
      dshape(0,0) =   0.0;
      dshape(1,0) =   2.0;
      dshape(2,0) = - 2.0;
   }
}

RefinedLinear2DFiniteElement::RefinedLinear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::TRIANGLE, 6, 5)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.5;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.5;
}

void RefinedLinear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   int i;

   double L0, L1, L2;
   L0 = 2.0 * ( 1. - ip.x - ip.y );
   L1 = 2.0 * ( ip.x );
   L2 = 2.0 * ( ip.y );

   // The reference triangle is split in 4 triangles as follows:
   //
   // T0 - 0,3,5
   // T1 - 1,3,4
   // T2 - 2,4,5
   // T3 - 3,4,5

   for (i = 0; i < 6; i++)
   {
      shape(i) = 0.0;
   }

   if (L0 >= 1.0)   // T0
   {
      shape(0) = L0 - 1.0;
      shape(3) =       L1;
      shape(5) =       L2;
   }
   else if (L1 >= 1.0)   // T1
   {
      shape(3) =       L0;
      shape(1) = L1 - 1.0;
      shape(4) =       L2;
   }
   else if (L2 >= 1.0)   // T2
   {
      shape(5) =       L0;
      shape(4) =       L1;
      shape(2) = L2 - 1.0;
   }
   else   // T3
   {
      shape(3) = 1.0 - L2;
      shape(4) = 1.0 - L0;
      shape(5) = 1.0 - L1;
   }
}

void RefinedLinear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   int i,j;

   double L0, L1, L2;
   L0 = 2.0 * ( 1. - ip.x - ip.y );
   L1 = 2.0 * ( ip.x );
   L2 = 2.0 * ( ip.y );

   double DL0[2], DL1[2], DL2[2];
   DL0[0] = -2.0; DL0[1] = -2.0;
   DL1[0] =  2.0; DL1[1] =  0.0;
   DL2[0] =  0.0; DL2[1] =  2.0;

   for (i = 0; i < 6; i++)
      for (j = 0; j < 2; j++)
      {
         dshape(i,j) = 0.0;
      }

   if (L0 >= 1.0)   // T0
   {
      for (j = 0; j < 2; j++)
      {
         dshape(0,j) = DL0[j];
         dshape(3,j) = DL1[j];
         dshape(5,j) = DL2[j];
      }
   }
   else if (L1 >= 1.0)   // T1
   {
      for (j = 0; j < 2; j++)
      {
         dshape(3,j) = DL0[j];
         dshape(1,j) = DL1[j];
         dshape(4,j) = DL2[j];
      }
   }
   else if (L2 >= 1.0)   // T2
   {
      for (j = 0; j < 2; j++)
      {
         dshape(5,j) = DL0[j];
         dshape(4,j) = DL1[j];
         dshape(2,j) = DL2[j];
      }
   }
   else   // T3
   {
      for (j = 0; j < 2; j++)
      {
         dshape(3,j) = - DL2[j];
         dshape(4,j) = - DL0[j];
         dshape(5,j) = - DL1[j];
      }
   }
}

RefinedLinear3DFiniteElement::RefinedLinear3DFiniteElement()
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, 10, 4)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 1.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 0.0;
   Nodes.IntPoint(6).x = 0.0;
   Nodes.IntPoint(6).y = 0.0;
   Nodes.IntPoint(6).z = 0.5;
   Nodes.IntPoint(7).x = 0.5;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(7).z = 0.0;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.5;
   Nodes.IntPoint(9).x = 0.0;
   Nodes.IntPoint(9).y = 0.5;
   Nodes.IntPoint(9).z = 0.5;
}

void RefinedLinear3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   int i;

   double L0, L1, L2, L3, L4, L5;
   L0 = 2.0 * ( 1. - ip.x - ip.y - ip.z );
   L1 = 2.0 * ( ip.x );
   L2 = 2.0 * ( ip.y );
   L3 = 2.0 * ( ip.z );
   L4 = 2.0 * ( ip.x + ip.y );
   L5 = 2.0 * ( ip.y + ip.z );

   // The reference tetrahedron is split in 8 tetrahedra as follows:
   //
   // T0 - 0,4,5,6
   // T1 - 1,4,7,8
   // T2 - 2,5,7,9
   // T3 - 3,6,8,9
   // T4 - 4,5,6,8
   // T5 - 4,5,7,8
   // T6 - 5,6,8,9
   // T7 - 5,7,8,9

   for (i = 0; i < 10; i++)
   {
      shape(i) = 0.0;
   }

   if (L0 >= 1.0)   // T0
   {
      shape(0) = L0 - 1.0;
      shape(4) =       L1;
      shape(5) =       L2;
      shape(6) =       L3;
   }
   else if (L1 >= 1.0)   // T1
   {
      shape(4) =       L0;
      shape(1) = L1 - 1.0;
      shape(7) =       L2;
      shape(8) =       L3;
   }
   else if (L2 >= 1.0)   // T2
   {
      shape(5) =       L0;
      shape(7) =       L1;
      shape(2) = L2 - 1.0;
      shape(9) =       L3;
   }
   else if (L3 >= 1.0)   // T3
   {
      shape(6) =       L0;
      shape(8) =       L1;
      shape(9) =       L2;
      shape(3) = L3 - 1.0;
   }
   else if ((L4 <= 1.0) && (L5 <= 1.0))   // T4
   {
      shape(4) = 1.0 - L5;
      shape(5) =       L2;
      shape(6) = 1.0 - L4;
      shape(8) = 1.0 - L0;
   }
   else if ((L4 >= 1.0) && (L5 <= 1.0))   // T5
   {
      shape(4) = 1.0 - L5;
      shape(5) = 1.0 - L1;
      shape(7) = L4 - 1.0;
      shape(8) =       L3;
   }
   else if ((L4 <= 1.0) && (L5 >= 1.0))   // T6
   {
      shape(5) = 1.0 - L3;
      shape(6) = 1.0 - L4;
      shape(8) =       L1;
      shape(9) = L5 - 1.0;
   }
   else if ((L4 >= 1.0) && (L5 >= 1.0))   // T7
   {
      shape(5) =       L0;
      shape(7) = L4 - 1.0;
      shape(8) = 1.0 - L2;
      shape(9) = L5 - 1.0;
   }
}

void RefinedLinear3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   int i,j;

   double L0, L1, L2, L3, L4, L5;
   L0 = 2.0 * ( 1. - ip.x - ip.y - ip.z );
   L1 = 2.0 * ( ip.x );
   L2 = 2.0 * ( ip.y );
   L3 = 2.0 * ( ip.z );
   L4 = 2.0 * ( ip.x + ip.y );
   L5 = 2.0 * ( ip.y + ip.z );

   double DL0[3], DL1[3], DL2[3], DL3[3], DL4[3], DL5[3];
   DL0[0] = -2.0; DL0[1] = -2.0; DL0[2] = -2.0;
   DL1[0] =  2.0; DL1[1] =  0.0; DL1[2] =  0.0;
   DL2[0] =  0.0; DL2[1] =  2.0; DL2[2] =  0.0;
   DL3[0] =  0.0; DL3[1] =  0.0; DL3[2] =  2.0;
   DL4[0] =  2.0; DL4[1] =  2.0; DL4[2] =  0.0;
   DL5[0] =  0.0; DL5[1] =  2.0; DL5[2] =  2.0;

   for (i = 0; i < 10; i++)
      for (j = 0; j < 3; j++)
      {
         dshape(i,j) = 0.0;
      }

   if (L0 >= 1.0)   // T0
   {
      for (j = 0; j < 3; j++)
      {
         dshape(0,j) = DL0[j];
         dshape(4,j) = DL1[j];
         dshape(5,j) = DL2[j];
         dshape(6,j) = DL3[j];
      }
   }
   else if (L1 >= 1.0)   // T1
   {
      for (j = 0; j < 3; j++)
      {
         dshape(4,j) = DL0[j];
         dshape(1,j) = DL1[j];
         dshape(7,j) = DL2[j];
         dshape(8,j) = DL3[j];
      }
   }
   else if (L2 >= 1.0)   // T2
   {
      for (j = 0; j < 3; j++)
      {
         dshape(5,j) = DL0[j];
         dshape(7,j) = DL1[j];
         dshape(2,j) = DL2[j];
         dshape(9,j) = DL3[j];
      }
   }
   else if (L3 >= 1.0)   // T3
   {
      for (j = 0; j < 3; j++)
      {
         dshape(6,j) = DL0[j];
         dshape(8,j) = DL1[j];
         dshape(9,j) = DL2[j];
         dshape(3,j) = DL3[j];
      }
   }
   else if ((L4 <= 1.0) && (L5 <= 1.0))   // T4
   {
      for (j = 0; j < 3; j++)
      {
         dshape(4,j) = - DL5[j];
         dshape(5,j) =   DL2[j];
         dshape(6,j) = - DL4[j];
         dshape(8,j) = - DL0[j];
      }
   }
   else if ((L4 >= 1.0) && (L5 <= 1.0))   // T5
   {
      for (j = 0; j < 3; j++)
      {
         dshape(4,j) = - DL5[j];
         dshape(5,j) = - DL1[j];
         dshape(7,j) =   DL4[j];
         dshape(8,j) =   DL3[j];
      }
   }
   else if ((L4 <= 1.0) && (L5 >= 1.0))   // T6
   {
      for (j = 0; j < 3; j++)
      {
         dshape(5,j) = - DL3[j];
         dshape(6,j) = - DL4[j];
         dshape(8,j) =   DL1[j];
         dshape(9,j) =   DL5[j];
      }
   }
   else if ((L4 >= 1.0) && (L5 >= 1.0))   // T7
   {
      for (j = 0; j < 3; j++)
      {
         dshape(5,j) =   DL0[j];
         dshape(7,j) =   DL4[j];
         dshape(8,j) = - DL2[j];
         dshape(9,j) =   DL5[j];
      }
   }
}


RefinedBiLinear2DFiniteElement::RefinedBiLinear2DFiniteElement()
   : NodalFiniteElement(2, Geometry::SQUARE, 9, 1, FunctionSpace::rQk)
{
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(8).x = 0.5;
   Nodes.IntPoint(8).y = 0.5;
}

void RefinedBiLinear2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                               Vector &shape) const
{
   int i;
   double x = ip.x, y = ip.y;
   double Lx, Ly;
   Lx = 2.0 * ( 1. - x );
   Ly = 2.0 * ( 1. - y );

   // The reference square is split in 4 squares as follows:
   //
   // T0 - 0,4,7,8
   // T1 - 1,4,5,8
   // T2 - 2,5,6,8
   // T3 - 3,6,7,8

   for (i = 0; i < 9; i++)
   {
      shape(i) = 0.0;
   }

   if ((x <= 0.5) && (y <= 0.5))   // T0
   {
      shape(0) = (Lx - 1.0) * (Ly - 1.0);
      shape(4) = (2.0 - Lx) * (Ly - 1.0);
      shape(8) = (2.0 - Lx) * (2.0 - Ly);
      shape(7) = (Lx - 1.0) * (2.0 - Ly);
   }
   else if ((x >= 0.5) && (y <= 0.5))   // T1
   {
      shape(4) =        Lx  * (Ly - 1.0);
      shape(1) = (1.0 - Lx) * (Ly - 1.0);
      shape(5) = (1.0 - Lx) * (2.0 - Ly);
      shape(8) =        Lx  * (2.0 - Ly);
   }
   else if ((x >= 0.5) && (y >= 0.5))   // T2
   {
      shape(8) =        Lx  *        Ly ;
      shape(5) = (1.0 - Lx) *        Ly ;
      shape(2) = (1.0 - Lx) * (1.0 - Ly);
      shape(6) =        Lx  * (1.0 - Ly);
   }
   else if ((x <= 0.5) && (y >= 0.5))   // T3
   {
      shape(7) = (Lx - 1.0) *        Ly ;
      shape(8) = (2.0 - Lx) *        Ly ;
      shape(6) = (2.0 - Lx) * (1.0 - Ly);
      shape(3) = (Lx - 1.0) * (1.0 - Ly);
   }
}

void RefinedBiLinear2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                                DenseMatrix &dshape) const
{
   int i,j;
   double x = ip.x, y = ip.y;
   double Lx, Ly;
   Lx = 2.0 * ( 1. - x );
   Ly = 2.0 * ( 1. - y );

   for (i = 0; i < 9; i++)
      for (j = 0; j < 2; j++)
      {
         dshape(i,j) = 0.0;
      }

   if ((x <= 0.5) && (y <= 0.5))   // T0
   {
      dshape(0,0) =  2.0 * (1.0 - Ly);
      dshape(0,1) =  2.0 * (1.0 - Lx);

      dshape(4,0) =  2.0 * (Ly - 1.0);
      dshape(4,1) = -2.0 * (2.0 - Lx);

      dshape(8,0) =  2.0 * (2.0 - Ly);
      dshape(8,1) =  2.0 * (2.0 - Lx);

      dshape(7,0) = -2.0 * (2.0 - Ly);
      dshape(7,0) =  2.0 * (Lx - 1.0);
   }
   else if ((x >= 0.5) && (y <= 0.5))   // T1
   {
      dshape(4,0) = -2.0 * (Ly - 1.0);
      dshape(4,1) = -2.0 * Lx;

      dshape(1,0) =  2.0 * (Ly - 1.0);
      dshape(1,1) = -2.0 * (1.0 - Lx);

      dshape(5,0) =  2.0 * (2.0 - Ly);
      dshape(5,1) =  2.0 * (1.0 - Lx);

      dshape(8,0) = -2.0 * (2.0 - Ly);
      dshape(8,1) =  2.0 * Lx;
   }
   else if ((x >= 0.5) && (y >= 0.5))   // T2
   {
      dshape(8,0) = -2.0 * Ly;
      dshape(8,1) = -2.0 * Lx;

      dshape(5,0) =  2.0 * Ly;
      dshape(5,1) = -2.0 * (1.0 - Lx);

      dshape(2,0) =  2.0 * (1.0 - Ly);
      dshape(2,1) =  2.0 * (1.0 - Lx);

      dshape(6,0) = -2.0 * (1.0 - Ly);
      dshape(6,1) =  2.0 * Lx;
   }
   else if ((x <= 0.5) && (y >= 0.5))   // T3
   {
      dshape(7,0) = -2.0 * Ly;
      dshape(7,1) = -2.0 * (Lx - 1.0);

      dshape(8,0) =  2.0 * Ly ;
      dshape(8,1) = -2.0 * (2.0 - Lx);

      dshape(6,0) = 2.0 * (1.0 - Ly);
      dshape(6,1) = 2.0 * (2.0 - Lx);

      dshape(3,0) = -2.0 * (1.0 - Ly);
      dshape(3,1) =  2.0 * (Lx - 1.0);
   }
}

RefinedTriLinear3DFiniteElement::RefinedTriLinear3DFiniteElement()
   : NodalFiniteElement(3, Geometry::CUBE, 27, 2, FunctionSpace::rQk)
{
   double I[27];
   double J[27];
   double K[27];
   // nodes
   I[ 0] = 0.0; J[ 0] = 0.0; K[ 0] = 0.0;
   I[ 1] = 1.0; J[ 1] = 0.0; K[ 1] = 0.0;
   I[ 2] = 1.0; J[ 2] = 1.0; K[ 2] = 0.0;
   I[ 3] = 0.0; J[ 3] = 1.0; K[ 3] = 0.0;
   I[ 4] = 0.0; J[ 4] = 0.0; K[ 4] = 1.0;
   I[ 5] = 1.0; J[ 5] = 0.0; K[ 5] = 1.0;
   I[ 6] = 1.0; J[ 6] = 1.0; K[ 6] = 1.0;
   I[ 7] = 0.0; J[ 7] = 1.0; K[ 7] = 1.0;
   // edges
   I[ 8] = 0.5; J[ 8] = 0.0; K[ 8] = 0.0;
   I[ 9] = 1.0; J[ 9] = 0.5; K[ 9] = 0.0;
   I[10] = 0.5; J[10] = 1.0; K[10] = 0.0;
   I[11] = 0.0; J[11] = 0.5; K[11] = 0.0;
   I[12] = 0.5; J[12] = 0.0; K[12] = 1.0;
   I[13] = 1.0; J[13] = 0.5; K[13] = 1.0;
   I[14] = 0.5; J[14] = 1.0; K[14] = 1.0;
   I[15] = 0.0; J[15] = 0.5; K[15] = 1.0;
   I[16] = 0.0; J[16] = 0.0; K[16] = 0.5;
   I[17] = 1.0; J[17] = 0.0; K[17] = 0.5;
   I[18] = 1.0; J[18] = 1.0; K[18] = 0.5;
   I[19] = 0.0; J[19] = 1.0; K[19] = 0.5;
   // faces
   I[20] = 0.5; J[20] = 0.5; K[20] = 0.0;
   I[21] = 0.5; J[21] = 0.0; K[21] = 0.5;
   I[22] = 1.0; J[22] = 0.5; K[22] = 0.5;
   I[23] = 0.5; J[23] = 1.0; K[23] = 0.5;
   I[24] = 0.0; J[24] = 0.5; K[24] = 0.5;
   I[25] = 0.5; J[25] = 0.5; K[25] = 1.0;
   // element
   I[26] = 0.5; J[26] = 0.5; K[26] = 0.5;

   for (int n = 0; n < 27; n++)
   {
      Nodes.IntPoint(n).x = I[n];
      Nodes.IntPoint(n).y = J[n];
      Nodes.IntPoint(n).z = K[n];
   }
}

void RefinedTriLinear3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                                Vector &shape) const
{
   int i, N[8];
   double Lx, Ly, Lz;
   double x = ip.x, y = ip.y, z = ip.z;

   for (i = 0; i < 27; i++)
   {
      shape(i) = 0.0;
   }

   if ((x <= 0.5) && (y <= 0.5) && (z <= 0.5))   // T0
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] =  0;
      N[1] =  8;
      N[2] = 20;
      N[3] = 11;
      N[4] = 16;
      N[5] = 21;
      N[6] = 26;
      N[7] = 24;
   }
   else if ((x >= 0.5) && (y <= 0.5) && (z <= 0.5))   // T1
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] =  8;
      N[1] =  1;
      N[2] =  9;
      N[3] = 20;
      N[4] = 21;
      N[5] = 17;
      N[6] = 22;
      N[7] = 26;
   }
   else if ((x <= 0.5) && (y >= 0.5) && (z <= 0.5))   // T2
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] = 20;
      N[1] =  9;
      N[2] =  2;
      N[3] = 10;
      N[4] = 26;
      N[5] = 22;
      N[6] = 18;
      N[7] = 23;
   }
   else if ((x >= 0.5) && (y >= 0.5) && (z <= 0.5))   // T3
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] = 11;
      N[1] = 20;
      N[2] = 10;
      N[3] =  3;
      N[4] = 24;
      N[5] = 26;
      N[6] = 23;
      N[7] = 19;
   }
   else if ((x <= 0.5) && (y <= 0.5) && (z >= 0.5))   // T4
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 16;
      N[1] = 21;
      N[2] = 26;
      N[3] = 24;
      N[4] =  4;
      N[5] = 12;
      N[6] = 25;
      N[7] = 15;
   }
   else if ((x >= 0.5) && (y <= 0.5) && (z >= 0.5))   // T5
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 21;
      N[1] = 17;
      N[2] = 22;
      N[3] = 26;
      N[4] = 12;
      N[5] =  5;
      N[6] = 13;
      N[7] = 25;
   }
   else if ((x <= 0.5) && (y >= 0.5) && (z >= 0.5))   // T6
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 26;
      N[1] = 22;
      N[2] = 18;
      N[3] = 23;
      N[4] = 25;
      N[5] = 13;
      N[6] =  6;
      N[7] = 14;
   }
   else   // T7
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 24;
      N[1] = 26;
      N[2] = 23;
      N[3] = 19;
      N[4] = 15;
      N[5] = 25;
      N[6] = 14;
      N[7] =  7;
   }

   shape(N[0]) = Lx       * Ly       * Lz;
   shape(N[1]) = (1 - Lx) * Ly       * Lz;
   shape(N[2]) = (1 - Lx) * (1 - Ly) * Lz;
   shape(N[3]) = Lx       * (1 - Ly) * Lz;
   shape(N[4]) = Lx       * Ly       * (1 - Lz);
   shape(N[5]) = (1 - Lx) * Ly       * (1 - Lz);
   shape(N[6]) = (1 - Lx) * (1 - Ly) * (1 - Lz);
   shape(N[7]) = Lx       * (1 - Ly) * (1 - Lz);
}

void RefinedTriLinear3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                                 DenseMatrix &dshape) const
{
   int i, j, N[8];
   double Lx, Ly, Lz;
   double x = ip.x, y = ip.y, z = ip.z;

   for (i = 0; i < 27; i++)
      for (j = 0; j < 3; j++)
      {
         dshape(i,j) = 0.0;
      }

   if ((x <= 0.5) && (y <= 0.5) && (z <= 0.5))   // T0
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] =  0;
      N[1] =  8;
      N[2] = 20;
      N[3] = 11;
      N[4] = 16;
      N[5] = 21;
      N[6] = 26;
      N[7] = 24;
   }
   else if ((x >= 0.5) && (y <= 0.5) && (z <= 0.5))   // T1
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] =  8;
      N[1] =  1;
      N[2] =  9;
      N[3] = 20;
      N[4] = 21;
      N[5] = 17;
      N[6] = 22;
      N[7] = 26;
   }
   else if ((x <= 0.5) && (y >= 0.5) && (z <= 0.5))   // T2
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] = 20;
      N[1] =  9;
      N[2] =  2;
      N[3] = 10;
      N[4] = 26;
      N[5] = 22;
      N[6] = 18;
      N[7] = 23;
   }
   else if ((x >= 0.5) && (y >= 0.5) && (z <= 0.5))   // T3
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 1.0 - 2.0 * z;

      N[0] = 11;
      N[1] = 20;
      N[2] = 10;
      N[3] =  3;
      N[4] = 24;
      N[5] = 26;
      N[6] = 23;
      N[7] = 19;
   }
   else if ((x <= 0.5) && (y <= 0.5) && (z >= 0.5))   // T4
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 16;
      N[1] = 21;
      N[2] = 26;
      N[3] = 24;
      N[4] =  4;
      N[5] = 12;
      N[6] = 25;
      N[7] = 15;
   }
   else if ((x >= 0.5) && (y <= 0.5) && (z >= 0.5))   // T5
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 1.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 21;
      N[1] = 17;
      N[2] = 22;
      N[3] = 26;
      N[4] = 12;
      N[5] =  5;
      N[6] = 13;
      N[7] = 25;
   }
   else if ((x <= 0.5) && (y >= 0.5) && (z >= 0.5))   // T6
   {
      Lx = 2.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 26;
      N[1] = 22;
      N[2] = 18;
      N[3] = 23;
      N[4] = 25;
      N[5] = 13;
      N[6] =  6;
      N[7] = 14;
   }
   else   // T7
   {
      Lx = 1.0 - 2.0 * x;
      Ly = 2.0 - 2.0 * y;
      Lz = 2.0 - 2.0 * z;

      N[0] = 24;
      N[1] = 26;
      N[2] = 23;
      N[3] = 19;
      N[4] = 15;
      N[5] = 25;
      N[6] = 14;
      N[7] =  7;
   }

   dshape(N[0],0) = -2.0 * Ly       * Lz      ;
   dshape(N[0],1) = -2.0 * Lx       * Lz      ;
   dshape(N[0],2) = -2.0 * Lx       * Ly      ;

   dshape(N[1],0) =  2.0 * Ly       * Lz      ;
   dshape(N[1],1) = -2.0 * (1 - Lx) * Lz      ;
   dshape(N[1],2) = -2.0 * (1 - Lx) * Ly      ;

   dshape(N[2],0) =  2.0 * (1 - Ly) * Lz      ;
   dshape(N[2],1) =  2.0 * (1 - Lx) * Lz      ;
   dshape(N[2],2) = -2.0 * (1 - Lx) * (1 - Ly);

   dshape(N[3],0) = -2.0 * (1 - Ly) * Lz      ;
   dshape(N[3],1) =  2.0 * Lx       * Lz      ;
   dshape(N[3],2) = -2.0 * Lx       * (1 - Ly);

   dshape(N[4],0) = -2.0 * Ly       * (1 - Lz);
   dshape(N[4],1) = -2.0 * Lx       * (1 - Lz);
   dshape(N[4],2) =  2.0 * Lx       * Ly      ;

   dshape(N[5],0) =  2.0 * Ly       * (1 - Lz);
   dshape(N[5],1) = -2.0 * (1 - Lx) * (1 - Lz);
   dshape(N[5],2) =  2.0 * (1 - Lx) * Ly      ;

   dshape(N[6],0) =  2.0 * (1 - Ly) * (1 - Lz);
   dshape(N[6],1) =  2.0 * (1 - Lx) * (1 - Lz);
   dshape(N[6],2) =  2.0 * (1 - Lx) * (1 - Ly);

   dshape(N[7],0) = -2.0 * (1 - Ly) * (1 - Lz);
   dshape(N[7],1) =  2.0 * Lx       * (1 - Lz);
   dshape(N[7],2) =  2.0 * Lx       * (1 - Ly);
}


Nedelec1HexFiniteElement::Nedelec1HexFiniteElement()
   : VectorFiniteElement(3, Geometry::CUBE, 12, 1, H_CURL, FunctionSpace::Qk)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(1).z = 0.0;

   Nodes.IntPoint(2).x = 0.5;
   Nodes.IntPoint(2).y = 1.0;
   Nodes.IntPoint(2).z = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.5;
   Nodes.IntPoint(3).z = 0.0;

   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 1.0;

   Nodes.IntPoint(5).x = 1.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 1.0;

   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 1.0;
   Nodes.IntPoint(6).z = 1.0;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(7).z = 1.0;

   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.5;

   Nodes.IntPoint(9).x = 1.0;
   Nodes.IntPoint(9).y = 0.0;
   Nodes.IntPoint(9).z = 0.5;

   Nodes.IntPoint(10).x= 1.0;
   Nodes.IntPoint(10).y= 1.0;
   Nodes.IntPoint(10).z= 0.5;

   Nodes.IntPoint(11).x= 0.0;
   Nodes.IntPoint(11).y= 1.0;
   Nodes.IntPoint(11).z= 0.5;
}

void Nedelec1HexFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                          DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;

   shape(0,0) = (1. - y) * (1. - z);
   shape(0,1) = 0.;
   shape(0,2) = 0.;

   shape(2,0) = y * (1. - z);
   shape(2,1) = 0.;
   shape(2,2) = 0.;

   shape(4,0) = z * (1. - y);
   shape(4,1) = 0.;
   shape(4,2) = 0.;

   shape(6,0) = y * z;
   shape(6,1) = 0.;
   shape(6,2) = 0.;

   shape(1,0) = 0.;
   shape(1,1) = x * (1. - z);
   shape(1,2) = 0.;

   shape(3,0) = 0.;
   shape(3,1) = (1. - x) * (1. - z);
   shape(3,2) = 0.;

   shape(5,0) = 0.;
   shape(5,1) = x * z;
   shape(5,2) = 0.;

   shape(7,0) = 0.;
   shape(7,1) = (1. - x) * z;
   shape(7,2) = 0.;

   shape(8,0) = 0.;
   shape(8,1) = 0.;
   shape(8,2) = (1. - x) * (1. - y);

   shape(9,0) = 0.;
   shape(9,1) = 0.;
   shape(9,2) = x * (1. - y);

   shape(10,0) = 0.;
   shape(10,1) = 0.;
   shape(10,2) = x * y;

   shape(11,0) = 0.;
   shape(11,1) = 0.;
   shape(11,2) = y * (1. - x);

}

void Nedelec1HexFiniteElement::CalcCurlShape(const IntegrationPoint &ip,
                                             DenseMatrix &curl_shape)
const
{
   double x = ip.x, y = ip.y, z = ip.z;

   curl_shape(0,0) = 0.;
   curl_shape(0,1) = y - 1.;
   curl_shape(0,2) = 1. - z;

   curl_shape(2,0) = 0.;
   curl_shape(2,1) = -y;
   curl_shape(2,2) = z - 1.;

   curl_shape(4,0) = 0;
   curl_shape(4,1) = 1. - y;
   curl_shape(4,2) = z;

   curl_shape(6,0) = 0.;
   curl_shape(6,1) = y;
   curl_shape(6,2) = -z;

   curl_shape(1,0) = x;
   curl_shape(1,1) = 0.;
   curl_shape(1,2) = 1. - z;

   curl_shape(3,0) = 1. - x;
   curl_shape(3,1) = 0.;
   curl_shape(3,2) = z - 1.;

   curl_shape(5,0) = -x;
   curl_shape(5,1) = 0.;
   curl_shape(5,2) = z;

   curl_shape(7,0) = x - 1.;
   curl_shape(7,1) = 0.;
   curl_shape(7,2) = -z;

   curl_shape(8,0) = x - 1.;
   curl_shape(8,1) = 1. - y;
   curl_shape(8,2) = 0.;

   curl_shape(9,0) = -x;
   curl_shape(9,1) = y - 1.;
   curl_shape(9,2) = 0;

   curl_shape(10,0) = x;
   curl_shape(10,1) = -y;
   curl_shape(10,2) = 0.;

   curl_shape(11,0) = 1. - x;
   curl_shape(11,1) = y;
   curl_shape(11,2) = 0.;
}

const double Nedelec1HexFiniteElement::tk[12][3] =
{
   {1,0,0}, {0,1,0}, {1,0,0}, {0,1,0},
   {1,0,0}, {0,1,0}, {1,0,0}, {0,1,0},
   {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}
};

void Nedelec1HexFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 12; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 12; j++)
      {
         double d = ( vshape(j,0)*tk[k][0] + vshape(j,1)*tk[k][1] +
                      vshape(j,2)*tk[k][2] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "Nedelec1HexFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear (more to have embedding?)
   const DenseMatrix &J = Trans.Jacobian();
   double vk[3];
   Vector xk (vk, 3);

   for (k = 0; k < 12; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2];
      CalcVShape (ip, vshape);
      //  vk = J tk
      vk[0] = J(0,0)*tk[k][0]+J(0,1)*tk[k][1]+J(0,2)*tk[k][2];
      vk[1] = J(1,0)*tk[k][0]+J(1,1)*tk[k][1]+J(1,2)*tk[k][2];
      vk[2] = J(2,0)*tk[k][0]+J(2,1)*tk[k][1]+J(2,2)*tk[k][2];
      for (j = 0; j < 12; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void Nedelec1HexFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[3];
   Vector xk (vk, 3);

   for (int k = 0; k < 12; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      const DenseMatrix &J = Trans.Jacobian();

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t J tk
      dofs(k) =
         vk[0] * ( J(0,0)*tk[k][0]+J(0,1)*tk[k][1]+J(0,2)*tk[k][2] ) +
         vk[1] * ( J(1,0)*tk[k][0]+J(1,1)*tk[k][1]+J(1,2)*tk[k][2] ) +
         vk[2] * ( J(2,0)*tk[k][0]+J(2,1)*tk[k][1]+J(2,2)*tk[k][2] );
   }
}


Nedelec1TetFiniteElement::Nedelec1TetFiniteElement()
   : VectorFiniteElement(3, Geometry::TETRAHEDRON, 6, 1, H_CURL)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;

   Nodes.IntPoint(1).x = 0.0;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(1).z = 0.0;

   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.5;

   Nodes.IntPoint(3).x = 0.5;
   Nodes.IntPoint(3).y = 0.5;
   Nodes.IntPoint(3).z = 0.0;

   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.5;

   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 0.5;
}

void Nedelec1TetFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                          DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;

   shape(0,0) = 1. - y - z;
   shape(0,1) = x;
   shape(0,2) = x;

   shape(1,0) = y;
   shape(1,1) = 1. - x - z;
   shape(1,2) = y;

   shape(2,0) = z;
   shape(2,1) = z;
   shape(2,2) = 1. - x - y;

   shape(3,0) = -y;
   shape(3,1) = x;
   shape(3,2) = 0.;

   shape(4,0) = -z;
   shape(4,1) = 0.;
   shape(4,2) = x;

   shape(5,0) = 0.;
   shape(5,1) = -z;
   shape(5,2) = y;
}

void Nedelec1TetFiniteElement::CalcCurlShape(const IntegrationPoint &ip,
                                             DenseMatrix &curl_shape)
const
{
   curl_shape(0,0) =  0.;
   curl_shape(0,1) = -2.;
   curl_shape(0,2) =  2.;

   curl_shape(1,0) =  2.;
   curl_shape(1,1) =  0.;
   curl_shape(1,2) = -2.;

   curl_shape(2,0) = -2.;
   curl_shape(2,1) =  2.;
   curl_shape(2,2) =  0.;

   curl_shape(3,0) = 0.;
   curl_shape(3,1) = 0.;
   curl_shape(3,2) = 2.;

   curl_shape(4,0) =  0.;
   curl_shape(4,1) = -2.;
   curl_shape(4,2) =  0.;

   curl_shape(5,0) = 2.;
   curl_shape(5,1) = 0.;
   curl_shape(5,2) = 0.;
}

const double Nedelec1TetFiniteElement::tk[6][3] =
{{1,0,0}, {0,1,0}, {0,0,1}, {-1,1,0}, {-1,0,1}, {0,-1,1}};

void Nedelec1TetFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 6; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 6; j++)
      {
         double d = ( vshape(j,0)*tk[k][0] + vshape(j,1)*tk[k][1] +
                      vshape(j,2)*tk[k][2] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "Nedelec1TetFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   const DenseMatrix &J = Trans.Jacobian();
   double vk[3];
   Vector xk (vk, 3);

   for (k = 0; k < 6; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2];
      CalcVShape (ip, vshape);
      //  vk = J tk
      vk[0] = J(0,0)*tk[k][0]+J(0,1)*tk[k][1]+J(0,2)*tk[k][2];
      vk[1] = J(1,0)*tk[k][0]+J(1,1)*tk[k][1]+J(1,2)*tk[k][2];
      vk[2] = J(2,0)*tk[k][0]+J(2,1)*tk[k][1]+J(2,2)*tk[k][2];
      for (j = 0; j < 6; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void Nedelec1TetFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[3];
   Vector xk (vk, 3);

   for (int k = 0; k < 6; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      const DenseMatrix &J = Trans.Jacobian();

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t J tk
      dofs(k) =
         vk[0] * ( J(0,0)*tk[k][0]+J(0,1)*tk[k][1]+J(0,2)*tk[k][2] ) +
         vk[1] * ( J(1,0)*tk[k][0]+J(1,1)*tk[k][1]+J(1,2)*tk[k][2] ) +
         vk[2] * ( J(2,0)*tk[k][0]+J(2,1)*tk[k][1]+J(2,2)*tk[k][2] );
   }
}

Nedelec1PentFiniteElement::Nedelec1PentFiniteElement()
   : VectorFiniteElement(4, Geometry::PENTATOPE, 10, 1, H_CURL)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 0.0;
   Nodes.IntPoint(1).y = 0.5;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.5;
   Nodes.IntPoint(2).t = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.0;
   Nodes.IntPoint(3).z = 0.0;
   Nodes.IntPoint(3).t = 0.5;

   Nodes.IntPoint(4).x = 0.5;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(4).t = 0.0;

   Nodes.IntPoint(5).x = 0.5;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 0.5;
   Nodes.IntPoint(5).t = 0.0;

   Nodes.IntPoint(6).x = 0.5;
   Nodes.IntPoint(6).y = 0.0;
   Nodes.IntPoint(6).z = 0.0;
   Nodes.IntPoint(6).t = 0.5;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.5;
   Nodes.IntPoint(7).z = 0.5;
   Nodes.IntPoint(7).t = 0.0;

   Nodes.IntPoint(8).x = 0.0;
   Nodes.IntPoint(8).y = 0.5;
   Nodes.IntPoint(8).z = 0.0;
   Nodes.IntPoint(8).t = 0.5;

   Nodes.IntPoint(9).x = 0.0;
   Nodes.IntPoint(9).y = 0.0;
   Nodes.IntPoint(9).z = 0.5;
   Nodes.IntPoint(9).t = 0.5;
}

void Nedelec1PentFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                           DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z, t = ip.t;

   shape(0,0) = 1. - y - z - t;
   shape(0,1) = x;
   shape(0,2) = x;
   shape(0,3) = x;

   shape(1,0) = y;
   shape(1,1) = 1. - x - z - t;
   shape(1,2) = y;
   shape(1,3) = y;

   shape(2,0) = z;
   shape(2,1) = z;
   shape(2,2) = 1. - x - y - t;
   shape(2,3) = z;

   shape(3,0) = t;
   shape(3,1) = t;
   shape(3,2) = t;
   shape(3,3) = 1. - x - y - z;

   shape(4,0) = -y;
   shape(4,1) = x;
   shape(4,2) = 0.;
   shape(4,3) = 0.;

   shape(5,0) = -z;
   shape(5,1) = 0.;
   shape(5,2) = x;
   shape(5,3) = 0.;

   shape(6,0) = -t;
   shape(6,1) = 0.;
   shape(6,2) = 0;
   shape(6,3) = x;

   shape(7,0) = 0.;
   shape(7,1) = -z;
   shape(7,2) = y;
   shape(7,3) = 0.;

   shape(8,0) = 0.;
   shape(8,1) = -t;
   shape(8,2) = 0.;
   shape(8,3) = y;

   shape(9,0) = 0.;
   shape(9,1) = 0.;
   shape(9,2) = -t;
   shape(9,3) = z;
}

void Nedelec1PentFiniteElement::CalcCurlShape(const IntegrationPoint &ip,
                                              DenseMatrix &curl_shape)
const
{
   curl_shape(0,0) =  2.;
   curl_shape(0,1) =  2.;
   curl_shape(0,2) =  2.;
   curl_shape(0,3) =  0.;
   curl_shape(0,4) =  0.;
   curl_shape(0,5) =  0.;

   curl_shape(1,0) =  -2.;
   curl_shape(1,1) =  0.;
   curl_shape(1,2) =  0.;
   curl_shape(1,3) =  2.;
   curl_shape(1,4) =  2.;
   curl_shape(1,5) =  0.;

   curl_shape(2,0) =  0.;
   curl_shape(2,1) =  -2.;
   curl_shape(2,2) =  0.;
   curl_shape(2,3) =  -2.;
   curl_shape(2,4) =  0.;
   curl_shape(2,5) =  2.;

   curl_shape(3,0) =  0.;
   curl_shape(3,1) =  0.;
   curl_shape(3,2) =  -2.;
   curl_shape(3,3) =  0.;
   curl_shape(3,4) =  -2.;
   curl_shape(3,5) =  -2.;

   curl_shape(4,0) =  2.;
   curl_shape(4,1) =  0.;
   curl_shape(4,2) =  0.;
   curl_shape(4,3) =  0.;
   curl_shape(4,4) =  0.;
   curl_shape(4,5) =  0.;

   curl_shape(5,0) =  0.;
   curl_shape(5,1) =  2.;
   curl_shape(5,2) =  0.;
   curl_shape(5,3) =  0.;
   curl_shape(5,4) =  0.;
   curl_shape(5,5) =  0.;

   curl_shape(6,0) =  0.;
   curl_shape(6,1) =  0.;
   curl_shape(6,2) =  2.;
   curl_shape(6,3) =  0.;
   curl_shape(6,4) =  0.;
   curl_shape(6,5) =  0.;

   curl_shape(7,0) =  0.;
   curl_shape(7,1) =  0.;
   curl_shape(7,2) =  0.;
   curl_shape(7,3) =  2.;
   curl_shape(7,4) =  0.;
   curl_shape(7,5) =  0.;

   curl_shape(8,0) =  0.;
   curl_shape(8,1) =  0.;
   curl_shape(8,2) =  0.;
   curl_shape(8,3) =  0.;
   curl_shape(8,4) =  2.;
   curl_shape(8,5) =  0.;

   curl_shape(9,0) =  0.;
   curl_shape(9,1) =  0.;
   curl_shape(9,2) =  0.;
   curl_shape(9,3) =  0.;
   curl_shape(9,4) =  0.;
   curl_shape(9,5) =  2.;

}

const double Nedelec1PentFiniteElement::tk[10][4] =
{{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {-1,1,0,0}, {-1,0,1,0}, {-1,0,0,1}, {0,-1,1,0}, {0,-1,0,1}, {0,0,-1,1}};

void Nedelec1PentFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(Dof, Dim);
#endif

#ifdef MFEM_DEBUG
   cout << "not implemented in NED_FE_PENT" << endl;
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = ip.t = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   const DenseMatrix &J = Trans.Jacobian();
   double vk[4];
   Vector xk (vk, 4);

   for (k = 0; k < 10; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2]; ip.t = vk[3];
      CalcVShape (ip, vshape);
      //  vk = J tk
      vk[0] = J(0,0)*tk[k][0]+J(0,1)*tk[k][1]+J(0,2)*tk[k][2]+J(0,3)*tk[k][3];
      vk[1] = J(1,0)*tk[k][0]+J(1,1)*tk[k][1]+J(1,2)*tk[k][2]+J(1,3)*tk[k][3];
      vk[2] = J(2,0)*tk[k][0]+J(2,1)*tk[k][1]+J(2,2)*tk[k][2]+J(2,3)*tk[k][3];
      vk[3] = J(3,0)*tk[k][0]+J(3,1)*tk[k][1]+J(3,2)*tk[k][2]+J(3,3)*tk[k][3];
      for (j = 0; j < 10; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2]+vshape(j,3)*vk[3])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void Nedelec1PentFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   dofs.SetSize(10); dofs = 0.0;
   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      const DenseMatrix &J = Trans.Jacobian();

      vc.Eval(vi, Trans, Nodes.IntPoint (k));

      J.Mult(tk[k],t);

      dofs(k) = vi * ti;
   }
}

void Nedelec1PentFiniteElement::Project(const FiniteElement &fe,
                                        ElementTransformation &Trans,
                                        DenseMatrix &Id) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   Id.SetSize(10,dim*dof); Id = 0.0;

   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   Vector shape(dof);

   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      fe.CalcShape(Nodes.IntPoint(k), shape);

      ti = 0.0;
      J.Mult(tk[k],t);

      for (int di=0; di<dim; di++)
      {
         for (int j=0; j<dof; j++)
         {
            vi = 0.0; v[di] = shape(j);

            Id(k, di*dof+j) = vi * ti;
         }
      }
   }
}

void Nedelec1PentFiniteElement::ProjectGrad(const FiniteElement &fe,
                                            ElementTransformation &Trans,
                                            DenseMatrix &grad) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   grad.SetSize(10,dof); grad = 0.0;

   DenseMatrix dshape(dof,dim);
   DenseMatrix gshape(dof,dim);
   DenseMatrix invJ(dim,dim);

   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      CalcInverse(J, invJ);

      fe.CalcDShape(Nodes.IntPoint(k), dshape);
      Mult(dshape, invJ, gshape);

      ti = 0.0;
      J.Mult(tk[k],t);

      for (int j=0; j<dof; j++)
      {
         v[0] = gshape(j,0);
         v[1] = gshape(j,1);
         v[2] = gshape(j,2);
         v[3] = gshape(j,3);

         grad(k, j) = vi * ti;
      }

   }
}

Nedelec1FullPentFiniteElement::Nedelec1FullPentFiniteElement()
   : VectorFiniteElement(4, Geometry::PENTATOPE, 20, 1, H_CURL)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(0).y = 0.0;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 1.0;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 0.0;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(2).t = 0.0;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(3).z = 0.0;
   Nodes.IntPoint(3).t = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.0;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(4).t = 0.0;

   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 1.0;
   Nodes.IntPoint(5).t = 0.0;

   Nodes.IntPoint(6).x = 0.0;
   Nodes.IntPoint(6).y = 0.0;
   Nodes.IntPoint(6).z = 0.0;
   Nodes.IntPoint(6).t = 0.0;

   Nodes.IntPoint(7).x = 0.0;
   Nodes.IntPoint(7).y = 0.0;
   Nodes.IntPoint(7).z = 0.0;
   Nodes.IntPoint(7).t = 1.0;

   Nodes.IntPoint(8).x = 1.0;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.0;
   Nodes.IntPoint(8).t = 0.0;

   Nodes.IntPoint(9).x = 0.0;
   Nodes.IntPoint(9).y = 1.0;
   Nodes.IntPoint(9).z = 0.0;
   Nodes.IntPoint(9).t = 0.0;

   Nodes.IntPoint(10).x = 1.0;
   Nodes.IntPoint(10).y = 0.0;
   Nodes.IntPoint(10).z = 0.0;
   Nodes.IntPoint(10).t = 0.0;

   Nodes.IntPoint(11).x = 0.0;
   Nodes.IntPoint(11).y = 0.0;
   Nodes.IntPoint(11).z = 1.0;
   Nodes.IntPoint(11).t = 0.0;

   Nodes.IntPoint(12).x = 1.0;
   Nodes.IntPoint(12).y = 0.0;
   Nodes.IntPoint(12).z = 0.0;
   Nodes.IntPoint(12).t = 0.0;

   Nodes.IntPoint(13).x = 0.0;
   Nodes.IntPoint(13).y = 0.0;
   Nodes.IntPoint(13).z = 0.0;
   Nodes.IntPoint(13).t = 1.0;

   Nodes.IntPoint(14).x = 0.0;
   Nodes.IntPoint(14).y = 1.0;
   Nodes.IntPoint(14).z = 0.0;
   Nodes.IntPoint(14).t = 0.0;

   Nodes.IntPoint(15).x = 0.0;
   Nodes.IntPoint(15).y = 0.0;
   Nodes.IntPoint(15).z = 1.0;
   Nodes.IntPoint(15).t = 0.0;

   Nodes.IntPoint(16).x = 0.0;
   Nodes.IntPoint(16).y = 1.0;
   Nodes.IntPoint(16).z = 0.0;
   Nodes.IntPoint(16).t = 0.0;

   Nodes.IntPoint(17).x = 0.0;
   Nodes.IntPoint(17).y = 0.0;
   Nodes.IntPoint(17).z = 0.0;
   Nodes.IntPoint(17).t = 1.0;

   Nodes.IntPoint(18).x = 0.0;
   Nodes.IntPoint(18).y = 0.0;
   Nodes.IntPoint(18).z = 1.0;
   Nodes.IntPoint(18).t = 0.0;

   Nodes.IntPoint(19).x = 0.0;
   Nodes.IntPoint(19).y = 0.0;
   Nodes.IntPoint(19).z = 0.0;
   Nodes.IntPoint(19).t = 1.0;
}

const double Nedelec1FullPentFiniteElement::tk[10][4] =
{{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {-1,1,0,0}, {-1,0,1,0}, {-1,0,0,1}, {0,-1,1,0}, {0,-1,0,1}, {0,0,-1,1}};

void Nedelec1FullPentFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                               DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z, t = ip.t;

   shape(0,0) = 1. -x - y - z - t;
   shape(0,1) = 0;
   shape(0,2) = 0;
   shape(0,3) = 0;

   shape(1,0) = x;
   shape(1,1) = x;
   shape(1,2) = x;
   shape(1,3) = x;



   shape(2,0) = 0;
   shape(2,1) = 1. - x - y - z - t;
   shape(2,2) = 0;
   shape(2,3) = 0;

   shape(3,0) = y;
   shape(3,1) = y;
   shape(3,2) = y;
   shape(3,3) = y;



   shape(4,0) = 0;
   shape(4,1) = 0;
   shape(4,2) = 1. - x - y - z - t;
   shape(4,3) = 0;

   shape(5,0) = z;
   shape(5,1) = z;
   shape(5,2) = z;
   shape(5,3) = z;



   shape(6,0) = 0;
   shape(6,1) = 0;
   shape(6,2) = 0;
   shape(6,3) = 1. - x - y - z - t;

   shape(7,0) = t;
   shape(7,1) = t;
   shape(7,2) = t;
   shape(7,3) = t;



   shape(8,0) = 0;
   shape(8,1) = x;
   shape(8,2) = 0.;
   shape(8,3) = 0.;

   shape(9,0) = -y;
   shape(9,1) = 0;
   shape(9,2) = 0;
   shape(9,3) = 0;



   shape(10,0) = 0;
   shape(10,1) = 0.;
   shape(10,2) = x;
   shape(10,3) = 0.;

   shape(11,0) = -z;
   shape(11,1) = 0;
   shape(11,2) = 0;
   shape(11,3) = 0;



   shape(12,0) = 0;
   shape(12,1) = 0.;
   shape(12,2) = 0;
   shape(12,3) = x;

   shape(13,0) = -t;
   shape(13,1) = 0;
   shape(13,2) = 0;
   shape(13,3) = 0;



   shape(14,0) = 0.;
   shape(14,1) = 0;
   shape(14,2) = y;
   shape(14,3) = 0.;

   shape(15,0) = 0;
   shape(15,1) = -z;
   shape(15,2) = 0;
   shape(15,3) = 0;



   shape(16,0) = 0.;
   shape(16,1) = 0;
   shape(16,2) = 0.;
   shape(16,3) = y;

   shape(17,0) = 0;
   shape(17,1) = -t;
   shape(17,2) = 0;
   shape(17,3) = 0;



   shape(18,0) = 0.;
   shape(18,1) = 0.;
   shape(18,2) = 0;
   shape(18,3) = z;

   shape(19,0) = 0;
   shape(19,1) = 0;
   shape(19,2) = -t;
   shape(19,3) = 0;
}

void Nedelec1FullPentFiniteElement::CalcCurlShape(const IntegrationPoint &ip,
                                                  DenseMatrix &curl_shape)
const
{
   curl_shape(0,0) = 1.0;
   curl_shape(0,1) = 1.0;
   curl_shape(0,2) = 1.0;
   curl_shape(0,3) = 0.;
   curl_shape(0,4) = 0.;
   curl_shape(0,5) = 0.;

   curl_shape(1,0) = 1.0;
   curl_shape(1,1) = 1.0;
   curl_shape(1,2) = 1.0;
   curl_shape(1,3) = 0.;
   curl_shape(1,4) = 0.;
   curl_shape(1,5) = 0.;

   curl_shape(2,0) = -1.0;
   curl_shape(2,1) = 0.;
   curl_shape(2,2) = 0.;
   curl_shape(2,3) = 1.0;
   curl_shape(2,4) = 1.0;
   curl_shape(2,5) = 0.;

   curl_shape(3,0) = -1.0;
   curl_shape(3,1) = 0.;
   curl_shape(3,2) = 0.;
   curl_shape(3,3) = 1.0;
   curl_shape(3,4) = 1.0;
   curl_shape(3,5) = 0.;

   curl_shape(4,0) = 0.;
   curl_shape(4,1) = -1.0;
   curl_shape(4,2) = 0.;
   curl_shape(4,3) = -1.0;
   curl_shape(4,4) = 0.;
   curl_shape(4,5) = 1.0;

   curl_shape(5,0) = 0.;
   curl_shape(5,1) = -1.0;
   curl_shape(5,2) = 0.;
   curl_shape(5,3) = -1.0;
   curl_shape(5,4) = 0.;
   curl_shape(5,5) = 1.0;

   curl_shape(6,0) = 0.;
   curl_shape(6,1) = 0.;
   curl_shape(6,2) = -1.0;
   curl_shape(6,3) = 0.;
   curl_shape(6,4) = -1.0;
   curl_shape(6,5) = -1.0;

   curl_shape(7,0) = 0.;
   curl_shape(7,1) = 0.;
   curl_shape(7,2) = -1.0;
   curl_shape(7,3) = 0.;
   curl_shape(7,4) = -1.0;
   curl_shape(7,5) = -1.0;

   curl_shape(8,0) = 1.0;
   curl_shape(8,1) = 0.;
   curl_shape(8,2) = 0.;
   curl_shape(8,3) = 0.;
   curl_shape(8,4) = 0.;
   curl_shape(8,5) = 0.;

   curl_shape(9,0) = 1.0;
   curl_shape(9,1) = 0.;
   curl_shape(9,2) = 0.;
   curl_shape(9,3) = 0.;
   curl_shape(9,4) = 0.;
   curl_shape(9,5) = 0.;

   curl_shape(10,0) = 0.;
   curl_shape(10,1) = 1.0;
   curl_shape(10,2) = 0.;
   curl_shape(10,3) = 0.;
   curl_shape(10,4) = 0.;
   curl_shape(10,5) = 0.;

   curl_shape(11,0) = 0.;
   curl_shape(11,1) = 1.0;
   curl_shape(11,2) = 0.;
   curl_shape(11,3) = 0.;
   curl_shape(11,4) = 0.;
   curl_shape(11,5) = 0.;

   curl_shape(12,0) = 0.;
   curl_shape(12,1) = 0.;
   curl_shape(12,2) = 1.0;
   curl_shape(12,3) = 0.;
   curl_shape(12,4) = 0.;
   curl_shape(12,5) = 0.;

   curl_shape(13,0) = 0.;
   curl_shape(13,1) = 0.;
   curl_shape(13,2) = 1.0;
   curl_shape(13,3) = 0.;
   curl_shape(13,4) = 0.;
   curl_shape(13,5) = 0.;

   curl_shape(14,0) = 0.;
   curl_shape(14,1) = 0.;
   curl_shape(14,2) = 0.;
   curl_shape(14,3) = 1.0;
   curl_shape(14,4) = 0.;
   curl_shape(14,5) = 0.;

   curl_shape(15,0) = 0.;
   curl_shape(15,1) = 0.;
   curl_shape(15,2) = 0.;
   curl_shape(15,3) = 1.0;
   curl_shape(15,4) = 0.;
   curl_shape(15,5) = 0.;

   curl_shape(16,0) = 0.;
   curl_shape(16,1) = 0.;
   curl_shape(16,2) = 0.;
   curl_shape(16,3) = 0.;
   curl_shape(16,4) = 1.0;
   curl_shape(16,5) = 0.;

   curl_shape(17,0) = 0.;
   curl_shape(17,1) = 0.;
   curl_shape(17,2) = 0.;
   curl_shape(17,3) = 0.;
   curl_shape(17,4) = 1.0;
   curl_shape(17,5) = 0.;

   curl_shape(18,0) = 0.;
   curl_shape(18,1) = 0.;
   curl_shape(18,2) = 0.;
   curl_shape(18,3) = 0.;
   curl_shape(18,4) = 0.;
   curl_shape(18,5) = 1.0;

   curl_shape(19,0) = 0.;
   curl_shape(19,1) = 0.;
   curl_shape(19,2) = 0.;
   curl_shape(19,3) = 0.;
   curl_shape(19,4) = 0.;
   curl_shape(19,5) = 1.0;

}

void Nedelec1FullPentFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   dofs.SetSize(20); dofs = 0.0;
   for (int k = 0; k < 10; k++)
   {
      for (int ld=0; ld<2; ld++)
      {
         int dofID=2*k+ld;
         Trans.SetIntPoint(&Nodes.IntPoint(dofID));
         const DenseMatrix &J = Trans.Jacobian();
         vc.Eval(vi, Trans, Nodes.IntPoint(dofID));
         J.Mult(tk[k],t);
         dofs(dofID) = vi * ti;
      }
   }
}

void Nedelec1FullPentFiniteElement::Project(const FiniteElement &fe,
                                            ElementTransformation &Trans,
                                            DenseMatrix &Id) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   Id.SetSize(20,dim*dof); Id = 0.0;

   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   Vector shape(dof);

   for (int k = 0; k < 10; k++)
   {
      for (int ld=0; ld<2; ld++)
      {
         int dofID=2*k+ld;

         Trans.SetIntPoint(&Nodes.IntPoint(dofID));
         const DenseMatrix &J = Trans.Jacobian();
         fe.CalcShape(Nodes.IntPoint(dofID), shape);

         J.Mult(tk[k],t);

         for (int di=0; di<dim; di++)
         {
            for (int j=0; j<dof; j++)
            {
               vi = 0.0; v[di] = shape(j);

               Id(dofID, di*dof+j) = vi * ti;
            }
         }
      }
   }
}

void Nedelec1FullPentFiniteElement::ProjectGrad(const FiniteElement &fe,
                                                ElementTransformation &Trans,
                                                DenseMatrix &grad) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   grad.SetSize(20,dof); grad = 0.0;

   DenseMatrix dshape(dof,dim);
   DenseMatrix gshape(dof,dim);
   DenseMatrix invJ(dim,dim);

   double v[4]; Vector vi(v, 4);
   double t[4]; Vector ti(t, 4);

   for (int k = 0; k < 10; k++)
   {
      for (int ld=0; ld<2; ld++)
      {
         int dofID=2*k+ld;
         Trans.SetIntPoint(&Nodes.IntPoint(dofID));
         const DenseMatrix &J = Trans.Jacobian();
         CalcInverse(J, invJ);

         fe.CalcDShape(Nodes.IntPoint(dofID), dshape);
         Mult(dshape, invJ, gshape);

         J.Mult(tk[k],t);

         for (int j=0; j<dof; j++)
         {
            v[0] = gshape(j,0);
            v[1] = gshape(j,1);
            v[2] = gshape(j,2);
            v[3] = gshape(j,3);

            grad(dofID, j) = vi * ti;
         }
      }

   }
}

DivSkew1PentFiniteElement::DivSkew1PentFiniteElement()
   : VectorFiniteElement(4, Geometry::PENTATOPE, 10, 1, H_DIV_SKEW)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(0).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(0).z = 0.0;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(1).t = 0.0;

   Nodes.IntPoint(2).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.0;
   Nodes.IntPoint(2).t = 0.3333333333333333333333333333333333;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(3).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(3).t = 0.0;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(4).z = 0.0;
   Nodes.IntPoint(4).t = 0.3333333333333333333333333333333333;

   Nodes.IntPoint(5).x = 0.0;
   Nodes.IntPoint(5).y = 0.0;
   Nodes.IntPoint(5).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(5).t = 0.3333333333333333333333333333333333;

   Nodes.IntPoint(6).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(6).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(6).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(6).t = 0.0;

   Nodes.IntPoint(7).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(7).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(7).z = 0.0;
   Nodes.IntPoint(7).t = 0.3333333333333333333333333333333333;

   Nodes.IntPoint(8).x = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(8).y = 0.0;
   Nodes.IntPoint(8).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(8).t = 0.3333333333333333333333333333333333;

   Nodes.IntPoint(9).x = 0.0;
   Nodes.IntPoint(9).y = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(9).z = 0.3333333333333333333333333333333333;
   Nodes.IntPoint(9).t = 0.3333333333333333333333333333333333;
}

const double DivSkew1PentFiniteElement::tk1[10][4] =
{{1,0,0,0},{1,0,0,0},{1,0,0,0},{0,1,0,0},{0,1,0,0},{0,0,1,0},{-1,1,0,0},{-1,1,0,0},{-1,0,1,0},{0,-1,1,0}};

const double DivSkew1PentFiniteElement::tk2[10][4] =
{{-1,1,0,0},{-1,0,1,0},{-1,0,0,1},{0,-1,1,0},{0,-1,0,1},{0,0,-1,1},{0,-1,1,0},{0,-1,0,1},{0,0,-1,1},{0,0,-1,1}};

void DivSkew1PentFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                           DenseMatrix &shape) const
{
   double x1 = ip.x, x2 = ip.y, x3 = ip.z, x4 = ip.t;

   shape(0,0) = 0.;
   shape(0,1) = 0.;
   shape(0,2) = x1;
   shape(0,3) = -1.*x1;
   shape(0,4) = 0.;
   shape(0,5) = 0.;
   shape(0,6) = x2;
   shape(0,7) = -1.*x2;
   shape(0,8) = -1.*x1;
   shape(0,9) = -1.*x2;
   shape(0,10) = 0.;
   shape(0,11) = 1. - 1.*x3 - 1.*x4;
   shape(0,12) = x1;
   shape(0,13) = x2;
   shape(0,14) = -1. + x3 + x4;
   shape(0,15) = 0.;

   shape(1,0) = 0.;
   shape(1,1) = -1.*x1;
   shape(1,2) = 0.;
   shape(1,3) = x1;
   shape(1,4) = x1;
   shape(1,5) = 0.;
   shape(1,6) = x3;
   shape(1,7) = -1. + x2 + x4;
   shape(1,8) = 0.;
   shape(1,9) = -1.*x3;
   shape(1,10) = 0.;
   shape(1,11) = x3;
   shape(1,12) = -1.*x1;
   shape(1,13) = 1. - 1.*x2 - 1.*x4;
   shape(1,14) = -1.*x3;
   shape(1,15) = 0.;

   shape(2,0) = 0.;
   shape(2,1) = x1;
   shape(2,2) = -1.*x1;
   shape(2,3) = 0.;
   shape(2,4) = -1.*x1;
   shape(2,5) = 0.;
   shape(2,6) = 1. - 1.*x2 - 1.*x3;
   shape(2,7) = -1.*x4;
   shape(2,8) = x1;
   shape(2,9) = -1. + x2 + x3;
   shape(2,10) = 0.;
   shape(2,11) = x4;
   shape(2,12) = 0.;
   shape(2,13) = x4;
   shape(2,14) = -1.*x4;
   shape(2,15) = 0.;

   shape(3,0) = 0.;
   shape(3,1) = -1.*x2;
   shape(3,2) = -1.*x3;
   shape(3,3) = 1. - 1.*x1 - 1.*x4;
   shape(3,4) = x2;
   shape(3,5) = 0.;
   shape(3,6) = 0.;
   shape(3,7) = -1.*x2;
   shape(3,8) = x3;
   shape(3,9) = 0.;
   shape(3,10) = 0.;
   shape(3,11) = -1.*x3;
   shape(3,12) = -1. + x1 + x4;
   shape(3,13) = x2;
   shape(3,14) = x3;
   shape(3,15) = 0.;

   shape(4,0) = 0.;
   shape(4,1) = x2;
   shape(4,2) = -1. + x1 + x3;
   shape(4,3) = x4;
   shape(4,4) = -1.*x2;
   shape(4,5) = 0.;
   shape(4,6) = x2;
   shape(4,7) = 0.;
   shape(4,8) = 1. - 1.*x1 - 1.*x3;
   shape(4,9) = -1.*x2;
   shape(4,10) = 0.;
   shape(4,11) = -1.*x4;
   shape(4,12) = -1.*x4;
   shape(4,13) = 0.;
   shape(4,14) = x4;
   shape(4,15) = 0.;

   shape(5,0) = 0.;
   shape(5,1) = 1. - 1.*x1 - 1.*x2;
   shape(5,2) = -1.*x3;
   shape(5,3) = -1.*x4;
   shape(5,4) = -1. + x1 + x2;
   shape(5,5) = 0.;
   shape(5,6) = x3;
   shape(5,7) = x4;
   shape(5,8) = x3;
   shape(5,9) = -1.*x3;
   shape(5,10) = 0.;
   shape(5,11) = 0.;
   shape(5,12) = x4;
   shape(5,13) = -1.*x4;
   shape(5,14) = 0.;
   shape(5,15) = 0.;

   shape(6,0) = 0.;
   shape(6,1) = 0.;
   shape(6,2) = 0.;
   shape(6,3) = x1;
   shape(6,4) = 0.;
   shape(6,5) = 0.;
   shape(6,6) = 0.;
   shape(6,7) = x2;
   shape(6,8) = 0.;
   shape(6,9) = 0.;
   shape(6,10) = 0.;
   shape(6,11) = x3;
   shape(6,12) = -1.*x1;
   shape(6,13) = -1.*x2;
   shape(6,14) = -1.*x3;
   shape(6,15) = 0.;

   shape(7,0) = 0.;
   shape(7,1) = 0.;
   shape(7,2) = -1.*x1;
   shape(7,3) = 0.;
   shape(7,4) = 0.;
   shape(7,5) = 0.;
   shape(7,6) = -1.*x2;
   shape(7,7) = 0.;
   shape(7,8) = x1;
   shape(7,9) = x2;
   shape(7,10) = 0.;
   shape(7,11) = x4;
   shape(7,12) = 0.;
   shape(7,13) = 0.;
   shape(7,14) = -1.*x4;
   shape(7,15) = 0.;

   shape(8,0) = 0.;
   shape(8,1) = x1;
   shape(8,2) = 0.;
   shape(8,3) = 0.;
   shape(8,4) = -1.*x1;
   shape(8,5) = 0.;
   shape(8,6) = -1.*x3;
   shape(8,7) = -1.*x4;
   shape(8,8) = 0.;
   shape(8,9) = x3;
   shape(8,10) = 0.;
   shape(8,11) = 0.;
   shape(8,12) = 0.;
   shape(8,13) = x4;
   shape(8,14) = 0.;
   shape(8,15) = 0.;

   shape(9,0) = 0.;
   shape(9,1) = x2;
   shape(9,2) = x3;
   shape(9,3) = x4;
   shape(9,4) = -1.*x2;
   shape(9,5) = 0.;
   shape(9,6) = 0.;
   shape(9,7) = 0.;
   shape(9,8) = -1.*x3;
   shape(9,9) = 0.;
   shape(9,10) = 0.;
   shape(9,11) = 0.;
   shape(9,12) = -1.*x4;
   shape(9,13) = 0.;
   shape(9,14) = 0.;
   shape(9,15) = 0.;

}

void DivSkew1PentFiniteElement::CalcDivSkewShape(const IntegrationPoint &ip,
                                                 DenseMatrix &divSkew_shape) const
{
   divSkew_shape(0,0) = 0.;
   divSkew_shape(0,1) = 0.;
   divSkew_shape(0,2) = -3.;
   divSkew_shape(0,3) = 3.;

   divSkew_shape(1,0) = 0.;
   divSkew_shape(1,1) = 3.;
   divSkew_shape(1,2) = 0.;
   divSkew_shape(1,3) = -3.;

   divSkew_shape(2,0) = 0.;
   divSkew_shape(2,1) = -3.;
   divSkew_shape(2,2) = 3.;
   divSkew_shape(2,3) = 0.;

   divSkew_shape(3,0) = -3.;
   divSkew_shape(3,1) = 0.;
   divSkew_shape(3,2) = 0.;
   divSkew_shape(3,3) = 3.;

   divSkew_shape(4,0) = 3.;
   divSkew_shape(4,1) = 0.;
   divSkew_shape(4,2) = -3.;
   divSkew_shape(4,3) = 0.;

   divSkew_shape(5,0) = -3.;
   divSkew_shape(5,1) = 3.;
   divSkew_shape(5,2) = 0.;
   divSkew_shape(5,3) = 0.;

   divSkew_shape(6,0) = 0.;
   divSkew_shape(6,1) = 0.;
   divSkew_shape(6,2) = 0.;
   divSkew_shape(6,3) = -3.;

   divSkew_shape(7,0) = 0.;
   divSkew_shape(7,1) = 0.;
   divSkew_shape(7,2) = 3.;
   divSkew_shape(7,3) = 0.;

   divSkew_shape(8,0) = 0.;
   divSkew_shape(8,1) = -3.;
   divSkew_shape(8,2) = 0.;
   divSkew_shape(8,3) = 0.;

   divSkew_shape(9,0) = 3.;
   divSkew_shape(9,1) = 0.;
   divSkew_shape(9,2) = 0.;
   divSkew_shape(9,3) = 0.;
}

void DivSkew1PentFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   Vector v(6);
   double t1[4]; Vector t1i(t1, 4);
   double t2[4]; Vector t2i(t2, 4);
   Vector Mt(4);
   DenseMatrix mat(4,4); mat = 0.0;

   dofs.SetSize(10); dofs = 0.0;
   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      const DenseMatrix &J = Trans.Jacobian();

      vc.Eval(v, Trans, Nodes.IntPoint (k));

      mat(0,1) =  v(5); mat(0,2) = -v(4); mat(0,3) =  v(3);
      mat(1,0) = -v(5);                   mat(1,2) =  v(2); mat(1,3) = -v(1);
      mat(2,0) =  v(4); mat(2,1) = -v(2);                   mat(2,3) =  v(0);
      mat(3,0) = -v(3); mat(3,1) =  v(1); mat(3,2) = -v(0);

      J.Mult(tk1[k],t1);
      J.Mult(tk2[k],t2);

      mat.Mult(t2i, Mt);

      dofs(k) = t1i * Mt;
   }
}

void DivSkew1PentFiniteElement::Project(const FiniteElement &fe,
                                        ElementTransformation &Trans,
                                        DenseMatrix &Id) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   Id.SetSize(10,6*dof); Id = 0.0;

   double t1[4]; Vector t1i(t1, 4);
   double t2[4]; Vector t2i(t2, 4);
   Vector Mt(4);
   DenseMatrix mat(4,4); mat = 0.0;

   int indI[6] = {0,0,0,1,1,2};
   int indJ[6] = {1,2,3,2,3,3};

   Vector shape(dof);
   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      fe.CalcShape(Nodes.IntPoint(k), shape);

      J.Mult(tk1[k],t1);
      J.Mult(tk2[k],t2);

      for (int di=0; di<6; di++)
      {
         for (int j=0; j<dof; j++)
         {
            mat = 0.0;
            mat(indI[di], indJ[di]) =  shape(j);
            mat(indJ[di], indI[di]) = -shape(j);

            mat.Mult(t2i, Mt);

            Id(k, di*dof+j) = t1i * Mt;
         }
      }
   }
}


void DivSkew1PentFiniteElement::ProjectCurl(const FiniteElement &fe,
                                            ElementTransformation &Trans,
                                            DenseMatrix &curl) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   curl.SetSize(10,dof); curl = 0.0;

   double t1[4]; Vector t1i(t1, 4);
   double t2[4]; Vector t2i(t2, 4);
   Vector Mt(4);
   DenseMatrix mat(4,4); mat = 0.0;

   DenseMatrix tSh(4,4);
   DenseMatrix trShTemp(4,4);
   DenseMatrix invJ(4,4);
   DenseMatrix invJtr(invJ);
   DenseMatrix curlShape(dof, 6);

   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      CalcInverse(J, invJ);
      invJtr = invJ; invJtr.Transpose();

      J.Mult(tk1[k],t1);
      J.Mult(tk2[k],t2);

      //compute the curl
      fe.CalcCurlShape(Nodes.IntPoint(k), curlShape);
      for (int j=0; j<dof; j++)
      {
         tSh(0,1) =  curlShape(j,0); tSh(0,2) =  curlShape(j,1);
         tSh(0,3) =  curlShape(j,2);
         tSh(1,0) = -curlShape(j,0);
         tSh(1,2) =  curlShape(j,3); tSh(1,3) =  curlShape(j,4);
         tSh(2,0) = -curlShape(j,1); tSh(2,1) = -curlShape(j,3);
         tSh(2,3) =  curlShape(j,5);
         tSh(3,0) = -curlShape(j,2); tSh(3,1) = -curlShape(j,4);
         tSh(3,2) = -curlShape(j,5);

         Mult(tSh, invJ, trShTemp);
         Mult(invJtr, trShTemp, mat);

         mat.Mult(t2i, Mt);
         curl(k, j) = t1i * Mt;
      }
   }
}

void DivSkew1PentFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{

   //cerr << "DivSkew1PentFiniteElement::GetLocalInterpolation (...)  is not implemented yet!!! \n";

   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim*dim);
   DenseMatrix J(dim);
#else
   vshape.SetSize(dof, dim*dim);
   J.SetSize(dim);
#endif

   DenseMatrix mat(dim,dim); mat = 0.0;
   Vector t1i(4);
   Vector t2i(4);
   Vector Mt(4);

#ifdef MFEM_DEBUG
   for (k = 0; k < 10; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);

      for (j = 0; j < 10; j++)
      {
         mat(0,0) =  0.0;          mat(0,1) =  vshape(j,11);
         mat(0,2) =  vshape(j,13); mat(0,3) =  vshape(j,6);
         mat(1,0) =  vshape(j,14); mat(1,1) =  0.0;
         mat(1,2) =  vshape(j,3);  mat(1,3) =  vshape(j,8);
         mat(2,0) =  vshape(j,7);  mat(2,1) =  vshape(j,12); mat(2,2) =  0.0;
         mat(2,3) =  vshape(j,1);
         mat(3,0) =  vshape(j,9);  mat(3,1) =  vshape(j,2);
         mat(3,2) =  vshape(j,4);  mat(3,3) =  0.0;

         for ( int dim = 0; dim < 4; ++dim)
         {
             t1i[dim] = tk1[k][dim];
             t2i[dim] = tk2[k][dim];
         }
         mat.Mult(t2i, Mt);
         double d = t1i * Mt;

         //double d = ( vshape(j,0)*nk[k][0] + vshape(j,1)*nk[k][1] +
                      //vshape(j,2)*nk[k][2] + vshape(j,3)*nk[k][3] );
         if (j == k) { d -= 1.0; }

         if (fabs(d) > 1.0e-12)
         {
            cerr << "DivSkew1PentFiniteElement::GetLocalInterpolation (...)\n"
                 " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = ip.t = 0.0;
   Trans.SetIntPoint (&ip);
   J = Trans.Jacobian();
   Vector vk1(4), vk2(4);
   Vector xk(4);

   for (k = 0; k < 10; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = xk[0]; ip.y = xk[1]; ip.z = xk[2]; ip.t = xk[3];
      CalcVShape (ip, vshape);

      //  vk1 = J tk1, vk2 =  J tk2
      for ( int dim = 0; dim < 4; ++dim)
      {
          t1i[dim] = tk1[k][dim];
          t2i[dim] = tk2[k][dim];
      }

      J.Mult(t1i, vk1);
      J.Mult(t2i, vk2);

      for (j = 0; j < 10; j++)
      {
         mat(0,0) =  0.0;          mat(0,1) =  vshape(j,11);
         mat(0,2) =  vshape(j,13); mat(0,3) =  vshape(j,6);
         mat(1,0) =  vshape(j,14); mat(1,1) =  0.0;
         mat(1,2) =  vshape(j,3);  mat(1,3) =  vshape(j,8);
         mat(2,0) =  vshape(j,7);  mat(2,1) =  vshape(j,12); mat(2,2) =  0.0;
         mat(2,3) =  vshape(j,1);
         mat(3,0) =  vshape(j,9);  mat(3,1) =  vshape(j,2);
         mat(3,2) =  vshape(j,4);  mat(3,3) =  0.0;

         mat.Mult(vk2, Mt);
         I(k,j) = vk1 * Mt;

         if (fabs (I(k,j)) < 1.0e-12)
         {
             I(k,j) = 0.0;
         }
      }
   }
}


RT0HexFiniteElement::RT0HexFiniteElement()
   : VectorFiniteElement(3, Geometry::CUBE, 6, 1, H_DIV, FunctionSpace::Qk)
{
   // not real nodes ...
   // z = 0, y = 0, x = 1, y = 1, x = 0, z = 1
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.5;
   Nodes.IntPoint(0).z = 0.0;

   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.5;

   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 0.5;
   Nodes.IntPoint(2).z = 0.5;

   Nodes.IntPoint(3).x = 0.5;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(3).z = 0.5;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(4).z = 0.5;

   Nodes.IntPoint(5).x = 0.5;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 1.0;
}

void RT0HexFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;
   // z = 0
   shape(0,0) = 0.;
   shape(0,1) = 0.;
   shape(0,2) = z - 1.;
   // y = 0
   shape(1,0) = 0.;
   shape(1,1) = y - 1.;
   shape(1,2) = 0.;
   // x = 1
   shape(2,0) = x;
   shape(2,1) = 0.;
   shape(2,2) = 0.;
   // y = 1
   shape(3,0) = 0.;
   shape(3,1) = y;
   shape(3,2) = 0.;
   // x = 0
   shape(4,0) = x - 1.;
   shape(4,1) = 0.;
   shape(4,2) = 0.;
   // z = 1
   shape(5,0) = 0.;
   shape(5,1) = 0.;
   shape(5,2) = z;
}

void RT0HexFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
   divshape(0) = 1.;
   divshape(1) = 1.;
   divshape(2) = 1.;
   divshape(3) = 1.;
   divshape(4) = 1.;
   divshape(5) = 1.;
}

const double RT0HexFiniteElement::nk[6][3] =
{{0,0,-1}, {0,-1,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,0,1}};

void RT0HexFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 6; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 6; j++)
      {
         double d = ( vshape(j,0)*nk[k][0] + vshape(j,1)*nk[k][1] +
                      vshape(j,2)*nk[k][2] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT0HexFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[3];
   Vector xk (vk, 3);

   for (k = 0; k < 6; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2];
      vk[2] = Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2];
      for (j = 0; j < 6; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT0HexFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[3];
   Vector xk (vk, 3);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 6; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) =
         vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2] ) +
         vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2] ) +
         vk[2] * ( Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2] );
   }
}

RT1HexFiniteElement::RT1HexFiniteElement()
   : VectorFiniteElement(3, Geometry::CUBE, 36, 2, H_DIV, FunctionSpace::Qk)
{
   // z = 0
   Nodes.IntPoint(2).x  = 1./3.;
   Nodes.IntPoint(2).y  = 1./3.;
   Nodes.IntPoint(2).z  = 0.0;
   Nodes.IntPoint(3).x  = 2./3.;
   Nodes.IntPoint(3).y  = 1./3.;
   Nodes.IntPoint(3).z  = 0.0;
   Nodes.IntPoint(0).x  = 1./3.;
   Nodes.IntPoint(0).y  = 2./3.;
   Nodes.IntPoint(0).z  = 0.0;
   Nodes.IntPoint(1).x  = 2./3.;
   Nodes.IntPoint(1).y  = 2./3.;
   Nodes.IntPoint(1).z  = 0.0;
   // y = 0
   Nodes.IntPoint(4).x  = 1./3.;
   Nodes.IntPoint(4).y  = 0.0;
   Nodes.IntPoint(4).z  = 1./3.;
   Nodes.IntPoint(5).x  = 2./3.;
   Nodes.IntPoint(5).y  = 0.0;
   Nodes.IntPoint(5).z  = 1./3.;
   Nodes.IntPoint(6).x  = 1./3.;
   Nodes.IntPoint(6).y  = 0.0;
   Nodes.IntPoint(6).z  = 2./3.;
   Nodes.IntPoint(7).x  = 2./3.;
   Nodes.IntPoint(7).y  = 0.0;
   Nodes.IntPoint(7).z  = 2./3.;
   // x = 1
   Nodes.IntPoint(8).x  = 1.0;
   Nodes.IntPoint(8).y  = 1./3.;
   Nodes.IntPoint(8).z  = 1./3.;
   Nodes.IntPoint(9).x  = 1.0;
   Nodes.IntPoint(9).y  = 2./3.;
   Nodes.IntPoint(9).z  = 1./3.;
   Nodes.IntPoint(10).x = 1.0;
   Nodes.IntPoint(10).y = 1./3.;
   Nodes.IntPoint(10).z = 2./3.;
   Nodes.IntPoint(11).x = 1.0;
   Nodes.IntPoint(11).y = 2./3.;
   Nodes.IntPoint(11).z = 2./3.;
   // y = 1
   Nodes.IntPoint(13).x = 1./3.;
   Nodes.IntPoint(13).y = 1.0;
   Nodes.IntPoint(13).z = 1./3.;
   Nodes.IntPoint(12).x = 2./3.;
   Nodes.IntPoint(12).y = 1.0;
   Nodes.IntPoint(12).z = 1./3.;
   Nodes.IntPoint(15).x = 1./3.;
   Nodes.IntPoint(15).y = 1.0;
   Nodes.IntPoint(15).z = 2./3.;
   Nodes.IntPoint(14).x = 2./3.;
   Nodes.IntPoint(14).y = 1.0;
   Nodes.IntPoint(14).z = 2./3.;
   // x = 0
   Nodes.IntPoint(17).x = 0.0;
   Nodes.IntPoint(17).y = 1./3.;
   Nodes.IntPoint(17).z = 1./3.;
   Nodes.IntPoint(16).x = 0.0;
   Nodes.IntPoint(16).y = 2./3.;
   Nodes.IntPoint(16).z = 1./3.;
   Nodes.IntPoint(19).x = 0.0;
   Nodes.IntPoint(19).y = 1./3.;
   Nodes.IntPoint(19).z = 2./3.;
   Nodes.IntPoint(18).x = 0.0;
   Nodes.IntPoint(18).y = 2./3.;
   Nodes.IntPoint(18).z = 2./3.;
   // z = 1
   Nodes.IntPoint(20).x = 1./3.;
   Nodes.IntPoint(20).y = 1./3.;
   Nodes.IntPoint(20).z = 1.0;
   Nodes.IntPoint(21).x = 2./3.;
   Nodes.IntPoint(21).y = 1./3.;
   Nodes.IntPoint(21).z = 1.0;
   Nodes.IntPoint(22).x = 1./3.;
   Nodes.IntPoint(22).y = 2./3.;
   Nodes.IntPoint(22).z = 1.0;
   Nodes.IntPoint(23).x = 2./3.;
   Nodes.IntPoint(23).y = 2./3.;
   Nodes.IntPoint(23).z = 1.0;
   // x = 0.5 (interior)
   Nodes.IntPoint(24).x = 0.5;
   Nodes.IntPoint(24).y = 1./3.;
   Nodes.IntPoint(24).z = 1./3.;
   Nodes.IntPoint(25).x = 0.5;
   Nodes.IntPoint(25).y = 1./3.;
   Nodes.IntPoint(25).z = 2./3.;
   Nodes.IntPoint(26).x = 0.5;
   Nodes.IntPoint(26).y = 2./3.;
   Nodes.IntPoint(26).z = 1./3.;
   Nodes.IntPoint(27).x = 0.5;
   Nodes.IntPoint(27).y = 2./3.;
   Nodes.IntPoint(27).z = 2./3.;
   // y = 0.5 (interior)
   Nodes.IntPoint(28).x = 1./3.;
   Nodes.IntPoint(28).y = 0.5;
   Nodes.IntPoint(28).z = 1./3.;
   Nodes.IntPoint(29).x = 1./3.;
   Nodes.IntPoint(29).y = 0.5;
   Nodes.IntPoint(29).z = 2./3.;
   Nodes.IntPoint(30).x = 2./3.;
   Nodes.IntPoint(30).y = 0.5;
   Nodes.IntPoint(30).z = 1./3.;
   Nodes.IntPoint(31).x = 2./3.;
   Nodes.IntPoint(31).y = 0.5;
   Nodes.IntPoint(31).z = 2./3.;
   // z = 0.5 (interior)
   Nodes.IntPoint(32).x = 1./3.;
   Nodes.IntPoint(32).y = 1./3.;
   Nodes.IntPoint(32).z = 0.5;
   Nodes.IntPoint(33).x = 1./3.;
   Nodes.IntPoint(33).y = 2./3.;
   Nodes.IntPoint(33).z = 0.5;
   Nodes.IntPoint(34).x = 2./3.;
   Nodes.IntPoint(34).y = 1./3.;
   Nodes.IntPoint(34).z = 0.5;
   Nodes.IntPoint(35).x = 2./3.;
   Nodes.IntPoint(35).y = 2./3.;
   Nodes.IntPoint(35).z = 0.5;
}

void RT1HexFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   double x = ip.x, y = ip.y, z = ip.z;
   // z = 0
   shape(2,0)  = 0.;
   shape(2,1)  = 0.;
   shape(2,2)  = -(1. - 3.*z + 2.*z*z)*( 2. - 3.*x)*( 2. - 3.*y);
   shape(3,0)  = 0.;
   shape(3,1)  = 0.;
   shape(3,2)  = -(1. - 3.*z + 2.*z*z)*(-1. + 3.*x)*( 2. - 3.*y);
   shape(0,0)  = 0.;
   shape(0,1)  = 0.;
   shape(0,2)  = -(1. - 3.*z + 2.*z*z)*( 2. - 3.*x)*(-1. + 3.*y);
   shape(1,0)  = 0.;
   shape(1,1)  = 0.;
   shape(1,2)  = -(1. - 3.*z + 2.*z*z)*(-1. + 3.*x)*(-1. + 3.*y);
   // y = 0
   shape(4,0)  = 0.;
   shape(4,1)  = -(1. - 3.*y + 2.*y*y)*( 2. - 3.*x)*( 2. - 3.*z);
   shape(4,2)  = 0.;
   shape(5,0)  = 0.;
   shape(5,1)  = -(1. - 3.*y + 2.*y*y)*(-1. + 3.*x)*( 2. - 3.*z);
   shape(5,2)  = 0.;
   shape(6,0)  = 0.;
   shape(6,1)  = -(1. - 3.*y + 2.*y*y)*( 2. - 3.*x)*(-1. + 3.*z);
   shape(6,2)  = 0.;
   shape(7,0)  = 0.;
   shape(7,1)  = -(1. - 3.*y + 2.*y*y)*(-1. + 3.*x)*(-1. + 3.*z);
   shape(7,2)  = 0.;
   // x = 1
   shape(8,0)  = (-x + 2.*x*x)*( 2. - 3.*y)*( 2. - 3.*z);
   shape(8,1)  = 0.;
   shape(8,2)  = 0.;
   shape(9,0)  = (-x + 2.*x*x)*(-1. + 3.*y)*( 2. - 3.*z);
   shape(9,1)  = 0.;
   shape(9,2)  = 0.;
   shape(10,0) = (-x + 2.*x*x)*( 2. - 3.*y)*(-1. + 3.*z);
   shape(10,1) = 0.;
   shape(10,2) = 0.;
   shape(11,0) = (-x + 2.*x*x)*(-1. + 3.*y)*(-1. + 3.*z);
   shape(11,1) = 0.;
   shape(11,2) = 0.;
   // y = 1
   shape(13,0) = 0.;
   shape(13,1) = (-y + 2.*y*y)*( 2. - 3.*x)*( 2. - 3.*z);
   shape(13,2) = 0.;
   shape(12,0) = 0.;
   shape(12,1) = (-y + 2.*y*y)*(-1. + 3.*x)*( 2. - 3.*z);
   shape(12,2) = 0.;
   shape(15,0) = 0.;
   shape(15,1) = (-y + 2.*y*y)*( 2. - 3.*x)*(-1. + 3.*z);
   shape(15,2) = 0.;
   shape(14,0) = 0.;
   shape(14,1) = (-y + 2.*y*y)*(-1. + 3.*x)*(-1. + 3.*z);
   shape(14,2) = 0.;
   // x = 0
   shape(17,0) = -(1. - 3.*x + 2.*x*x)*( 2. - 3.*y)*( 2. - 3.*z);
   shape(17,1) = 0.;
   shape(17,2) = 0.;
   shape(16,0) = -(1. - 3.*x + 2.*x*x)*(-1. + 3.*y)*( 2. - 3.*z);
   shape(16,1) = 0.;
   shape(16,2) = 0.;
   shape(19,0) = -(1. - 3.*x + 2.*x*x)*( 2. - 3.*y)*(-1. + 3.*z);
   shape(19,1) = 0.;
   shape(19,2) = 0.;
   shape(18,0) = -(1. - 3.*x + 2.*x*x)*(-1. + 3.*y)*(-1. + 3.*z);
   shape(18,1) = 0.;
   shape(18,2) = 0.;
   // z = 1
   shape(20,0) = 0.;
   shape(20,1) = 0.;
   shape(20,2) = (-z + 2.*z*z)*( 2. - 3.*x)*( 2. - 3.*y);
   shape(21,0) = 0.;
   shape(21,1) = 0.;
   shape(21,2) = (-z + 2.*z*z)*(-1. + 3.*x)*( 2. - 3.*y);
   shape(22,0) = 0.;
   shape(22,1) = 0.;
   shape(22,2) = (-z + 2.*z*z)*( 2. - 3.*x)*(-1. + 3.*y);
   shape(23,0) = 0.;
   shape(23,1) = 0.;
   shape(23,2) = (-z + 2.*z*z)*(-1. + 3.*x)*(-1. + 3.*y);
   // x = 0.5 (interior)
   shape(24,0) = (4.*x - 4.*x*x)*( 2. - 3.*y)*( 2. - 3.*z);
   shape(24,1) = 0.;
   shape(24,2) = 0.;
   shape(25,0) = (4.*x - 4.*x*x)*( 2. - 3.*y)*(-1. + 3.*z);
   shape(25,1) = 0.;
   shape(25,2) = 0.;
   shape(26,0) = (4.*x - 4.*x*x)*(-1. + 3.*y)*( 2. - 3.*z);
   shape(26,1) = 0.;
   shape(26,2) = 0.;
   shape(27,0) = (4.*x - 4.*x*x)*(-1. + 3.*y)*(-1. + 3.*z);
   shape(27,1) = 0.;
   shape(27,2) = 0.;
   // y = 0.5 (interior)
   shape(28,0) = 0.;
   shape(28,1) = (4.*y - 4.*y*y)*( 2. - 3.*x)*( 2. - 3.*z);
   shape(28,2) = 0.;
   shape(29,0) = 0.;
   shape(29,1) = (4.*y - 4.*y*y)*( 2. - 3.*x)*(-1. + 3.*z);
   shape(29,2) = 0.;
   shape(30,0) = 0.;
   shape(30,1) = (4.*y - 4.*y*y)*(-1. + 3.*x)*( 2. - 3.*z);
   shape(30,2) = 0.;
   shape(31,0) = 0.;
   shape(31,1) = (4.*y - 4.*y*y)*(-1. + 3.*x)*(-1. + 3.*z);
   shape(31,2) = 0.;
   // z = 0.5 (interior)
   shape(32,0) = 0.;
   shape(32,1) = 0.;
   shape(32,2) = (4.*z - 4.*z*z)*( 2. - 3.*x)*( 2. - 3.*y);
   shape(33,0) = 0.;
   shape(33,1) = 0.;
   shape(33,2) = (4.*z - 4.*z*z)*( 2. - 3.*x)*(-1. + 3.*y);
   shape(34,0) = 0.;
   shape(34,1) = 0.;
   shape(34,2) = (4.*z - 4.*z*z)*(-1. + 3.*x)*( 2. - 3.*y);
   shape(35,0) = 0.;
   shape(35,1) = 0.;
   shape(35,2) = (4.*z - 4.*z*z)*(-1. + 3.*x)*(-1. + 3.*y);
}

void RT1HexFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
   double x = ip.x, y = ip.y, z = ip.z;
   // z = 0
   divshape(2)  = -(-3. + 4.*z)*( 2. - 3.*x)*( 2. - 3.*y);
   divshape(3)  = -(-3. + 4.*z)*(-1. + 3.*x)*( 2. - 3.*y);
   divshape(0)  = -(-3. + 4.*z)*( 2. - 3.*x)*(-1. + 3.*y);
   divshape(1)  = -(-3. + 4.*z)*(-1. + 3.*x)*(-1. + 3.*y);
   // y = 0
   divshape(4)  = -(-3. + 4.*y)*( 2. - 3.*x)*( 2. - 3.*z);
   divshape(5)  = -(-3. + 4.*y)*(-1. + 3.*x)*( 2. - 3.*z);
   divshape(6)  = -(-3. + 4.*y)*( 2. - 3.*x)*(-1. + 3.*z);
   divshape(7)  = -(-3. + 4.*y)*(-1. + 3.*x)*(-1. + 3.*z);
   // x = 1
   divshape(8)  = (-1. + 4.*x)*( 2. - 3.*y)*( 2. - 3.*z);
   divshape(9)  = (-1. + 4.*x)*(-1. + 3.*y)*( 2. - 3.*z);
   divshape(10) = (-1. + 4.*x)*( 2. - 3.*y)*(-1. + 3.*z);
   divshape(11) = (-1. + 4.*x)*(-1. + 3.*y)*(-1. + 3.*z);
   // y = 1
   divshape(13) = (-1. + 4.*y)*( 2. - 3.*x)*( 2. - 3.*z);
   divshape(12) = (-1. + 4.*y)*(-1. + 3.*x)*( 2. - 3.*z);
   divshape(15) = (-1. + 4.*y)*( 2. - 3.*x)*(-1. + 3.*z);
   divshape(14) = (-1. + 4.*y)*(-1. + 3.*x)*(-1. + 3.*z);
   // x = 0
   divshape(17) = -(-3. + 4.*x)*( 2. - 3.*y)*( 2. - 3.*z);
   divshape(16) = -(-3. + 4.*x)*(-1. + 3.*y)*( 2. - 3.*z);
   divshape(19) = -(-3. + 4.*x)*( 2. - 3.*y)*(-1. + 3.*z);
   divshape(18) = -(-3. + 4.*x)*(-1. + 3.*y)*(-1. + 3.*z);
   // z = 1
   divshape(20) = (-1. + 4.*z)*( 2. - 3.*x)*( 2. - 3.*y);
   divshape(21) = (-1. + 4.*z)*(-1. + 3.*x)*( 2. - 3.*y);
   divshape(22) = (-1. + 4.*z)*( 2. - 3.*x)*(-1. + 3.*y);
   divshape(23) = (-1. + 4.*z)*(-1. + 3.*x)*(-1. + 3.*y);
   // x = 0.5 (interior)
   divshape(24) = ( 4. - 8.*x)*( 2. - 3.*y)*( 2. - 3.*z);
   divshape(25) = ( 4. - 8.*x)*( 2. - 3.*y)*(-1. + 3.*z);
   divshape(26) = ( 4. - 8.*x)*(-1. + 3.*y)*( 2. - 3.*z);
   divshape(27) = ( 4. - 8.*x)*(-1. + 3.*y)*(-1. + 3.*z);
   // y = 0.5 (interior)
   divshape(28) = ( 4. - 8.*y)*( 2. - 3.*x)*( 2. - 3.*z);
   divshape(29) = ( 4. - 8.*y)*( 2. - 3.*x)*(-1. + 3.*z);
   divshape(30) = ( 4. - 8.*y)*(-1. + 3.*x)*( 2. - 3.*z);
   divshape(31) = ( 4. - 8.*y)*(-1. + 3.*x)*(-1. + 3.*z);
   // z = 0.5 (interior)
   divshape(32) = ( 4. - 8.*z)*( 2. - 3.*x)*( 2. - 3.*y);
   divshape(33) = ( 4. - 8.*z)*( 2. - 3.*x)*(-1. + 3.*y);
   divshape(34) = ( 4. - 8.*z)*(-1. + 3.*x)*( 2. - 3.*y);
   divshape(35) = ( 4. - 8.*z)*(-1. + 3.*x)*(-1. + 3.*y);
}

const double RT1HexFiniteElement::nk[36][3] =
{
   {0, 0,-1}, {0, 0,-1}, {0, 0,-1}, {0, 0,-1},
   {0,-1, 0}, {0,-1, 0}, {0,-1, 0}, {0,-1, 0},
   {1, 0, 0}, {1, 0, 0}, {1, 0, 0}, {1, 0, 0},
   {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0},
   {-1,0, 0}, {-1,0, 0}, {-1,0, 0}, {-1,0, 0},
   {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1},
   {1, 0, 0}, {1, 0, 0}, {1, 0, 0}, {1, 0, 0},
   {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0},
   {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}
};

void RT1HexFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 36; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 36; j++)
      {
         double d = ( vshape(j,0)*nk[k][0] + vshape(j,1)*nk[k][1] +
                      vshape(j,2)*nk[k][2] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT0HexFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[3];
   Vector xk (vk, 3);

   for (k = 0; k < 36; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2];
      vk[2] = Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2];
      for (j = 0; j < 36; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT1HexFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[3];
   Vector xk (vk, 3);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 36; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) =
         vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2] ) +
         vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2] ) +
         vk[2] * ( Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2] );
   }
}

RT0TetFiniteElement::RT0TetFiniteElement()
   : VectorFiniteElement(3, Geometry::TETRAHEDRON, 4, 1, H_DIV)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.33333333333333333333;
   Nodes.IntPoint(0).y = 0.33333333333333333333;
   Nodes.IntPoint(0).z = 0.33333333333333333333;

   Nodes.IntPoint(1).x = 0.0;
   Nodes.IntPoint(1).y = 0.33333333333333333333;
   Nodes.IntPoint(1).z = 0.33333333333333333333;

   Nodes.IntPoint(2).x = 0.33333333333333333333;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.33333333333333333333;

   Nodes.IntPoint(3).x = 0.33333333333333333333;
   Nodes.IntPoint(3).y = 0.33333333333333333333;
   Nodes.IntPoint(3).z = 0.0;
}

void RT0TetFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   double x2 = 2.0*ip.x, y2 = 2.0*ip.y, z2 = 2.0*ip.z;

   shape(0,0) = x2;
   shape(0,1) = y2;
   shape(0,2) = z2;

   shape(1,0) = x2 - 2.0;
   shape(1,1) = y2;
   shape(1,2) = z2;

   shape(2,0) = x2;
   shape(2,1) = y2 - 2.0;
   shape(2,2) = z2;

   shape(3,0) = x2;
   shape(3,1) = y2;
   shape(3,2) = z2 - 2.0;
}

void RT0TetFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
   divshape(0) = 6.0;
   divshape(1) = 6.0;
   divshape(2) = 6.0;
   divshape(3) = 6.0;
}

const double RT0TetFiniteElement::nk[4][3] =
{{.5,.5,.5}, {-.5,0,0}, {0,-.5,0}, {0,0,-.5}};

void RT0TetFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{
   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(Dof, Dim);
   DenseMatrix Jinv(Dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 4; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 4; j++)
      {
         double d = ( vshape(j,0)*nk[k][0] + vshape(j,1)*nk[k][1] +
                      vshape(j,2)*nk[k][2] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            mfem::err << "RT0TetFiniteElement::GetLocalInterpolation (...)\n"
                      " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[3];
   Vector xk (vk, 3);

   for (k = 0; k < 4; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2];
      vk[2] = Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2];
      for (j = 0; j < 4; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+
                             vshape(j,2)*vk[2])) < 1.0e-12)
         {
            I(k,j) = 0.0;
         }
   }
}

void RT0TetFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[3];
   Vector xk (vk, 3);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(Dim);
#endif

   for (int k = 0; k < 4; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose (Trans.Jacobian(), Jinv);

      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) =
         vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2] ) +
         vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2] ) +
         vk[2] * ( Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2] );
   }
}


RT0PentFiniteElement::RT0PentFiniteElement()
   : VectorFiniteElement(4, Geometry::PENTATOPE, 5, 1, H_DIV)
{
   // not real nodes ...
   Nodes.IntPoint(0).x = 0.25;
   Nodes.IntPoint(0).y = 0.25;
   Nodes.IntPoint(0).z = 0.25;
   Nodes.IntPoint(0).t = 0.0;

   Nodes.IntPoint(1).x = 0.25;
   Nodes.IntPoint(1).y = 0.25;
   Nodes.IntPoint(1).z = 0.0;
   Nodes.IntPoint(1).t = 0.25;

   Nodes.IntPoint(2).x = 0.25;
   Nodes.IntPoint(2).y = 0.0;
   Nodes.IntPoint(2).z = 0.25;
   Nodes.IntPoint(2).t = 0.25;

   Nodes.IntPoint(3).x = 0.0;
   Nodes.IntPoint(3).y = 0.25;
   Nodes.IntPoint(3).z = 0.25;
   Nodes.IntPoint(3).t = 0.25;

   Nodes.IntPoint(4).x = 0.25;
   Nodes.IntPoint(4).y = 0.25;
   Nodes.IntPoint(4).z = 0.25;
   Nodes.IntPoint(4).t = 0.25;
}

void RT0PentFiniteElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   double x3 = 3.0*ip.x, y3 = 3.0*ip.y, z3 = 3.0*ip.z, t3 = 3.0*ip.t;

   shape(0,0) = x3;
   shape(0,1) = y3;
   shape(0,2) = z3;
   shape(0,3) = t3 - 3.0;

   shape(1,0) = x3;
   shape(1,1) = y3;
   shape(1,2) = z3 - 3.0;
   shape(1,3) = t3;

   shape(2,0) = x3;
   shape(2,1) = y3 - 3.0;
   shape(2,2) = z3;
   shape(2,3) = t3;

   shape(3,0) = x3 - 3.0;
   shape(3,1) = y3;
   shape(3,2) = z3;
   shape(3,3) = t3;

   shape(4,0) = x3;
   shape(4,1) = y3;
   shape(4,2) = z3;
   shape(4,3) = t3;
    
    std::ofstream shape_file("shape_matrix.txt");
    shape.PrintMatlab(shape_file);
}

void RT0PentFiniteElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   divshape(0) = 12.0;
   divshape(1) = 12.0;
   divshape(2) = 12.0;
   divshape(3) = 12.0;
   divshape(4) = 12.0;
}

const double RT0PentFiniteElement::nk[5][4] =
{{0.0,0.0,0.0,-1./3.}, {0.0,0.0,-1./3.,0.0}, {0.0,-1./3.,0.0,0.0}, {-1./3.,0.0,0.0,0.0}, {1./3.,1./3.,1./3.,1./3.}};

void RT0PentFiniteElement::GetLocalInterpolation (
   ElementTransformation &Trans, DenseMatrix &I) const
{

   cerr << "RT0PentFiniteElement::GetLocalInterpolation (...) implementation not tested yet!!! \n";

   int k, j;
#ifdef MFEM_THREAD_SAFE
   DenseMatrix vshape(dof, dim);
   DenseMatrix Jinv(dim);
#endif

#ifdef MFEM_DEBUG
   for (k = 0; k < 5; k++)
   {
      CalcVShape (Nodes.IntPoint(k), vshape);
      for (j = 0; j < 5; j++)
      {
         double d = ( vshape(j,0)*nk[k][0] + vshape(j,1)*nk[k][1] +
                      vshape(j,2)*nk[k][2] + vshape(j,3)*nk[k][3] );
         if (j == k) { d -= 1.0; }
         if (fabs(d) > 1.0e-12)
         {
            cerr << "RT0PentFiniteElement::GetLocalInterpolation (...)\n"
                 " k = " << k << ", j = " << j << ", d = " << d << endl;
            mfem_error();
         }
      }
   }
#endif

   IntegrationPoint ip;
   ip.x = ip.y = ip.z = 0.0;
   Trans.SetIntPoint (&ip);
   // Trans must be linear
   // set Jinv = |J| J^{-t} = adj(J)^t
   CalcAdjugateTranspose (Trans.Jacobian(), Jinv);
   double vk[4];
   Vector xk(vk, 4);

   for (k = 0; k < 5; k++)
   {
      Trans.Transform (Nodes.IntPoint (k), xk);
      ip.x = vk[0]; ip.y = vk[1]; ip.z = vk[2]; ip.t = vk[3];
      CalcVShape (ip, vshape);
      //  vk = |J| J^{-t} nk
      vk[0] = Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2]+Jinv(0,
                                                                            3)*nk[k][3];
      vk[1] = Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2]+Jinv(1,
                                                                            3)*nk[k][3];
      vk[2] = Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2]+Jinv(2,
                                                                            3)*nk[k][3];
      vk[3] = Jinv(3,0)*nk[k][0]+Jinv(3,1)*nk[k][1]+Jinv(3,2)*nk[k][2]+Jinv(3,
                                                                            3)*nk[k][3];

      for (j = 0; j < 5; j++)
         if (fabs (I(k,j) = (vshape(j,0)*vk[0]+vshape(j,1)*vk[1]+vshape(j,
                                                                        2)*vk[2]+vshape(j,3)*vk[3])) < 1.0e-12) { I(k,j) = 0.0; }
   }
}

void RT0PentFiniteElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   double vk[4];
   Vector xk (vk, 4);
#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim);
#endif

   for (int k = 0; k < 5; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      // set Jinv = |J| J^{-t} = adj(J)^t
      CalcAdjugateTranspose(Trans.Jacobian(), Jinv);
      double w = Trans.Weight();


      vc.Eval (xk, Trans, Nodes.IntPoint (k));
      //  xk^t |J| J^{-t} nk
      dofs(k) =
         vk[0] * ( Jinv(0,0)*nk[k][0]+Jinv(0,1)*nk[k][1]+Jinv(0,2)*nk[k][2]+Jinv(0,
                                                                                 3)*nk[k][3] ) +
         vk[1] * ( Jinv(1,0)*nk[k][0]+Jinv(1,1)*nk[k][1]+Jinv(1,2)*nk[k][2]+Jinv(1,
                                                                                 3)*nk[k][3] ) +
         vk[2] * ( Jinv(2,0)*nk[k][0]+Jinv(2,1)*nk[k][1]+Jinv(2,2)*nk[k][2]+Jinv(2,
                                                                                 3)*nk[k][3] ) +
         vk[3] * ( Jinv(3,0)*nk[k][0]+Jinv(3,1)*nk[k][1]+Jinv(3,2)*nk[k][2]+Jinv(3,
                                                                                 3)*nk[k][3] );

      if (w < 0) { dofs(k) *= -1.0; }
   }
}

void RT0PentFiniteElement::Project(const FiniteElement &fe,
                                   ElementTransformation &Trans,
                                   DenseMatrix &Id) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   Id.SetSize(5,dim*dof); Id = 0.0;

   double n[4]; Vector ni(n, 4);
   Vector vecF(4);

   Vector shape(dof);
   for (int k = 0; k < 5; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      CalcAdjugateTranspose(J, Jinv);
      fe.CalcShape(Nodes.IntPoint(k), shape);

      Jinv.Mult(nk[k],n);

      for (int di=0; di<dim; di++)
      {
         for (int j=0; j<dof; j++)
         {
            vecF = 0.0; vecF(di) = shape(j);

            Id(k, di*dof+j) = vecF * ni;
         }
      }
   }
}

void RT0PentFiniteElement::ProjectDivSkew(const FiniteElement &fe,
                                          ElementTransformation &Trans,
                                          DenseMatrix &DivSkew) const
{
   int dim = fe.GetDim();
   int dof = fe.GetDof();

   DivSkew.SetSize(5,dof); DivSkew = 0.0;

   double n[4]; Vector ni(n, 4);
   Vector vecF(4);

   DenseMatrix DivSkewshape(dof,4);
   DenseMatrix DivSkew_dFt(dof,4);
   for (int k = 0; k < 5; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      CalcAdjugateTranspose(J, Jinv);

      fe.CalcDivSkewShape(Nodes.IntPoint(k), DivSkewshape);
      MultABt(DivSkewshape, J, DivSkew_dFt);
      DivSkew_dFt *= (1.0 / Trans.Weight());

      Jinv.Mult(nk[k],n);

      for (int j=0; j<dof; j++)
      {
         vecF(0) = DivSkew_dFt(j,0);
         vecF(1) = DivSkew_dFt(j,1);
         vecF(2) = DivSkew_dFt(j,2);
         vecF(3) = DivSkew_dFt(j,3);

         DivSkew(k, j) = vecF * ni;
      }
   }
}





RotTriLinearHexFiniteElement::RotTriLinearHexFiniteElement()
   : NodalFiniteElement(3, Geometry::CUBE, 6, 2, FunctionSpace::Qk)
{
   Nodes.IntPoint(0).x = 0.5;
   Nodes.IntPoint(0).y = 0.5;
   Nodes.IntPoint(0).z = 0.0;

   Nodes.IntPoint(1).x = 0.5;
   Nodes.IntPoint(1).y = 0.0;
   Nodes.IntPoint(1).z = 0.5;

   Nodes.IntPoint(2).x = 1.0;
   Nodes.IntPoint(2).y = 0.5;
   Nodes.IntPoint(2).z = 0.5;

   Nodes.IntPoint(3).x = 0.5;
   Nodes.IntPoint(3).y = 1.0;
   Nodes.IntPoint(3).z = 0.5;

   Nodes.IntPoint(4).x = 0.0;
   Nodes.IntPoint(4).y = 0.5;
   Nodes.IntPoint(4).z = 0.5;

   Nodes.IntPoint(5).x = 0.5;
   Nodes.IntPoint(5).y = 0.5;
   Nodes.IntPoint(5).z = 1.0;
}

void RotTriLinearHexFiniteElement::CalcShape(const IntegrationPoint &ip,
                                             Vector &shape) const
{
   double x = 2. * ip.x - 1.;
   double y = 2. * ip.y - 1.;
   double z = 2. * ip.z - 1.;
   double f5 = x * x - y * y;
   double f6 = y * y - z * z;

   shape(0) = (1./6.) * (1. - 3. * z -      f5 - 2. * f6);
   shape(1) = (1./6.) * (1. - 3. * y -      f5 +      f6);
   shape(2) = (1./6.) * (1. + 3. * x + 2. * f5 +      f6);
   shape(3) = (1./6.) * (1. + 3. * y -      f5 +      f6);
   shape(4) = (1./6.) * (1. - 3. * x + 2. * f5 +      f6);
   shape(5) = (1./6.) * (1. + 3. * z -      f5 - 2. * f6);
}

void RotTriLinearHexFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                              DenseMatrix &dshape) const
{
   const double a = 2./3.;

   double xt = a * (1. - 2. * ip.x);
   double yt = a * (1. - 2. * ip.y);
   double zt = a * (1. - 2. * ip.z);

   dshape(0,0) = xt;
   dshape(0,1) = yt;
   dshape(0,2) = -1. - 2. * zt;

   dshape(1,0) = xt;
   dshape(1,1) = -1. - 2. * yt;
   dshape(1,2) = zt;

   dshape(2,0) = 1. - 2. * xt;
   dshape(2,1) = yt;
   dshape(2,2) = zt;

   dshape(3,0) = xt;
   dshape(3,1) = 1. - 2. * yt;
   dshape(3,2) = zt;

   dshape(4,0) = -1. - 2. * xt;
   dshape(4,1) = yt;
   dshape(4,2) = zt;

   dshape(5,0) = xt;
   dshape(5,1) = yt;
   dshape(5,2) = 1. - 2. * zt;
}


Poly_1D::Basis::Basis(const int p, const double *nodes, EvalType etype)
   : etype(etype)
{
   switch (etype)
   {
      case ChangeOfBasis:
      {
         x.SetSize(p + 1);
         w.SetSize(p + 1);
         DenseMatrix A(p + 1);
         for (int i = 0; i <= p; i++)
         {
            CalcBasis(p, nodes[i], A.GetColumn(i));
         }
         Ai.Factor(A);
         // mfem::out << "Poly_1D::Basis(" << p << ",...) : "; Ai.TestInversion();
         break;
      }
      case Barycentric:
      {
         x.SetSize(p + 1);
         w.SetSize(p + 1);
         x = nodes;
         w = 1.0;
         for (int i = 0; i <= p; i++)
         {
            for (int j = 0; j < i; j++)
            {
               double xij = x(i) - x(j);
               w(i) *=  xij;
               w(j) *= -xij;
            }
         }
         for (int i = 0; i <= p; i++)
         {
            w(i) = 1.0/w(i);
         }

#ifdef MFEM_DEBUG
         // Make sure the nodes are increasing
         for (int i = 0; i < p; i++)
         {
            if (x(i) >= x(i+1))
            {
               mfem_error("Poly_1D::Basis::Basis : nodes are not increasing!");
            }
         }
#endif
         break;
      }
      case Positive:
         x.SetDataAndSize(NULL, p + 1); // use x to store (p + 1)
         break;

      default: break;
   }
}

void Poly_1D::Basis::Eval(const double y, Vector &u) const
{
   switch (etype)
   {
      case ChangeOfBasis:
      {
         CalcBasis(Ai.Width() - 1, y, x);
         Ai.Mult(x, u);
         break;
      }
      case Barycentric:
      {
         int i, k, p = x.Size() - 1;
         double l, lk;

         if (p == 0)
         {
            u(0) = 1.0;
            return;
         }

         lk = 1.0;
         for (k = 0; k < p; k++)
         {
            if (y >= (x(k) + x(k+1))/2)
            {
               lk *= y - x(k);
            }
            else
            {
               for (i = k+1; i <= p; i++)
               {
                  lk *= y - x(i);
               }
               break;
            }
         }
         l = lk * (y - x(k));

         for (i = 0; i < k; i++)
         {
            u(i) = l * w(i) / (y - x(i));
         }
         u(k) = lk * w(k);
         for (i++; i <= p; i++)
         {
            u(i) = l * w(i) / (y - x(i));
         }
         break;
      }
      case Positive:
         CalcBernstein(x.Size() - 1, y, u);
         break;

      default: break;
   }
}

void Poly_1D::Basis::Eval(const double y, Vector &u, Vector &d) const
{
   switch (etype)
   {
      case ChangeOfBasis:
      {
         CalcBasis(Ai.Width() - 1, y, x, w);
         Ai.Mult(x, u);
         Ai.Mult(w, d);
         break;
      }
      case Barycentric:
      {
         int i, k, p = x.Size() - 1;
         double l, lp, lk, sk, si;

         if (p == 0)
         {
            u(0) = 1.0;
            d(0) = 0.0;
            return;
         }

         lk = 1.0;
         for (k = 0; k < p; k++)
         {
            if (y >= (x(k) + x(k+1))/2)
            {
               lk *= y - x(k);
            }
            else
            {
               for (i = k+1; i <= p; i++)
               {
                  lk *= y - x(i);
               }
               break;
            }
         }
         l = lk * (y - x(k));

         sk = 0.0;
         for (i = 0; i < k; i++)
         {
            si = 1.0/(y - x(i));
            sk += si;
            u(i) = l * si * w(i);
         }
         u(k) = lk * w(k);
         for (i++; i <= p; i++)
         {
            si = 1.0/(y - x(i));
            sk += si;
            u(i) = l * si * w(i);
         }
         lp = l * sk + lk;

         for (i = 0; i < k; i++)
         {
            d(i) = (lp * w(i) - u(i))/(y - x(i));
         }
         d(k) = sk * u(k);
         for (i++; i <= p; i++)
         {
            d(i) = (lp * w(i) - u(i))/(y - x(i));
         }
         break;
      }
      case Positive:
         CalcBernstein(x.Size() - 1, y, u, d);
         break;

      default: break;
   }
}

void Poly_1D::Basis::Eval(const double y, Vector &u, Vector &d,
                          Vector &d2) const
{
   MFEM_VERIFY(etype == Barycentric,
               "Basis::Eval with second order derivatives not implemented for"
               " etype = " << etype);
   switch (etype)
   {
      case ChangeOfBasis:
      {
         CalcBasis(Ai.Width() - 1, y, x, w);
         Ai.Mult(x, u);
         Ai.Mult(w, d);
         // set d2 (not implemented yet)
         break;
      }
      case Barycentric:
      {
         int i, k, p = x.Size() - 1;
         double l, lp, lp2, lk, sk, si, sk2;

         if (p == 0)
         {
            u(0) = 1.0;
            d(0) = 0.0;
            d2(0) = 0.0;
            return;
         }

         lk = 1.0;
         for (k = 0; k < p; k++)
         {
            if (y >= (x(k) + x(k+1))/2)
            {
               lk *= y - x(k);
            }
            else
            {
               for (i = k+1; i <= p; i++)
               {
                  lk *= y - x(i);
               }
               break;
            }
         }
         l = lk * (y - x(k));

         sk = 0.0;
         sk2 = 0.0;
         for (i = 0; i < k; i++)
         {
            si = 1.0/(y - x(i));
            sk += si;
            sk2 -= si * si;
            u(i) = l * si * w(i);
         }
         u(k) = lk * w(k);
         for (i++; i <= p; i++)
         {
            si = 1.0/(y - x(i));
            sk += si;
            sk2 -= si * si;
            u(i) = l * si * w(i);
         }
         lp = l * sk + lk;
         lp2 = lp * sk + l * sk2 + sk * lk;

         for (i = 0; i < k; i++)
         {
            d(i) = (lp * w(i) - u(i))/(y - x(i));
            d2(i) = (lp2 * w(i) - 2 * d(i))/(y - x(i));
         }
         d(k) = sk * u(k);
         d2(k) = sk2 * u(k) + sk * d(k);
         for (i++; i <= p; i++)
         {
            d(i) = (lp * w(i) - u(i))/(y - x(i));
            d2(i) = (lp2 * w(i) - 2 * d(i))/(y - x(i));
         }
         break;
      }
      case Positive:
         CalcBernstein(x.Size() - 1, y, u, d);
         break;

      default: break;
   }
}

const int *Poly_1D::Binom(const int p)
{
   if (binom.NumCols() <= p)
   {
      binom.SetSize(p + 1, p + 1);
      for (int i = 0; i <= p; i++)
      {
         binom(i,0) = binom(i,i) = 1;
         for (int j = 1; j < i; j++)
         {
            binom(i,j) = binom(i-1,j) + binom(i-1,j-1);
         }
      }
   }
   return binom[p];
}

void Poly_1D::ChebyshevPoints(const int p, double *x)
{
   for (int i = 0; i <= p; i++)
   {
      // x[i] = 0.5*(1. + cos(M_PI*(p - i + 0.5)/(p + 1)));
      double s = sin(M_PI_2*(i + 0.5)/(p + 1));
      x[i] = s*s;
   }
}

void Poly_1D::CalcMono(const int p, const double x, double *u)
{
   double xn;
   u[0] = xn = 1.;
   for (int n = 1; n <= p; n++)
   {
      u[n] = (xn *= x);
   }
}

void Poly_1D::CalcMono(const int p, const double x, double *u, double *d)
{
   double xn;
   u[0] = xn = 1.;
   d[0] = 0.;
   for (int n = 1; n <= p; n++)
   {
      d[n] = n * xn;
      u[n] = (xn *= x);
   }
}

void Poly_1D::CalcBinomTerms(const int p, const double x, const double y,
                             double *u)
{
   if (p == 0)
   {
      u[0] = 1.;
   }
   else
   {
      int i;
      const int *b = Binom(p);
      double z = x;

      for (i = 1; i < p; i++)
      {
         u[i] = b[i]*z;
         z *= x;
      }
      u[p] = z;
      z = y;
      for (i--; i > 0; i--)
      {
         u[i] *= z;
         z *= y;
      }
      u[0] = z;
   }
}

void Poly_1D::CalcBinomTerms(const int p, const double x, const double y,
                             double *u, double *d)
{
   if (p == 0)
   {
      u[0] = 1.;
      d[0] = 0.;
   }
   else
   {
      int i;
      const int *b = Binom(p);
      const double xpy = x + y, ptx = p*x;
      double z = 1.;

      for (i = 1; i < p; i++)
      {
         d[i] = b[i]*z*(i*xpy - ptx);
         z *= x;
         u[i] = b[i]*z;
      }
      d[p] = p*z;
      u[p] = z*x;
      z = 1.;
      for (i--; i > 0; i--)
      {
         d[i] *= z;
         z *= y;
         u[i] *= z;
      }
      d[0] = -p*z;
      u[0] = z*y;
   }
}

void Poly_1D::CalcDBinomTerms(const int p, const double x, const double y,
                              double *d)
{
   if (p == 0)
   {
      d[0] = 0.;
   }
   else
   {
      int i;
      const int *b = Binom(p);
      const double xpy = x + y, ptx = p*x;
      double z = 1.;

      for (i = 1; i < p; i++)
      {
         d[i] = b[i]*z*(i*xpy - ptx);
         z *= x;
      }
      d[p] = p*z;
      z = 1.;
      for (i--; i > 0; i--)
      {
         d[i] *= z;
         z *= y;
      }
      d[0] = -p*z;
   }
}

void Poly_1D::CalcLegendre(const int p, const double x, double *u)
{
   // use the recursive definition for [-1,1]:
   // (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
   double z;
   u[0] = 1.;
   if (p == 0) { return; }
   u[1] = z = 2.*x - 1.;
   for (int n = 1; n < p; n++)
   {
      u[n+1] = ((2*n + 1)*z*u[n] - n*u[n-1])/(n + 1);
   }
}

void Poly_1D::CalcLegendre(const int p, const double x, double *u, double *d)
{
   // use the recursive definition for [-1,1]:
   // (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
   // for the derivative use, z in [-1,1]:
   // P'_{n+1}(z) = (2*n+1)*P_n(z)+P'_{n-1}(z)
   double z;
   u[0] = 1.;
   d[0] = 0.;
   if (p == 0) { return; }
   u[1] = z = 2.*x - 1.;
   d[1] = 2.;
   for (int n = 1; n < p; n++)
   {
      u[n+1] = ((2*n + 1)*z*u[n] - n*u[n-1])/(n + 1);
      d[n+1] = (4*n + 2)*u[n] + d[n-1];
   }
}

void Poly_1D::CalcLegendre(const int &p, const double &x, std::vector<double> &u, std::vector<double> &d)
{
   // use the recursive definition for [-1,1]:
   // (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
   // for the derivative use, z in [-1,1]:
   // P'_{n+1}(z) = (2*n+1)*P_n(z)+P'_{n-1}(z)
   double z;
   u.push_back(1.0);
   d.push_back(0.0);
   if (p == 0) { return; }
   z = 2.*x - 1.;
   u.push_back(z);
   d.push_back(2.0);
   for (int n = 1; n < p; n++)
   {
      u.push_back(((2*n + 1)*z*u[n] - n*u[n-1])/(n + 1));
      d.push_back((4*n + 2)*u[n] + d[n-1]);
   }
}

void Poly_1D::CalcLegendreShifted(const int &p, const double &x, const double &y, std::vector<double> &u)
{
    // use recursive definiton for [0,1] shifted Legendre Polynomials
    // (n+1)*P_{n+1} = (2x-1)(2*n+1)*z*P_n(z)-n*P_{n-1}(z)
    // ** Note the scaled and Shifted polynomials are otained when y != 1
    double z;
    u.push_back(1.0);
    if (p == 0) { return; }
    u.push_back(2.*x - y);
    z = 2.*x - y;
    for (int n = 1; n < p; n++)
    {
        u.push_back(((2*n + 1)*z*u[n] - n*(y*y)*u[n-1])/(n + 1));
    }
}

void Poly_1D::CalcScaledLegendreDerivative(const int &p, const double &x, const double &t,
                                        std::vector<double> &u,
                                        std::vector<double> &dudx, std::vector<double> &dudt)
{
   MFEM_ASSERT(p >= 0, "Polynomial order must be zero or larger");
   if (t > 0.0)
   {
      Poly_1D::CalcLegendre(p, x / t, u, dudx);
      dudx[0] = 0.0;
      dudt.push_back( - dudx[0] * x / t);
      for (int i = 1; i <= p; i++)
      {
         u[i]    *= pow(t, i);
         dudx[i] *= pow(t, i - 1);
         dudt.push_back( (u[i] * i - dudx[i] * x) / t);
      }
   }
   else
   {
      // This assumes x = 0 as well as t = 0 since x \in [0,t]
      u[0]    = 1.0;
      dudx[0] = 0.0;
      dudt.push_back(0.0);
      if (p >=1)
      {
         u[1]    =  0.0;
         dudx[1] =  2.0;
         dudt.push_back(-1.0);
      }
      for (int i = 2; i <= p; i++)
      {
         u[i] = 0.0;
         dudx[i] = 0.0;
         dudt.push_back(0.0);
      }
   }
}

void Poly_1D::CalcIntLegendre(const int &p, const double &x, const double &y, std::vector<double> &u)
{
    // use recursive definiton for [0,1] shifted, scaled and integrated Legendre Polynomials are
    // 2(2n+1)*L_{n+1} = P_{n+1}-P_{n-1}
    
    u.push_back(x);
    if (p == 1) { return; }
    for (int n = 1; n < p; n++)
    {
        double i = n+1;
        std::vector<double> P_i;
        std::vector<double> P_i_2;
        CalcLegendreShifted(i, x, y, P_i);
        CalcLegendreShifted(i-2, x, y, P_i_2);
        u.push_back((P_i[P_i.size()-1] - (y*y)*P_i_2[P_i_2.size()-1])/(2.*(2.*i-1)));
    }
}

void Poly_1D::CalcRLegendre(const int &p, const double &x, const double &y, std::vector<double> &u)
{
    //
    
    u.push_back(0);
    if (p == 0) { return; }
    for (int n = 1; n <= p; n++)
    {
        double i = n;
        std::vector<double> P_i;
        std::vector<double> P_i_1;
        CalcLegendreShifted(i, x, y, P_i);
        CalcLegendreShifted(i-1, x, y, P_i_1);
        u.push_back(-0.5*(P_i[P_i.size()-1] + y*P_i_1[P_i_1.size()-1]));
    }
}


void Poly_1D::CalcJacobi(const int &p, const double &x, const double &y, const double &alpha, std::vector<double> &u)
{
    
    //use the recursive definition for [0,1]:
    // a*P_i = b(c*(2*x-y)+alpha^2*y)P_(i-1)-d*y^2*P_(i-2)
    // ** Note the scaled and Shifted polynomials are otained when y != 1
    
    //u[0] = 1.;
    u.push_back(1.0);
    if (p == 0) { return; }
    //u[1] = 2.*x +alpha*x - 1.;
    u.push_back(2.*x -y +alpha*x);
    for (int n = 1; n < p; n++)
    {
        double i = n+1;
        double a = 2*i*(i+alpha)*(2*i+alpha-2);
        double b = 2*i+alpha-1;
        double c = (2*i+alpha)*(2*i+alpha-2);
        double d = 2*(i+alpha-1)*(i-1)*(2*i+alpha);
        u.push_back((b*((c)*(2*x-y)+(alpha*alpha)*y)*u[n] - d*(y*y)*u[n-1])/(a));
    }
}

void Poly_1D::CalcRJacobi(const int &p, const double &x, const double &y, const double &alpha, std::vector<double> &u)
{
    //
    
    u.push_back(0);
    if (p == 0) { return; }
    for (int n = 1; n <= p; n++)
    {
        double i = n;
        std::vector<double> P_i;
        std::vector<double> P_i_1;
        CalcJacobi(i, x, y, alpha, P_i);
        CalcJacobi(i-1, x, y, alpha, P_i_1);
        u.push_back((-1*(i/(2*i+alpha)))*(P_i[P_i.size()-1] + y*P_i_1[P_i_1.size()-1]));
    }
}

void Poly_1D::CalcScaledJacobiDerivative(const int &p, const double &alpha,
                                      const double &x, const double &t,
                                      std::vector<double> &u, std::vector<double> &dudx, std::vector<double> &dudt)
{
   MFEM_ASSERT(p >= 0, "Polynomial order must be zero or larger");

   u.push_back(1.0);
   dudx.push_back(0.0);
   dudt.push_back(0.0);
   if (p >= 1)
   {
      u.push_back( (2.0 + alpha) * x - t);
      dudx.push_back(2.0 + alpha);
      dudt.push_back(-1.0);
   }
   for (int i = 2; i <= p; i++)
   {
      double a = 2.0 * i * (alpha + i) * (2.0 * i + alpha - 2.0);
      double b = 2.0 * i + alpha - 1.0;
      double c = (2.0 * i + alpha) * (2.0 * i + alpha - 2.0);
      double d = 2.0 * (alpha + i - 1.0) * (i - 1) * (2.0 * i + alpha);
      u.push_back((b * (c * (2.0 * x - t) + alpha * alpha * t) * u[i - 1]
              - d * t * t * u[i - 2]) / a);
      dudx.push_back((b * ((c * (2.0 * x - t) + alpha * alpha * t) * dudx[i - 1] +
                      2.0 * c * u[i - 1])
                 - d * t * t * dudx[i - 2]) / a);
      dudt.push_back( (b * ((c * (2.0 * x - t) + alpha * alpha * t) * dudt[i - 1] +
                      (alpha * alpha - c) * u[i - 1])
                 - d * t * t * dudt[i - 2] - 2.0 * d * t * u[i - 2]) / a);
   }
}

void Poly_1D::CalcIntJacobi(const int &p, const double &x, const double &y, const double &alpha, std::vector<double> &u)
{
    
    // use the recursive definition for [0,1]:
    // L_i = a*P_i+b*y*P_(i-1)-c*y^2*P_(i-2)
    // ** Note the scaled and Shifted polynomials are otained when y != 1
    
    u.push_back(x);
    if (p==1) {
        return;
     }
    for (int n = 1; n < p; n++)
    {
        double i = n+1;
        std::vector<double> P_i;
        std::vector<double> P_i_1;
        std::vector<double> P_i_2;
        CalcJacobi(i, x, y, alpha, P_i);
        CalcJacobi(i-1, x, y, alpha, P_i_1);
        CalcJacobi(i-2, x, y, alpha, P_i_2);
        double a = (i+alpha)/((2*i+alpha-1)*(2*i+alpha));
        double b = (alpha)/((2*i+alpha-2)*(2*i+alpha));
        double c = (i-1)/((2*i+alpha-2)*(2*i+alpha-1));
        double nextval = a*P_i[P_i.size()-1]+b*y*P_i_1[P_i_1.size()-1]-c*(y*y)*P_i_2[P_i_2.size()-1];
        u.push_back(nextval);
    }
}

void Poly_1D::CalcChebyshev(const int p, const double x, double *u)
{
   // recursive definition, z in [-1,1]
   // T_0(z) = 1,  T_1(z) = z
   // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
   double z;
   u[0] = 1.;
   if (p == 0) { return; }
   u[1] = z = 2.*x - 1.;
   for (int n = 1; n < p; n++)
   {
      u[n+1] = 2*z*u[n] - u[n-1];
   }
}

void Poly_1D::CalcChebyshev(const int p, const double x, double *u, double *d)
{
   // recursive definition, z in [-1,1]
   // T_0(z) = 1,  T_1(z) = z
   // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
   // T'_n(z) = n*U_{n-1}(z)
   // U_0(z) = 1  U_1(z) = 2*z
   // U_{n+1}(z) = 2*z*U_n(z) - U_{n-1}(z)
   // U_n(z) = z*U_{n-1}(z) + T_n(z) = z*T'_n(z)/n + T_n(z)
   // T'_{n+1}(z) = (n + 1)*(z*T'_n(z)/n + T_n(z))
   double z;
   u[0] = 1.;
   d[0] = 0.;
   if (p == 0) { return; }
   u[1] = z = 2.*x - 1.;
   d[1] = 2.;
   for (int n = 1; n < p; n++)
   {
      u[n+1] = 2*z*u[n] - u[n-1];
      d[n+1] = (n + 1)*(z*d[n]/n + 2*u[n]);
   }
}

void Poly_1D::CalcChebyshev(const int p, const double x, double *u, double *d,
                            double *dd)
{
   // recursive definition, z in [-1,1]
   // T_0(z) = 1,  T_1(z) = z
   // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
   // T'_n(z) = n*U_{n-1}(z)
   // U_0(z) = 1  U_1(z) = 2*z
   // U_{n+1}(z) = 2*z*U_n(z) - U_{n-1}(z)
   // U_n(z) = z*U_{n-1}(z) + T_n(z) = z*T'_n(z)/n + T_n(z)
   // T'_{n+1}(z) = (n + 1)*(z*T'_n(z)/n + T_n(z))
   // T''_{n+1}(z) = (n + 1)*(2*(n + 1)*T'_n(z) + z*T''_n(z)) / n
   double z;
   u[0] = 1.;
   d[0] = 0.;
   dd[0]= 0.;
   if (p == 0) { return; }
   u[1] = z = 2.*x - 1.;
   d[1] = 2.;
   dd[1] = 0;
   for (int n = 1; n < p; n++)
   {
      u[n+1] = 2*z*u[n] - u[n-1];
      d[n+1] = (n + 1)*(z*d[n]/n + 2*u[n]);
      dd[n+1] = (n + 1)*(2.*(n + 1)*d[n] + z*dd[n])/n;
   }
}

const double *Poly_1D::GetPoints(const int p, const int btype)
{
   BasisType::Check(btype);
   const int qtype = BasisType::GetQuadrature1D(btype);

   if (qtype == Quadrature1D::Invalid) { return NULL; }

   if (points_container.find(btype) == points_container.end())
   {
      points_container[btype] = new Array<double*>(h_mt);
   }
   Array<double*> &pts = *points_container[btype];
   if (pts.Size() <= p)
   {
      pts.SetSize(p + 1, NULL);
   }
   if (pts[p] == NULL)
   {
      pts[p] = new double[p + 1];
      quad_func.GivePolyPoints(p+1, pts[p], qtype);
   }
   return pts[p];
}

Poly_1D::Basis &Poly_1D::GetBasis(const int p, const int btype)
{
   BasisType::Check(btype);

   if ( bases_container.find(btype) == bases_container.end() )
   {
      // we haven't been asked for basis or points of this type yet
      bases_container[btype] = new Array<Basis*>(h_mt);
   }
   Array<Basis*> &bases = *bases_container[btype];
   if (bases.Size() <= p)
   {
      bases.SetSize(p + 1, NULL);
   }
   if (bases[p] == NULL)
   {
      EvalType etype = (btype == BasisType::Positive) ? Positive : Barycentric;
      bases[p] = new Basis(p, GetPoints(p, btype), etype);
   }
   return *bases[p];
}

Poly_1D::~Poly_1D()
{
   for (PointsMap::iterator it = points_container.begin();
        it != points_container.end() ; ++it)
   {
      Array<double*>& pts = *it->second;
      for ( int i = 0 ; i < pts.Size() ; ++i )
      {
         delete [] pts[i];
      }
      delete it->second;
   }

   for (BasisMap::iterator it = bases_container.begin();
        it != bases_container.end() ; ++it)
   {
      Array<Basis*>& bases = *it->second;
      for ( int i = 0 ; i < bases.Size() ; ++i )
      {
         delete bases[i];
      }
      delete it->second;
   }
}

Array2D<int> Poly_1D::binom;
Poly_1D poly1d;


TensorBasisElement::TensorBasisElement(const int dims, const int p,
                                       const int btype, const DofMapType dmtype)
   : b_type(btype),
     basis1d(poly1d.GetBasis(p, b_type))
{
   if (dmtype == H1_DOF_MAP || dmtype == Sr_DOF_MAP)
   {
      switch (dims)
      {
         case 1:
         {
            dof_map.SetSize(p + 1);
            dof_map[0] = 0;
            dof_map[p] = 1;
            for (int i = 1; i < p; i++)
            {
               dof_map[i] = i+1;
            }
            break;
         }
         case 2:
         {
            const int p1 = p + 1;
            dof_map.SetSize(p1*p1);

            // vertices
            dof_map[0 + 0*p1] = 0;
            dof_map[p + 0*p1] = 1;
            dof_map[p + p*p1] = 2;
            dof_map[0 + p*p1] = 3;

            // edges
            int o = 4;
            for (int i = 1; i < p; i++)
            {
               dof_map[i + 0*p1] = o++;
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[p + i*p1] = o++;
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[(p-i) + p*p1] = o++;
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[0 + (p-i)*p1] = o++;
            }

            // interior
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[i + j*p1] = o++;
               }
            }
            break;
         }
         case 3:
         {
            const int p1 = p + 1;
            dof_map.SetSize(p1*p1*p1);

            // vertices
            dof_map[0 + (0 + 0*p1)*p1] = 0;
            dof_map[p + (0 + 0*p1)*p1] = 1;
            dof_map[p + (p + 0*p1)*p1] = 2;
            dof_map[0 + (p + 0*p1)*p1] = 3;
            dof_map[0 + (0 + p*p1)*p1] = 4;
            dof_map[p + (0 + p*p1)*p1] = 5;
            dof_map[p + (p + p*p1)*p1] = 6;
            dof_map[0 + (p + p*p1)*p1] = 7;

            // edges (see Hexahedron::edges in mesh/hexahedron.cpp).
            // edges (see Constants<Geometry::CUBE>::Edges in fem/geom.cpp).
            int o = 8;
            for (int i = 1; i < p; i++)
            {
               dof_map[i + (0 + 0*p1)*p1] = o++;   // (0,1)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[p + (i + 0*p1)*p1] = o++;   // (1,2)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[i + (p + 0*p1)*p1] = o++;   // (3,2)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[0 + (i + 0*p1)*p1] = o++;   // (0,3)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[i + (0 + p*p1)*p1] = o++;   // (4,5)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[p + (i + p*p1)*p1] = o++;   // (5,6)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[i + (p + p*p1)*p1] = o++;   // (7,6)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[0 + (i + p*p1)*p1] = o++;   // (4,7)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[0 + (0 + i*p1)*p1] = o++;   // (0,4)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[p + (0 + i*p1)*p1] = o++;   // (1,5)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[p + (p + i*p1)*p1] = o++;   // (2,6)
            }
            for (int i = 1; i < p; i++)
            {
               dof_map[0 + (p + i*p1)*p1] = o++;   // (3,7)
            }

            // faces (see Mesh::GenerateFaces in mesh/mesh.cpp)
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[i + ((p-j) + 0*p1)*p1] = o++;   // (3,2,1,0)
               }
            }
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[i + (0 + j*p1)*p1] = o++;   // (0,1,5,4)
               }
            }
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[p + (i + j*p1)*p1] = o++;   // (1,2,6,5)
               }
            }
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[(p-i) + (p + j*p1)*p1] = o++;   // (2,3,7,6)
               }
            }
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[0 + ((p-i) + j*p1)*p1] = o++;   // (3,0,4,7)
               }
            }
            for (int j = 1; j < p; j++)
            {
               for (int i = 1; i < p; i++)
               {
                  dof_map[i + (j + p*p1)*p1] = o++;   // (4,5,6,7)
               }
            }

            // interior
            for (int k = 1; k < p; k++)
            {
               for (int j = 1; j < p; j++)
               {
                  for (int i = 1; i < p; i++)
                  {
                     dof_map[i + (j + k*p1)*p1] = o++;
                  }
               }
            }
            break;
         }
         default:
            MFEM_ABORT("invalid dimension: " << dims);
            break;
      }
   }
   else if (dmtype == L2_DOF_MAP)
   {
      // leave dof_map empty, indicating that the dofs are ordered
      // lexicographically, i.e. the dof_map is identity
   }
   else
   {
      MFEM_ABORT("invalid DofMapType: " << dmtype);
   }
}


NodalTensorFiniteElement::NodalTensorFiniteElement(const int dims,
                                                   const int p,
                                                   const int btype,
                                                   const DofMapType dmtype)
   : NodalFiniteElement(dims, GetTensorProductGeometry(dims), Pow(p + 1, dims),
                        p, dims > 1 ? FunctionSpace::Qk : FunctionSpace::Pk),
     TensorBasisElement(dims, p, VerifyNodal(btype), dmtype)
{
   lex_ordering = dof_map;
}


PositiveTensorFiniteElement::PositiveTensorFiniteElement(
   const int dims, const int p, const DofMapType dmtype)
   : PositiveFiniteElement(dims, GetTensorProductGeometry(dims),
                           Pow(p + 1, dims), p,
                           dims > 1 ? FunctionSpace::Qk : FunctionSpace::Pk),
     TensorBasisElement(dims, p, BasisType::Positive, dmtype) { }

VectorTensorFiniteElement::VectorTensorFiniteElement(const int dims,
                                                     const int d,
                                                     const int p,
                                                     const int cbtype,
                                                     const int obtype,
                                                     const int M,
                                                     const DofMapType dmtype)
   : VectorFiniteElement(dims, GetTensorProductGeometry(dims), d,
                         p, M, FunctionSpace::Qk),
     TensorBasisElement(dims, p, VerifyNodal(cbtype), dmtype),
     cbasis1d(poly1d.GetBasis(p, VerifyClosed(cbtype))),
     obasis1d(poly1d.GetBasis(p - 1, VerifyOpen(obtype))) { }

H1_SegmentElement::H1_SegmentElement(const int p, const int btype)
   : NodalTensorFiniteElement(1, p, VerifyClosed(btype), H1_DOF_MAP)
{
   const double *cp = poly1d.ClosedPoints(p, b_type);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p+1);
   dshape_x.SetSize(p+1);
   d2shape_x.SetSize(p+1);
#endif

   Nodes.IntPoint(0).x = cp[0];
   Nodes.IntPoint(1).x = cp[p];
   for (int i = 1; i < p; i++)
   {
      Nodes.IntPoint(i+1).x = cp[i];
   }
}

void H1_SegmentElement::CalcShape(const IntegrationPoint &ip,
                                  Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1);
#endif

   basis1d.Eval(ip.x, shape_x);

   shape(0) = shape_x(0);
   shape(1) = shape_x(p);
   for (int i = 1; i < p; i++)
   {
      shape(i+1) = shape_x(i);
   }
}

void H1_SegmentElement::CalcDShape(const IntegrationPoint &ip,
                                   DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), dshape_x(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x);

   dshape(0,0) = dshape_x(0);
   dshape(1,0) = dshape_x(p);
   for (int i = 1; i < p; i++)
   {
      dshape(i+1,0) = dshape_x(i);
   }
}

void H1_SegmentElement::CalcHessian(const IntegrationPoint &ip,
                                    DenseMatrix &Hessian) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), dshape_x(p+1), d2shape_x(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x, d2shape_x);

   Hessian(0,0) = d2shape_x(0);
   Hessian(1,0) = d2shape_x(p);
   for (int i = 1; i < p; i++)
   {
      Hessian(i+1,0) = d2shape_x(i);
   }
}

void H1_SegmentElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *cp = poly1d.ClosedPoints(p, b_type);

   switch (vertex)
   {
      case 0:
         dofs(0) = poly1d.CalcDelta(p, (1.0 - cp[0]));
         dofs(1) = poly1d.CalcDelta(p, (1.0 - cp[p]));
         for (int i = 1; i < p; i++)
         {
            dofs(i+1) = poly1d.CalcDelta(p, (1.0 - cp[i]));
         }
         break;

      case 1:
         dofs(0) = poly1d.CalcDelta(p, cp[0]);
         dofs(1) = poly1d.CalcDelta(p, cp[p]);
         for (int i = 1; i < p; i++)
         {
            dofs(i+1) = poly1d.CalcDelta(p, cp[i]);
         }
         break;
   }
}


H1_QuadrilateralElement::H1_QuadrilateralElement(const int p, const int btype)
   : NodalTensorFiniteElement(2, p, VerifyClosed(btype), H1_DOF_MAP)
{
   const double *cp = poly1d.ClosedPoints(p, b_type);

#ifndef MFEM_THREAD_SAFE
   const int p1 = p + 1;

   shape_x.SetSize(p1);
   shape_y.SetSize(p1);
   dshape_x.SetSize(p1);
   dshape_y.SetSize(p1);
   d2shape_x.SetSize(p1);
   d2shape_y.SetSize(p1);
#endif

   int o = 0;
   for (int j = 0; j <= p; j++)
   {
      for (int i = 0; i <= p; i++)
      {
         Nodes.IntPoint(dof_map[o++]).Set2(cp[i], cp[j]);
      }
   }
}

void H1_QuadrilateralElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   basis1d.Eval(ip.x, shape_x);
   basis1d.Eval(ip.y, shape_y);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         shape(dof_map[o++]) = shape_x(i)*shape_y(j);
      }
}

void H1_QuadrilateralElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), dshape_x(p+1), dshape_y(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y);

   for (int o = 0, j = 0; j <= p; j++)
   {
      for (int i = 0; i <= p; i++)
      {
         dshape(dof_map[o],0) = dshape_x(i)* shape_y(j);
         dshape(dof_map[o],1) =  shape_x(i)*dshape_y(j);  o++;
      }
   }
}

void H1_QuadrilateralElement::CalcHessian(const IntegrationPoint &ip,
                                          DenseMatrix &Hessian) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), dshape_x(p+1), dshape_y(p+1),
          d2shape_x(p+1), d2shape_y(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x, d2shape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y, d2shape_y);

   for (int o = 0, j = 0; j <= p; j++)
   {
      for (int i = 0; i <= p; i++)
      {
         Hessian(dof_map[o],0) = d2shape_x(i)*  shape_y(j);
         Hessian(dof_map[o],1) =  dshape_x(i)* dshape_y(j);
         Hessian(dof_map[o],2) =   shape_x(i)*d2shape_y(j);  o++;
      }
   }
}

void H1_QuadrilateralElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *cp = poly1d.ClosedPoints(p, b_type);

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   for (int i = 0; i <= p; i++)
   {
      shape_x(i) = poly1d.CalcDelta(p, (1.0 - cp[i]));
      shape_y(i) = poly1d.CalcDelta(p, cp[i]);
   }

   switch (vertex)
   {
      case 0:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs(dof_map[o++]) = shape_x(i)*shape_x(j);
            }
         break;
      case 1:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs(dof_map[o++]) = shape_y(i)*shape_x(j);
            }
         break;
      case 2:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs(dof_map[o++]) = shape_y(i)*shape_y(j);
            }
         break;
      case 3:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs(dof_map[o++]) = shape_x(i)*shape_y(j);
            }
         break;
   }
}


H1_HexahedronElement::H1_HexahedronElement(const int p, const int btype)
   : NodalTensorFiniteElement(3, p, VerifyClosed(btype), H1_DOF_MAP)
{
   const double *cp = poly1d.ClosedPoints(p, b_type);

#ifndef MFEM_THREAD_SAFE
   const int p1 = p + 1;

   shape_x.SetSize(p1);
   shape_y.SetSize(p1);
   shape_z.SetSize(p1);
   dshape_x.SetSize(p1);
   dshape_y.SetSize(p1);
   dshape_z.SetSize(p1);
   d2shape_x.SetSize(p1);
   d2shape_y.SetSize(p1);
   d2shape_z.SetSize(p1);
#endif

   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            Nodes.IntPoint(dof_map[o++]).Set3(cp[i], cp[j], cp[k]);
         }
}

void H1_HexahedronElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), shape_z(p+1);
#endif

   basis1d.Eval(ip.x, shape_x);
   basis1d.Eval(ip.y, shape_y);
   basis1d.Eval(ip.z, shape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            shape(dof_map[o++]) = shape_x(i)*shape_y(j)*shape_z(k);
         }
}

void H1_HexahedronElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1),  shape_y(p+1),  shape_z(p+1);
   Vector dshape_x(p+1), dshape_y(p+1), dshape_z(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y);
   basis1d.Eval(ip.z, shape_z, dshape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dshape(dof_map[o],0) = dshape_x(i)* shape_y(j)* shape_z(k);
            dshape(dof_map[o],1) =  shape_x(i)*dshape_y(j)* shape_z(k);
            dshape(dof_map[o],2) =  shape_x(i)* shape_y(j)*dshape_z(k);  o++;
         }
}

void H1_HexahedronElement::CalcHessian(const IntegrationPoint &ip,
                                       DenseMatrix &Hessian) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1),  shape_y(p+1),  shape_z(p+1);
   Vector dshape_x(p+1), dshape_y(p+1), dshape_z(p+1);
   Vector d2shape_x(p+1), d2shape_y(p+1), d2shape_z(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x, d2shape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y, d2shape_y);
   basis1d.Eval(ip.z, shape_z, dshape_z, d2shape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            Hessian(dof_map[o],0) = d2shape_x(i)*  shape_y(j)*  shape_z(k);
            Hessian(dof_map[o],1) =  dshape_x(i)* dshape_y(j)*  shape_z(k);
            Hessian(dof_map[o],2) =  dshape_x(i)*  shape_y(j)* dshape_z(k);
            Hessian(dof_map[o],3) =   shape_x(i)*d2shape_y(j)*  shape_z(k);
            Hessian(dof_map[o],4) =   shape_x(i)* dshape_y(j)* dshape_z(k);
            Hessian(dof_map[o],5) =   shape_x(i)*  shape_y(j)*d2shape_z(k);
            o++;
         }
}

void H1_HexahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *cp = poly1d.ClosedPoints(p,b_type);

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   for (int i = 0; i <= p; i++)
   {
      shape_x(i) = poly1d.CalcDelta(p, (1.0 - cp[i]));
      shape_y(i) = poly1d.CalcDelta(p, cp[i]);
   }

   switch (vertex)
   {
      case 0:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_x(i)*shape_x(j)*shape_x(k);
               }
         break;
      case 1:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_y(i)*shape_x(j)*shape_x(k);
               }
         break;
      case 2:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_y(i)*shape_y(j)*shape_x(k);
               }
         break;
      case 3:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_x(i)*shape_y(j)*shape_x(k);
               }
         break;
      case 4:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_x(i)*shape_x(j)*shape_y(k);
               }
         break;
      case 5:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_y(i)*shape_x(j)*shape_y(k);
               }
         break;
      case 6:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_y(i)*shape_y(j)*shape_y(k);
               }
         break;
      case 7:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs(dof_map[o++]) = shape_x(i)*shape_y(j)*shape_y(k);
               }
         break;
   }
}


H1Pos_SegmentElement::H1Pos_SegmentElement(const int p)
   : PositiveTensorFiniteElement(1, p, H1_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   // thread private versions; see class header.
   shape_x.SetSize(p+1);
   dshape_x.SetSize(p+1);
#endif

   // Endpoints need to be first in the list, so reorder them.
   Nodes.IntPoint(0).x = 0.0;
   Nodes.IntPoint(1).x = 1.0;
   for (int i = 1; i < p; i++)
   {
      Nodes.IntPoint(i+1).x = double(i)/p;
   }
}

void H1Pos_SegmentElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData() );

   // Endpoints need to be first in the list, so reorder them.
   shape(0) = shape_x(0);
   shape(1) = shape_x(p);
   for (int i = 1; i < p; i++)
   {
      shape(i+1) = shape_x(i);
   }
}

void H1Pos_SegmentElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), dshape_x(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData(), dshape_x.GetData() );

   // Endpoints need to be first in the list, so reorder them.
   dshape(0,0) = dshape_x(0);
   dshape(1,0) = dshape_x(p);
   for (int i = 1; i < p; i++)
   {
      dshape(i+1,0) = dshape_x(i);
   }
}

void H1Pos_SegmentElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   dofs[vertex] = 1.0;
}


H1Pos_QuadrilateralElement::H1Pos_QuadrilateralElement(const int p)
   : PositiveTensorFiniteElement(2, p, H1_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   const int p1 = p + 1;

   shape_x.SetSize(p1);
   shape_y.SetSize(p1);
   dshape_x.SetSize(p1);
   dshape_y.SetSize(p1);
#endif

   int o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         Nodes.IntPoint(dof_map[o++]).Set2(double(i)/p, double(j)/p);
      }
}

void H1Pos_QuadrilateralElement::CalcShape(const IntegrationPoint &ip,
                                           Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData() );
   Poly_1D::CalcBernstein(p, ip.y, shape_y.GetData() );

   // Reorder so that vertices are at the beginning of the list
   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         shape(dof_map[o++]) = shape_x(i)*shape_y(j);
      }
}

void H1Pos_QuadrilateralElement::CalcDShape(const IntegrationPoint &ip,
                                            DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), dshape_x(p+1), dshape_y(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData(), dshape_x.GetData() );
   Poly_1D::CalcBernstein(p, ip.y, shape_y.GetData(), dshape_y.GetData() );

   // Reorder so that vertices are at the beginning of the list
   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         dshape(dof_map[o],0) = dshape_x(i)* shape_y(j);
         dshape(dof_map[o],1) =  shape_x(i)*dshape_y(j);  o++;
      }
}

void H1Pos_QuadrilateralElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   dofs[vertex] = 1.0;
}


H1Pos_HexahedronElement::H1Pos_HexahedronElement(const int p)
   : PositiveTensorFiniteElement(3, p, H1_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   const int p1 = p + 1;

   shape_x.SetSize(p1);
   shape_y.SetSize(p1);
   shape_z.SetSize(p1);
   dshape_x.SetSize(p1);
   dshape_y.SetSize(p1);
   dshape_z.SetSize(p1);
#endif

   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
            Nodes.IntPoint(dof_map[o++]).Set3(double(i)/p, double(j)/p,
                                              double(k)/p);
}

void H1Pos_HexahedronElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), shape_z(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData() );
   Poly_1D::CalcBernstein(p, ip.y, shape_y.GetData() );
   Poly_1D::CalcBernstein(p, ip.z, shape_z.GetData() );

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            shape(dof_map[o++]) = shape_x(i)*shape_y(j)*shape_z(k);
         }
}

void H1Pos_HexahedronElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1),  shape_y(p+1),  shape_z(p+1);
   Vector dshape_x(p+1), dshape_y(p+1), dshape_z(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x.GetData(), dshape_x.GetData() );
   Poly_1D::CalcBernstein(p, ip.y, shape_y.GetData(), dshape_y.GetData() );
   Poly_1D::CalcBernstein(p, ip.z, shape_z.GetData(), dshape_z.GetData() );

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dshape(dof_map[o],0) = dshape_x(i)* shape_y(j)* shape_z(k);
            dshape(dof_map[o],1) =  shape_x(i)*dshape_y(j)* shape_z(k);
            dshape(dof_map[o],2) =  shape_x(i)* shape_y(j)*dshape_z(k);  o++;
         }
}

void H1Pos_HexahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   dofs[vertex] = 1.0;
}


H1_TriangleElement::H1_TriangleElement(const int p, const int btype)
   : NodalFiniteElement(2, Geometry::TRIANGLE, ((p + 1)*(p + 2))/2, p,
                        FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyNodal(VerifyClosed(btype)));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof, (dim * (dim + 1)) / 2 );
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
#endif

   int p2p3 = 2*p + 3;
   auto idx = [p2p3](int i, int j) { return ((p2p3-j)*j)/2+i; };
   lex_ordering.SetSize(dof);

   // vertices
   lex_ordering[idx(0,0)] = 0;
   Nodes.IntPoint(0).Set2(cp[0], cp[0]);
   lex_ordering[idx(p,0)] = 1;
   Nodes.IntPoint(1).Set2(cp[p], cp[0]);
   lex_ordering[idx(0,p)] = 2;
   Nodes.IntPoint(2).Set2(cp[0], cp[p]);

   // edges
   int o = 3;
   for (int i = 1; i < p; i++)
   {
      lex_ordering[idx(i,0)] = o;
      Nodes.IntPoint(o++).Set2(cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)
   {
      lex_ordering[idx(p-i,i)] = o;
      Nodes.IntPoint(o++).Set2(cp[p-i], cp[i]);
   }
   for (int i = 1; i < p; i++)
   {
      lex_ordering[idx(0,p-i)] = o;
      Nodes.IntPoint(o++).Set2(cp[0], cp[p-i]);
   }

   // interior
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)
      {
         const double w = cp[i] + cp[j] + cp[p-i-j];
         lex_ordering[idx(i,j)] = o;
         Nodes.IntPoint(o++).Set2(cp[i]/w, cp[j]/w);
      }

   DenseMatrix T(dof);
   for (int k = 0; k < dof; k++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(k);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);

      o = 0;
      for (int j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            T(o++, k) = shape_x(i)*shape_y(j)*shape_l(p-i-j);
         }
   }

   Ti.Factor(T);
   // mfem::out << "H1_TriangleElement(" << p << ") : "; Ti.TestInversion();
}

void H1_TriangleElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1), u(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         u(o++) = shape_x(i)*shape_y(j)*shape_l(p-i-j);
      }

   Ti.Mult(u, shape);
}

void H1_TriangleElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_l(p + 1);
   DenseMatrix du(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l, dshape_l);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         int k = p - i - j;
         du(o,0) = ((dshape_x(i)* shape_l(k)) -
                    ( shape_x(i)*dshape_l(k)))*shape_y(j);
         du(o,1) = ((dshape_y(j)* shape_l(k)) -
                    ( shape_y(j)*dshape_l(k)))*shape_x(i);
         o++;
      }

   Ti.Mult(du, dshape);
}

void H1_TriangleElement::CalcHessian(const IntegrationPoint &ip,
                                     DenseMatrix &ddshape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector   shape_x(p + 1),   shape_y(p + 1),   shape_l(p + 1);
   Vector  dshape_x(p + 1),  dshape_y(p + 1),  dshape_l(p + 1);
   Vector ddshape_x(p + 1), ddshape_y(p + 1), ddshape_l(p + 1);
   DenseMatrix ddu(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x, ddshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y, ddshape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l, dshape_l, ddshape_l);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         int k = p - i - j;
         // u_xx, u_xy, u_yy
         ddu(o,0) = ((ddshape_x(i) * shape_l(k)) - 2. * (dshape_x(i) * dshape_l(k)) +
                     (shape_x(i) * ddshape_l(k))) * shape_y(j);
         ddu(o,1) = (((shape_x(i) * ddshape_l(k)) - dshape_x(i) * dshape_l(k)) * shape_y(
                        j)) + (((dshape_x(i) * shape_l(k)) - (shape_x(i) * dshape_l(k))) * dshape_y(j));
         ddu(o,2) = ((ddshape_y(j) * shape_l(k)) - 2. * (dshape_y(j) * dshape_l(k)) +
                     (shape_y(j) * ddshape_l(k))) * shape_x(i);
         o++;
      }

   Ti.Mult(ddu, ddshape);
}


H1_TetrahedronElement::H1_TetrahedronElement(const int p, const int btype)
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, ((p + 1)*(p + 2)*(p + 3))/6,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyNodal(VerifyClosed(btype)));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof, (dim * (dim + 1)) / 2);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
#endif

   auto tri = [](int k) { return (k*(k + 1))/2; };
   auto tet = [](int k) { return (k*(k + 1)*(k + 2))/6; };
   int ndof = tet(p+1);
   auto idx = [tri, tet, p, ndof](int i, int j, int k)
   {
      return ndof - tet(p - k) - tri(p + 1 - k - j) + i;
   };

   lex_ordering.SetSize(dof);

   // vertices
   lex_ordering[idx(0,0,0)] = 0;
   Nodes.IntPoint(0).Set3(cp[0], cp[0], cp[0]);
   lex_ordering[idx(p,0,0)] = 1;
   Nodes.IntPoint(1).Set3(cp[p], cp[0], cp[0]);
   lex_ordering[idx(0,p,0)] = 2;
   Nodes.IntPoint(2).Set3(cp[0], cp[p], cp[0]);
   lex_ordering[idx(0,0,p)] = 3;
   Nodes.IntPoint(3).Set3(cp[0], cp[0], cp[p]);

   // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
   int o = 4;
   for (int i = 1; i < p; i++)  // (0,1)
   {
      lex_ordering[idx(i,0,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,2)
   {
      lex_ordering[idx(0,i,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,3)
   {
      lex_ordering[idx(0,0,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (1,2)
   {
      lex_ordering[idx(p-i,i,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[p-i], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,3)
   {
      lex_ordering[idx(p-i,0,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[p-i], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (2,3)
   {
      lex_ordering[idx(0,p-i,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[p-i], cp[i]);
   }

   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (1,2,3)
      {
         lex_ordering[idx(p-i-j,i,j)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[p-i-j]/w, cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,3,2)
      {
         lex_ordering[idx(0,j,i)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[0], cp[j]/w, cp[i]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,1,3)
      {
         lex_ordering[idx(i,0,j)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,2,1)
      {
         lex_ordering[idx(j,i,0)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[j]/w, cp[i]/w, cp[0]);
      }

   // interior
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)
         {
            lex_ordering[idx(i,j,k)] = o;
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set3(cp[i]/w, cp[j]/w, cp[k]/w);
         }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);

      o = 0;
      for (int k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               T(o++, m) = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
            }
   }
    
//    ofstream logfile;
//    logfile.open ("logfile.txt");
//    for (int i =0; i<T.Size(); i++) {
//        for (int j = 0; j<T.Size(); j++) {
//            logfile << "T value @ " << i << "," << j << "= " << T(i,j) << std::endl;
//        }
//    }
//    logfile.close();

   Ti.Factor(T);
   // mfem::out << "H1_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void H1_TetrahedronElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   Vector u(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
         }

   Ti.Mult(u, shape);
}

void H1_TetrahedronElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   DenseMatrix du(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            int l = p - i - j - k;
            du(o,0) = ((dshape_x(i)* shape_l(l)) -
                       ( shape_x(i)*dshape_l(l)))*shape_y(j)*shape_z(k);
            du(o,1) = ((dshape_y(j)* shape_l(l)) -
                       ( shape_y(j)*dshape_l(l)))*shape_x(i)*shape_z(k);
            du(o,2) = ((dshape_z(k)* shape_l(l)) -
                       ( shape_z(k)*dshape_l(l)))*shape_x(i)*shape_y(j);
            o++;
         }

   Ti.Mult(du, dshape);
}

void H1_TetrahedronElement::CalcHessian(const IntegrationPoint &ip,
                                        DenseMatrix &ddshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector   shape_x(p + 1),   shape_y(p + 1),   shape_z(p + 1),   shape_l(p + 1);
   Vector  dshape_x(p + 1),  dshape_y(p + 1),  dshape_z(p + 1),  dshape_l(p + 1);
   Vector ddshape_x(p + 1), ddshape_y(p + 1), ddshape_z(p + 1), ddshape_l(p + 1);
   DenseMatrix ddu(dof, ((dim + 1) * dim) / 2);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x, ddshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y, ddshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z, ddshape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l, ddshape_l);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            // u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
            int l = p - i - j - k;
            ddu(o,0) = ((ddshape_x(i) * shape_l(l)) - 2. * (dshape_x(i) * dshape_l(l)) +
                        (shape_x(i) * ddshape_l(l))) * shape_y(j) * shape_z(k);
            ddu(o,1) = ((dshape_y(j) * ((dshape_x(i) * shape_l(l)) -
                                        (shape_x(i) * dshape_l(l)))) +
                        (shape_y(j) * ((ddshape_l(l) * shape_x(i)) -
                                       (dshape_x(i) * dshape_l(l)))))* shape_z(k);
            ddu(o,2) = ((dshape_z(k) * ((dshape_x(i) * shape_l(l)) -
                                        (shape_x(i) * dshape_l(l)))) +
                        (shape_z(k) * ((ddshape_l(l) * shape_x(i)) -
                                       (dshape_x(i) * dshape_l(l)))))* shape_y(j);
            ddu(o,3) = ((ddshape_y(j) * shape_l(l)) - 2. * (dshape_y(j) * dshape_l(l)) +
                        (shape_y(j) * ddshape_l(l))) * shape_x(i) * shape_z(k);
            ddu(o,4) = ((dshape_z(k) * ((dshape_y(j) * shape_l(l)) -
                                        (shape_y(j)*dshape_l(l))) ) +
                        (shape_z(k)* ((ddshape_l(l)*shape_y(j)) -
                                      (dshape_y(j) * dshape_l(l)) ) ) )* shape_x(i);
            ddu(o,5) = ((ddshape_z(k) * shape_l(l)) - 2. * (dshape_z(k) * dshape_l(l)) +
                        (shape_z(k) * ddshape_l(l))) * shape_y(j) * shape_x(i);
            o++;
         }
   Ti.Mult(ddu, ddshape);
}

H1_TetrahedronElement_Fuentes::H1_TetrahedronElement_Fuentes(const int p, const int btype)
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, ((p + 1)*(p + 2)*(p + 3))/6,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyNodal(VerifyClosed(btype)));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof, (dim * (dim + 1)) / 2);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
#endif

   auto tri = [](int k) { return (k*(k + 1))/2; };
   auto tet = [](int k) { return (k*(k + 1)*(k + 2))/6; };
   int ndof = tet(p+1);
   auto idx = [tri, tet, p, ndof](int i, int j, int k)
   {
      return ndof - tet(p - k) - tri(p + 1 - k - j) + i;
   };

   lex_ordering.SetSize(dof);

   // vertices
   lex_ordering[idx(0,0,0)] = 0;
   Nodes.IntPoint(0).Set3(cp[0], cp[0], cp[0]);
   lex_ordering[idx(p,0,0)] = 1;
   Nodes.IntPoint(1).Set3(cp[p], cp[0], cp[0]);
   lex_ordering[idx(0,p,0)] = 2;
   Nodes.IntPoint(2).Set3(cp[0], cp[p], cp[0]);
   lex_ordering[idx(0,0,p)] = 3;
   Nodes.IntPoint(3).Set3(cp[0], cp[0], cp[p]);

   // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
   int o = 4;
   for (int i = 1; i < p; i++)  // (0,1)
   {
      lex_ordering[idx(i,0,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,2)
   {
      lex_ordering[idx(0,i,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,3)
   {
      lex_ordering[idx(0,0,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (1,2)
   {
      lex_ordering[idx(p-i,i,0)] = o;
      Nodes.IntPoint(o++).Set3(cp[p-i], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,3)
   {
      lex_ordering[idx(p-i,0,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[p-i], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (2,3)
   {
      lex_ordering[idx(0,p-i,i)] = o;
      Nodes.IntPoint(o++).Set3(cp[0], cp[p-i], cp[i]);
   }

   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (1,2,3)
      {
         lex_ordering[idx(p-i-j,i,j)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[p-i-j]/w, cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,3,2)
      {
         lex_ordering[idx(0,j,i)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[0], cp[j]/w, cp[i]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,1,3)
      {
         lex_ordering[idx(i,0,j)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,2,1)
      {
         lex_ordering[idx(j,i,0)] = o;
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set3(cp[j]/w, cp[i]/w, cp[0]);
      }

   // interior
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)
         {
            lex_ordering[idx(i,j,k)] = o;
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set3(cp[i]/w, cp[j]/w, cp[k]/w);
         }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
       IntegrationPoint &ip = Nodes.IntPoint(m);
       
       std::vector<double> bary_vector{ip.x, ip.y, ip.z, (1.0 - ip.x - ip.y - ip.z)};
       
       int o = 0;
       
       double La, Lb, Lc, Ld;
       
       //Vertices
       for (int a=0; a<4; a++)
       {
           //compute barycentric function
           T(o++,m) = bary_vector[a];
       }
       
       //Edges
       for(int i=2; i<=p; i++)
       {
           for(int a=0; a<4; a++)
           {
               for(int b=0; b<4; b++)
               {
                   if(a<b)
                   {
                       La = bary_vector[a];
                       Lb = bary_vector[b];
                       
                       // compute polynomials
                       std::vector<double> int_Legendre_i;
                       double x = Lb;
                       double y = La + Lb;
                       poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                       // Add Basis Funcitons
                       T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1];
                       
                   }
               }
           }
       }// end of edges
       
       //Faces
       for(int i=2; i<=p;i++)
       {
           for(int j=1; j<=p;j++)
           {
               for(int a=0; a<4;a++)
               {
                   for(int b=0; b<4;b++)
                   {
                       for(int c=0; c<4;c++)
                       {
                           if((a<b)&&(b<c)&&((i+j)<=p))
                           {
                               La = bary_vector[a];
                               Lb = bary_vector[b];
                               Lc = bary_vector[c];
                               
                               // compute polynomials
                               double x = Lb;
                               double y = La + Lb;
                               std::vector<double> int_Legendre_i;
                               poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                               
                               double alpha = 2.0*i;
                               x = Lc;
                               y = La + Lb + Lc;
                               std::vector<double> int_Jacobi_j;
                               poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                               
                               // Add Basis functions
                               T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1];
                           }
                       }
                   }
               }
           }
       }
       // end of faces
       
      
      // Define Barycentric Coordinates
      La = bary_vector[0];
      Lb = bary_vector[1];
      Lc = bary_vector[2];
      Ld = bary_vector[3];

      for(int i=2; i<=p;i++)
      {
          for(int j=1; j<=p;j++)
          {
              for(int l=1; l<=p;l++)
              {
                  if((i+j+l)<=p)
                  {
                      
                      // compute polynomials
                      double x = Lb;
                      double y = La + Lb;
                      std::vector<double> int_Legendre_i;
                      poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                      
                      double alpha = 2.0*i;
                      x = Lc;
                      y = La + Lb + Lc;
                      std::vector<double> int_Jacobi_j;
                      poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                      
                      alpha = 2.0*(i+j);
                      x = Ld;
                      y = La + Lb + Lc + Ld;
                      std::vector<double> int_Jacobi_l;
                      poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                      
                      // Add Basis Functions
                      T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1];
                      
                  }
              }
          }
      } // end of Facets

   }
    
//    ofstream logfile;
//    logfile.open ("logfile.txt");
//    for (int i =0; i<T.Size(); i++) {
//        for (int j = 0; j<T.Size(); j++) {
//            logfile << "T value @ " << i << "," << j << "= " << T(i,j) << std::endl;
//        }
//    }
//    logfile.close();

   Ti.Factor(T);
   // mfem::out << "H1_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void H1_TetrahedronElement_Fuentes::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   Vector u(dof);
#endif

    
     //compute barycentric coordinates as function of ip
     double L[3];
     ip.Get(L,3);
    
     std::vector<double> bary_vector{L[0], L[1], L[2], (1.0 - L[0]-L[1]-L[2])};

     int o = 0;
     
     double La, Lb, Lc, Ld;
     
    //Vertices
     for (int a=0; a<4; a++)
     {
         //compute barycentric function
         u(o++) = bary_vector[a];
     }
     
    //Edges
     for(int i=2; i<=p; i++)
     {
         for(int a=0; a<4; a++)
         {
             for(int b=0; b<4; b++)
             {
                 if(a<b)
                 {
                     La = bary_vector[a];
                     Lb = bary_vector[b];
                     
                     // compute polynomials
                     std::vector<double> int_Legendre_i;
                     double x = Lb;
                     double y = La + Lb;
                     poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                     
                     // Add Basis Funcitons
                     u(o++) = int_Legendre_i[int_Legendre_i.size()-1];
                 }
             }
         }
     } // end of edges
     
    //Faces
    for(int i=2; i<=p;i++)
    {
        for(int j=1; j<=p;j++)
        {
            for(int a=0; a<4;a++)
            {
                for(int b=0; b<4;b++)
                {
                    for(int c=0; c<4;c++)
                    {
                        if((a<b)&&(b<c)&&((i+j)<=p))
                        {
                            La = bary_vector[a];
                            Lb = bary_vector[b];
                            Lc = bary_vector[c];
                            
                            // compute polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> int_Legendre_i;
                            poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                            
                            double alpha = 2.0*i;
                            x = Lc;
                            y = La + Lb + Lc;
                            std::vector<double> int_Jacobi_j;
                            poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                            
                            // Add Basis functions
                            u(o++) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1];
                        }
                    }
                }
            }
        }
    } // end of faces

     // interiors
     for(int i=2; i<=p;i++)
     {
         for(int j=1; j<=p;j++)
         {
             for(int l=1; l<=p;l++)
             {
                 if((i+j+l)<=p)
                 {
                     // Define Barycentric Coordinates
                     La = bary_vector[0];
                     Lb = bary_vector[1];
                     Lc = bary_vector[2];
                     Ld = bary_vector[3];
                     
                     // compute polynomials
                     double x = Lb;
                     double y = La + Lb;
                     std::vector<double> int_Legendre_i;
                     poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                     
                     double alpha = 2.0*i;
                     x = Lc;
                     y = La + Lb + Lc;
                     std::vector<double> int_Jacobi_j;
                     poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                     
                     alpha = 2.0*(i+j);
                     x = Ld;
                     y = La + Lb + Lc + Ld;
                     std::vector<double> int_Jacobi_l;
                     poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                     
                     // Add Basis Functions
                     u(o++) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1];

                 }
             }
         }
     }// end of interiors
    
    Ti.Mult(u, shape);
}

void H1_TetrahedronElement_Fuentes::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   DenseMatrix du(dof, dim);
#endif

    int o = 0;

    double L1, L2, L3, L4;
    //compute barycentric coordinates as function of ip
    double L[3];
    ip.Get(L,3);
    std::vector<double> bary_vector{L[0], L[1], L[2], (1.0-L[0]-L[1]-L[2])};
    double La, Lb, Lc, Ld;
    
    // compute grad of barycentric coordinates
    std::vector<std::vector<double>> grad_bary_vectors{{1.0,0,0},{0,1.0,0},{0,0,1.0},{-1.0,-1.0,-1.0}};
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld;
    
   //Vertices
    for (int a=0; a<4; a++)
    {
        //compute vertex
        du(o,0) = grad_bary_vectors[a][0];
        du(o,1) = grad_bary_vectors[a][1];
        du(o,2) = grad_bary_vectors[a][2];
        o++;
    }
    
   //Edges
    for(int i=2; i<=p; i++)
    {
        for(int a=0; a<4; a++)
        {
            for(int b=0; b<4; b++)
            {
                if(a<b)
                {
                    // Define Barycentric coordinates
                    La = bary_vector[a];
                    Lb = bary_vector[b];
                    
                    // compute polynomials
                    std::vector<double> Legendre_i; // Use i-1 polynomial
                    double x = Lb;
                    double y = La + Lb;
                
                    poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                    
                    std::vector<double> R_i;
                    x = Lb;
                    y = La + Lb;
                    poly1d.CalcRLegendre(i, x, y, R_i);
                    
                    
                    // Add Basis Funcitons
                    du(o,0) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][0] + R_i[R_i.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0]);
                    
                    du(o,1) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][1] + R_i[R_i.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1]);
                    
                    du(o,2) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][2] + R_i[R_i.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2]);
                                    
                    o++;
                }
            }
        }
    } // end of edges
    
   //Faces
   for(int i=2; i<=p;i++)
   {
       for(int j=1; j<=p;j++)
       {
           for(int a=0; a<5;a++)
           {
               for(int b=0; b<5;b++)
               {
                   for(int c=0; c<5;c++)
                   {
                       if((a<b)&&(b<c)&&((i+j)<=p))
                       {
                           // Define Barycentric coordinates
                           La = bary_vector[a];
                           Lb = bary_vector[b];
                           Lc = bary_vector[c];
                           
                           // compute polynomials
                           double x = Lb;
                           double y = La + Lb;
                           std::vector<double> int_Legendre_i;
                           poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                           
                           double alpha = 2.0*i;
                           x = Lc;
                           y = La + Lb + Lc;
                           std::vector<double> int_Jacobi_j;
                           poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                           
                           // compute polynomials for grad(L_i)
                           std::vector<double> Legendre_i; // Use i-1 polynomial
                           x = Lb;
                           y = La + Lb;
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           std::vector<double> R_i; // Use i-1 polynomial
                           x = Lb;
                           y = La + Lb;
                           poly1d.CalcRLegendre(i, x, y, R_i);
                           
                           // compute polynimials for grad(L_j^2i)
                           std::vector<double> Jacobi_j; // Use i-1 polynomial
                           x = Lc;
                           y = La + Lb + Lc;
                           alpha = 2.0*i;
                           poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                           
                           std::vector<double> R_j; // Use i-1 polynomial
                           x = Lc;
                           y = La + Lb + Lc;
                           poly1d.CalcRJacobi(j, x, y, alpha, R_j);
                                                      
                           
                           // Add Basis Funcitons
                           du(o,0) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][0] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][0] +grad_bary_vectors[b][0]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][0] + R_j[R_j.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0] + grad_bary_vectors[c][0]));
                           
                           du(o,1) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][1] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][1] +grad_bary_vectors[b][1]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][1] + R_j[R_j.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1] + grad_bary_vectors[c][1]));
                           
                           du(o,2) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][2] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][2] +grad_bary_vectors[b][2]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][2] + R_j[R_j.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2] + grad_bary_vectors[c][2]));
                                                      
                           o++;
                       }
                   }
               }
           }
       }
   } // end of faces
    
    
  
  // interiors
    
  // Define Barycentric Coordinates
  La = bary_vector[0];
  Lb = bary_vector[1];
  Lc = bary_vector[2];
  Ld = bary_vector[3];
  
  for(int i=2; i<=p;i++)
  {
      for(int j=1; j<=p;j++)
      {
          for(int l=1; l<=p;l++)
          {
              if((i+j+l)<=p)
              {
                  
                  // compute polynomials
                  double x = Lb;
                  double y = La + Lb;
                  std::vector<double> int_Legendre_i;
                  poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                  
                  double alpha = 2.0*i;
                  x = Lc;
                  y = La + Lb + Lc;
                  std::vector<double> int_Jacobi_j;
                  poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                                    
                  alpha = 2.0*(i+j);
                  x = Ld;
                  y = La + Lb + Lc + Ld;
                  std::vector<double> int_Jacobi_l;
                  poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                  
                  // compute polynomials for grad(L_i)
                  std::vector<double> Legendre_i; // Use i-1 polynomial
                  x = Lb;
                  y = La + Lb;
                  poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                  
                  std::vector<double> R_i; // Use i-1 polynomial
                  x = Lb;
                  y = La + Lb;
                  poly1d.CalcRLegendre(i, x, y, R_i);
                  
                  // compute polynimials for grad(L_j^2i)
                  std::vector<double> Jacobi_j; // Use i-1 polynomial
                  x = Lc;
                  y = La + Lb + Lc;
                  alpha = 2.0*i;
                  poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                  
                  std::vector<double> R_j; // Use i-1 polynomial
                  x = Lc;
                  y = La + Lb + Lc;
                  poly1d.CalcRJacobi(j, x, y, alpha, R_j);
                                    
                  // compute polynimials for grad(L_l^2(i+j))
                  std::vector<double> Jacobi_l; // Use i-1 polynomial
                  x = Ld;
                  y = La + Lb + Lc + Ld;
                  alpha = 2.0*(i+j);
                  poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                  
                  std::vector<double> R_l; // Use i-1 polynomial
                  x = Ld;
                  y = La + Lb + Lc + Ld;
                  poly1d.CalcRJacobi(l, x, y, alpha, R_l);
                  
                  // Add Basis Funcitons
                  du(o,0) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][0] +
                             R_i[R_i.size()-2] * (grad_bary_vectors[0][0] +grad_bary_vectors[1][0]))
                             * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                  + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][0] + R_j[R_j.size()-2] * (grad_bary_vectors[0][0] + grad_bary_vectors[1][0] + grad_bary_vectors[2][0])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                  + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][0] + R_l[R_l.size()-2] * (grad_bary_vectors[0][0] + grad_bary_vectors[1][0] + grad_bary_vectors[2][0] + grad_bary_vectors[3][0])));
                  
                  du(o,1) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][1] +
                             R_i[R_i.size()-2] * (grad_bary_vectors[0][1] +grad_bary_vectors[1][1]))
                             * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                  + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][1] + R_j[R_j.size()-2] * (grad_bary_vectors[0][1] + grad_bary_vectors[1][1] + grad_bary_vectors[2][1])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                  + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][1] + R_l[R_l.size()-2] * (grad_bary_vectors[0][1] + grad_bary_vectors[1][1] + grad_bary_vectors[2][1] + grad_bary_vectors[3][1])));
                  
                  du(o,2) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][2] +
                             R_i[R_i.size()-2] * (grad_bary_vectors[0][2] +grad_bary_vectors[1][2]))
                             * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                  + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][2] + R_j[R_j.size()-2] * (grad_bary_vectors[0][2] + grad_bary_vectors[1][2] + grad_bary_vectors[2][2])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                  + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][2] + R_l[R_l.size()-2] * (grad_bary_vectors[0][2] + grad_bary_vectors[1][2] + grad_bary_vectors[2][2] + grad_bary_vectors[3][2])));
                  
                  o++;
                      
              }
          }
      }
  } // end of interiors

   Ti.Mult(du, dshape);
}

void H1_TetrahedronElement_Fuentes::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
    //Do Nothing
    mfem_error("H1_TetrahedronElement_Fuentes::CalcHessian not implemented");
}


H1_PentatopeElement::H1_PentatopeElement(const int p, const int type)
   : NodalFiniteElement(4, Geometry::PENTATOPE,
                        ((p + 1)*(p + 2)*(p + 3)*(p + 4))/24,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyClosed(type));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_t.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof,dim*(dim+1)/2 );
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
#endif

   // vertices
   Nodes.IntPoint(0).Set4(cp[0], cp[0], cp[0], cp[0]);
   Nodes.IntPoint(1).Set4(cp[p], cp[0], cp[0], cp[0]);
   Nodes.IntPoint(2).Set4(cp[0], cp[p], cp[0], cp[0]);
   Nodes.IntPoint(3).Set4(cp[0], cp[0], cp[p], cp[0]);
   Nodes.IntPoint(4).Set4(cp[0], cp[0], cp[0], cp[p]);

   // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
   int o = 5;
   for (int i = 1; i < p; i++)  // (0,1)
   {
      Nodes.IntPoint(o++).Set4(cp[i], cp[0], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,2)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,3)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (1,2)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,3)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,4)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (2,3)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[p-i], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (2,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[p-i], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (3,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[p-i], cp[i]);
   }

   // planars (see Mesh::GeneratePlanars in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,2)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[0], cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (2,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[p-i-j]/w, cp[i]/w, cp[j]/w);
      }

   // face(volumes)s (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,2,3)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[0]);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,2,1,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[j]/w, cp[i]/w, cp[0], cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,3,2,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[0], cp[j]/w, cp[i]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (1,2,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[p-i-j-k]/w, cp[i]/w, cp[j]/w, cp[k]/w);
         }

   // interior
   for (int l = 1; l < p; l++)
      for (int k = 1; k + l < p; k++)
         for (int j = 1; j + k + l < p; j++)
            for (int i = 1; i + j + k + l < p; i++)
            {
               double w = cp[i] + cp[j] + cp[k] + cp[l] + cp[p-i-j-k-l];
               Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[l]/w);
            }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, ip.t, shape_t);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);
       


      o = 0;
      for (int l = 0; l <= p; l++)
         for (int k = 0; k + l <= p; k++)
            for (int j = 0; j + k +l <= p; j++)
               for (int i = 0; i + j + k + l <= p; i++)
               {
                  T(o++, m) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
                
               }
   }


   Ti.Factor(T);
   // cout << "H1_PentatopeElement(" << p << ") : "; Ti.TestInversion();
}

void H1_PentatopeElement::CalcShape(const IntegrationPoint &ip,
                                    Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
   Vector u(Dof);
#endif

    
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);
    int local_dof = 0;

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
                local_dof = o;
            }

   Ti.Mult(u, shape);
    


    
}

void H1_PentatopeElement::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_t(p+1),
           shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_t(p+1),
          dshape_l(p + 1);
   DenseMatrix du(Dof, Dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, ip.t, shape_t, dshape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l, dshape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               int m = p - i - j - k - l;
               du(o,0) = ((dshape_x(i)* shape_l(m)) -
                          ( shape_x(i)*dshape_l(m)))*shape_y(j)*shape_z(k)*shape_t(l);
               du(o,1) = ((dshape_y(j)* shape_l(m)) -
                          ( shape_y(j)*dshape_l(m)))*shape_x(i)*shape_z(k)*shape_t(l);
               du(o,2) = ((dshape_z(k)* shape_l(m)) -
                          ( shape_z(k)*dshape_l(m)))*shape_x(i)*shape_y(j)*shape_t(l);
               du(o,3) = ((dshape_t(l)* shape_l(m)) -
                          ( shape_t(l)*dshape_l(m)))*shape_x(i)*shape_y(j)*shape_z(k);
               o++;
            }

   Ti.Mult(du, dshape);
}

void H1_PentatopeElement::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector   shape_x(p + 1),   shape_y(p + 1),   shape_z(p + 1),   shape_t(p+1),
            shape_l(p + 1);
   Vector  dshape_x(p + 1),  dshape_y(p + 1),  dshape_z(p + 1),  dshape_t(p+1),
           dshape_l(p + 1);
   Vector ddshape_x(p + 1), ddshape_y(p + 1), ddshape_z(p + 1), ddshape_t(p+1),
          ddshape_l(p + 1);
   DenseMatrix ddu(Dof, ((Dim+1)*Dim)/2);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x, ddshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y, ddshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z, ddshape_z);
   poly1d.CalcBasis(p, ip.t, shape_t, dshape_t, ddshape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l, dshape_l,
                    ddshape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               // u_xx, u_xy, u_xz, u_xt, u_yy, u_yz, u_yt, u_zz, u_zt, u_tt
               int m = p - i - j - k - l;
               ddu(o,0) = ((ddshape_x(i)*shape_l(m)) - 2.* (dshape_x(i)*dshape_l(m)) +
                           (shape_x(i)*ddshape_l(m))) * shape_y(j) * shape_z(k) * shape_t(l);
               ddu(o,1) = ((dshape_y(j)* ( (dshape_x(i)*shape_l(m)) - (shape_x(i)*dshape_l(
                                                                          m))) ) + (shape_y(j)* ((ddshape_l(m)*shape_x(i)) - (dshape_x(i) * dshape_l(
                                                                                       m)) ) ) )* shape_z(k) * shape_t(l);
               ddu(o,2) = ((dshape_z(k)* ( (dshape_x(i)*shape_l(m)) - (shape_x(i)*dshape_l(
                                                                          m))) ) + (shape_z(k)* ((ddshape_l(m)*shape_x(i)) - (dshape_x(i) * dshape_l(
                                                                                       m)) ) ) )* shape_y(j) * shape_t(l);
               ddu(o,3) = ((dshape_t(l)* ( (dshape_x(i)*shape_l(m)) - (shape_x(i)*dshape_l(
                                                                          m))) ) + (shape_t(l)* ((ddshape_l(m)*shape_x(i)) - (dshape_x(i) * dshape_l(
                                                                                       m)) ) ) )* shape_y(j) * shape_z(k);
               ddu(o,4) = ((ddshape_y(j)*shape_l(m)) - 2.* (dshape_y(j)*dshape_l(m)) +
                           (shape_y(j)*ddshape_l(m))) * shape_x(i) * shape_z(k) * shape_t(l);
               ddu(o,5) = ((dshape_z(k)* ( (dshape_y(j)*shape_l(m)) - (shape_y(j)*dshape_l(
                                                                          m))) ) + (shape_z(k)* ((ddshape_l(m)*shape_y(j)) - (dshape_y(j) * dshape_l(
                                                                                       m)) ) ) )* shape_x(i) * shape_t(l);
               ddu(o,6) = ((dshape_t(l)* ( (dshape_y(j)*shape_l(m)) - (shape_y(j)*dshape_l(
                                                                          m))) ) + (shape_t(l)* ((ddshape_l(m)*shape_y(j)) - (dshape_y(j) * dshape_l(
                                                                                       m)) ) ) )* shape_x(i) * shape_z(k);
               ddu(o,7) = ((ddshape_z(k)*shape_l(m)) - 2.* (dshape_z(k)*dshape_l(m)) +
                           (shape_z(k)*ddshape_l(m))) * shape_y(j) * shape_x(i) * shape_t(l);
               ddu(o,8) = ((dshape_t(l)* ( (dshape_z(k)*shape_l(m)) - (shape_z(k)*dshape_l(
                                                                          m))) ) + (shape_t(l)* ((ddshape_l(m)*shape_z(k)) - (dshape_z(k) * dshape_l(
                                                                                       m)) ) ) )* shape_x(i) * shape_y(j);
               ddu(o,9) = ((ddshape_t(l)*shape_l(m)) - 2.* (dshape_t(l)*dshape_l(m)) +
                           (shape_t(l)*ddshape_l(m))) * shape_y(j) * shape_x(i) * shape_z(k);
               o++;
            }
   Ti.Mult(ddu, ddshape);
}

H1_PentatopeElement_Barycentric::H1_PentatopeElement_Barycentric(const int p, const int type)
   : NodalFiniteElement(4, Geometry::PENTATOPE,
                        ((p + 1)*(p + 2)*(p + 3)*(p + 4))/24,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyClosed(type));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_t.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof,dim*(dim+1)/2 );
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
#endif

   // vertices
   Nodes.IntPoint(0).Set4(cp[0], cp[0], cp[0], cp[0]);
   Nodes.IntPoint(1).Set4(cp[p], cp[0], cp[0], cp[0]);
   Nodes.IntPoint(2).Set4(cp[0], cp[p], cp[0], cp[0]);
   Nodes.IntPoint(3).Set4(cp[0], cp[0], cp[p], cp[0]);
   Nodes.IntPoint(4).Set4(cp[0], cp[0], cp[0], cp[p]);

   // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
   int o = 5;
   for (int i = 1; i < p; i++)  // (0,1)
   {
      Nodes.IntPoint(o++).Set4(cp[i], cp[0], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,2)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,3)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (0,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (1,2)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[i], cp[0], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,3)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[0], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (1,4)
   {
      Nodes.IntPoint(o++).Set4(cp[p-i], cp[0], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (2,3)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[p-i], cp[i], cp[0]);
   }
   for (int i = 1; i < p; i++)  // (2,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[p-i], cp[0], cp[i]);
   }
   for (int i = 1; i < p; i++)  // (3,4)
   {
      Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[p-i], cp[i]);
   }

   // faces (see Mesh::GeneratePlanars in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,2)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[0], cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (2,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[p-i-j]/w, cp[i]/w, cp[j]/w);
      }

   // facets (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,2,3)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[0]);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,2,1,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[j]/w, cp[i]/w, cp[0], cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,3,2,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[0], cp[j]/w, cp[i]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (1,2,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[p-i-j-k]/w, cp[i]/w, cp[j]/w, cp[k]/w);
         }

   // interior bubbles
   for (int l = 1; l < p; l++)
      for (int k = 1; k + l < p; k++)
         for (int j = 1; j + k + l < p; j++)
            for (int i = 1; i + j + k + l < p; i++)
            {
               double w = cp[i] + cp[j] + cp[k] + cp[l] + cp[p-i-j-k-l];
               Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[l]/w);
            }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
       IntegrationPoint &ip = Nodes.IntPoint(m);
       
       //std::cout << "Dof = " << m << " " << ip.x << "," << ip.y << "," << ip.z << "," << ip.t << std::endl;
       
       std::vector<double> bary_vector{(1.0 - ip.x - ip.y - ip.z - ip.t), ip.x, ip.y, ip.z, ip.t};
       
       int o = 0;
       
       double La, Lb, Lc, Ld, Le;
       
       //Vertices
       for (int a=0; a<5; a++)
       {
           //compute barycentric function
           T(o++,m) = bary_vector[a];
       }
       
       //Edges
       for(int i=2; i<=p; i++)
       {
           for(int a=0; a<5; a++)
           {
               for(int b=0; b<5; b++)
               {
                   if(a<b)
                   {
                       La = bary_vector[a];
                       Lb = bary_vector[b];
                       
                       // compute polynomials
                       std::vector<double> int_Legendre_i;
                       double x = Lb;
                       double y = La + Lb;
                       poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                       // Add Basis Funcitons
                       T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1];
                       
                   }
               }
           }
       }// end of edges
       
       //Faces
       for(int i=2; i<=p;i++)
       {
           for(int j=1; j<=p;j++)
           {
               for(int a=0; a<5;a++)
               {
                   for(int b=0; b<5;b++)
                   {
                       for(int c=0; c<5;c++)
                       {
                           if((a<b)&&(b<c)&&((i+j)<=p))
                           {
                               La = bary_vector[a];
                               Lb = bary_vector[b];
                               Lc = bary_vector[c];
                               
                               // compute polynomials
                               double x = Lb;
                               double y = La + Lb;
                               std::vector<double> int_Legendre_i;
                               poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                               
                               double alpha = 2.0*i;
                               x = Lc;
                               y = La + Lb + Lc;
                               std::vector<double> int_Jacobi_j;
                               poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                               
                               // Add Basis functions
                               T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1];
                           }
                       }
                   }
               }
           }
       } 
       // end of faces
       
      //Facets
      int a;
      int b;
      int c;
      int d;

      for (int f=0; f<5; f++)
      {
          /*
          // Define each facet
          if (f==0)
          {
              a = 0;
              b = 1;
              c = 2;
              d = 3;
          }
          else if(f==1)
          {
              // Define Facet
              a = 0;
              b = 2;
              c = 1;
              d = 4;
          }
          // Define each facet
          else if (f==2)
          {
              // Define Facet
              a = 0;
              b = 1;
              c = 3;
              d = 4;
          }
          // Define each facet
          else if (f==3)
          {
              // Define Facet
              a = 0;
              b = 3;
              c = 2;
              d = 4;
          }
          // Define each facet
          else if (f==4)
          {
              // Define Facet
              a = 1;
              b = 2;
              c = 3;
              d = 4;
          }
          else
          {
              mfem_error("Invaild facet");
          }*/
          
          // Define each facet
          if (f==0)
          {
              a = 0;
              b = 1;
              c = 3;
              d = 2;
          }
          else if(f==1)
          {
              // Define Facet
              a = 0;
              b = 1;
              c = 4;
              d = 2;
          }
          // Define each facet
          else if (f==2)
          {
              // Define Facet
              a = 0;
              b = 1;
              c = 3;
              d = 4;
          }
          // Define each facet
          else if (f==3)
          {
              // Define Facet
              a = 0;
              b = 2;
              c = 3;
              d = 4;
          }
          // Define each facet
          else if (f==4)
          {
              // Define Facet
              a = 1;
              b = 2;
              c = 3;
              d = 4;
          }
          else
          {
              mfem_error("Invaild facet");
          }

          // Define Barycentric Coordinates
          La = bary_vector[a];
          Lb = bary_vector[b];
          Lc = bary_vector[c];
          Ld = bary_vector[d];

          for(int i=2; i<=p;i++)
          {
              for(int j=1; j<=p;j++)
              {
                  for(int l=1; l<=p;l++)
                  {
                      if((i+j+l)<=p)
                      {

                          // compute polynomials
                          double x = Lb;
                          double y = La + Lb;
                          std::vector<double> int_Legendre_i;
                          poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);

                          double alpha = 2.0*i;
                          x = Lc;
                          y = La + Lb + Lc;
                          std::vector<double> int_Jacobi_j;
                          poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);

                          alpha = 2.0*(i+j);
                          x = Ld;
                          y = La + Lb + Lc + Ld;
                          std::vector<double> int_Jacobi_l;
                          poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);

                          // Add Basis Functions
                          T(o++, m) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1];

                      }
                  }
              }
          }
      }// end of Facets

       // define lamda
       La = bary_vector[0];
       Lb = bary_vector[1];
       Lc = bary_vector[2];
       Ld = bary_vector[3];
       Le = bary_vector[4];
       
      //Interiors
      for(int i=2; i<=p;i++)
      {
          for(int j=1; j<=p;j++)
          {
              for(int l=1; l<=p;l++)
              {
                  for(int q=1; q<=p;q++)
                  {
                      if ((i+j+l+q)<=p)
                      {

                          // compute polynomials
                          double x = Lb;
                          double y = La + Lb;
                          std::vector<double> int_Legendre_i;
                          poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                          
                          double alpha = 2.0*i;
                          x = Lc;
                          y = La + Lb + Lc;
                          std::vector<double> int_Jacobi_j;
                          poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                          
                          alpha = 2.0*(i+j);
                          x = Ld;
                          y = La + Lb + Lc + Ld;
                          std::vector<double> int_Jacobi_l;
                          poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                          
                          alpha = 2.0*(i + j + l);
                          x = Le;
                          y = 1.0;
                          std::vector<double> int_Jacobi_q;
                          poly1d.CalcIntJacobi(q, x, y, alpha, int_Jacobi_q);
                          
                          //Add Basis Functions
                          T(o++, m) = (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_q[int_Jacobi_q.size()-1]);
                          //T(o++, m) = (2.0*La -1.0)*(2.0*Lb -1.0)*(2.0*Lc -1.0)*(2.0*Ld -1.0)*(2.0*Le-1.0);
                          //T(o++, m) = La*Lb*Lc*Ld*Le;

                      }
                  }
              }
          }
      } //end of Interiors
       
//       
//       ofstream logfile;
//       logfile.open ("logfile.txt");
//       for (int i =0; i<T.Size(); i++) {
//           for (int j = 0; j<T.Size(); j++) {
//               logfile << "T value @ " << i << "," << j << "= " << T(i,j) << std::endl;
//           }
//       }
//       logfile.close();
//       
//       std::ofstream A_file("T_matrix.txt");
//
//       T.PrintMatlab(A_file);


       //mfem_error("Stop for debug");


   }

   Ti.Factor(T);
    
    //std::ofstream AA_file("Ti_matrix.txt");

    //Ti.PrintMatlab(AA_file);
   cout << "H1_PentatopeElement(" << p << ") : "; Ti.TestInversion();
}


void H1_PentatopeElement_Barycentric::CalcShape(const IntegrationPoint &ip,
                                                Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
   Vector u(Dof);
#endif
    
    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);
    
    //ofstream logfile;
    //logfile.open ("logfile_shape.txt", "a");
    //ofstream logfile("logfile_shape.txt", std::ios_base::app);
    //logfile << "pts = " <<  ip.x  << ", " << ip.y  << ", " << ip.z  << ", " << ip.t << std::endl;
    
    std::vector<double> bary_vector{(1.0 - L[0]-L[1]-L[2]-L[3]), L[0], L[1], L[2], L[3]};
    
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;
    
    int local_dof = 0;
    
   //Vertices
    for (int a=0; a<5; a++)
    {
        //compute barycentric function
        u(o++) = bary_vector[a];
    }
    
   //Edges
    for(int i=2; i<=p; i++)
    {
        for(int a=0; a<5; a++)
        {
            for(int b=0; b<5; b++)
            {
                if(a<b)
                {
                    La = bary_vector[a];
                    Lb = bary_vector[b];
                    
                    // compute polynomials
                    std::vector<double> int_Legendre_i;
                    double x = Lb;
                    double y = La + Lb;
                    poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                                        
                    // Add Basis Funcitons
                    u(o++) = int_Legendre_i[int_Legendre_i.size()-1];
                }
            }
        }
    } // end of edges

    
   //Faces
   for(int i=2; i<=p;i++)
   {
       for(int j=1; j<=p;j++)
       {
           for(int a=0; a<5;a++)
           {
               for(int b=0; b<5;b++)
               {
                   for(int c=0; c<5;c++)
                   {
                       if((a<b)&&(b<c)&&((i+j)<=p))
                       {
                           La = bary_vector[a];
                           Lb = bary_vector[b];
                           Lc = bary_vector[c];
                           
                           // compute polynomials
                           double x = Lb;
                           double y = La + Lb;
                           std::vector<double> int_Legendre_i;
                           poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                           
                           double alpha = 2.0*i;
                           x = Lc;
                           y = La + Lb + Lc;
                           std::vector<double> int_Jacobi_j;
                           poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                           
                           // Add Basis functions
                           u(o++) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1];
                       }
                   }
               }
           }
       }
   } // end of faces
   
    
    //Facets
    int a;
    int b;
    int c;
    int d;
      
    for (int f=0; f<5; f++)
    {
        /*// Define each facet
        if (f==0)
        {
            a = 0;
            b = 1;
            c = 2;
            d = 3;
        }
        else if(f==1)
        {
            // Define Facet
            a = 0;
            b = 2;
            c = 1;
            d = 4;
        }
        // Define each facet
        else if (f==2)
        {
            // Define Facet
            a = 0;
            b = 1;
            c = 3;
            d = 4;
        }
        // Define each facet
        else if (f==3)
        {
            // Define Facet
            a = 0;
            b = 3;
            c = 2;
            d = 4;
        }
        // Define each facet
        else if (f==4)
        {
            // Define Facet
            a = 1;
            b = 2;
            c = 3;
            d = 4;
        }
        else
        {
            mfem_error("Invaild facet");
        }*/
        
        // Define each facet
        if (f==0)
        {
            a = 0;
            b = 1;
            c = 3;
            d = 2;
        }
        else if(f==1)
        {
            // Define Facet
            a = 0;
            b = 1;
            c = 4;
            d = 2;
        }
        // Define each facet
        else if (f==2)
        {
            // Define Facet
            a = 0;
            b = 1;
            c = 3;
            d = 4;
        }
        // Define each facet
        else if (f==3)
        {
            // Define Facet
            a = 0;
            b = 2;
            c = 3;
            d = 4;
        }
        // Define each facet
        else if (f==4)
        {
            // Define Facet
            a = 1;
            b = 2;
            c = 3;
            d = 4;
        }
        else
        {
            mfem_error("Invaild facet");
        }
        
        // Define Barycentric Coordinates
        La = bary_vector[a];
        Lb = bary_vector[b];
        Lc = bary_vector[c];
        Ld = bary_vector[d];
        
        for(int i=2; i<=p;i++)
        {
            for(int j=1; j<=p;j++)
            {
                for(int l=1; l<=p;l++)
                {
                    if((i+j+l)<=p)
                    {
                        
                        // compute polynomials
                        double x = Lb;
                        double y = La + Lb;
                        std::vector<double> int_Legendre_i;
                        poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                        
                        double alpha = 2.0*i;
                        x = Lc;
                        y = La + Lb + Lc;
                        std::vector<double> int_Jacobi_j;
                        poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                        
                        alpha = 2.0*(i+j);
                        x = Ld;
                        y = La + Lb + Lc + Ld;
                        std::vector<double> int_Jacobi_l;
                        poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                        
                        // Add Basis Functions
                        u(o++) = int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1];
                        local_dof = o;


                    }
                }
            }
        }
    }// end of Facets

    La = bary_vector[0];
    Lb = bary_vector[1];
    Lc = bary_vector[2];
    Ld = bary_vector[3];
    Le = bary_vector[4];
           
    //Interiors
   for(int i=2; i<=p;i++)
   {
       for(int j=1; j<=p;j++)
       {
           for(int l=1; l<=p;l++)
           {
               for(int q=1; q<=p;q++)
               {
                   if ((i+j+l+q)<=p)
                   {
                       // define lamda
                       
//                       std::cout << "La= " << La << std::endl;
//                       std::cout << "Lb= " << Lb << std::endl;
//                       std::cout << "Lc= " << Lc << std::endl;
//                       std::cout << "Ld= " << Ld << std::endl;
//                       std::cout << "Le= " << Le << std::endl;
                       
                        //compute polynomials
                       double x = Lb;
                       double y = La + Lb;
                       std::vector<double> int_Legendre_i;
                       poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                       
                       //std::cout << "o = " << o << std::endl;
                       //std::cout << "x(x-y) :" << x*(x-y) << std::endl;
                       
                       double alpha = 2.0*i;
                       x = Lc;
                       y = La + Lb + Lc;
                       std::vector<double> int_Jacobi_j;
                       poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                       
                       alpha = 2.0*(i+j);
                       x = Ld;
                       y = La + Lb + Lc + Ld;
                       std::vector<double> int_Jacobi_l;
                       poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                       
                       alpha = 2.0*(i + j + l);
                       x = Le;
                       y = 1.0;
                       std::vector<double> int_Jacobi_q;
                       poly1d.CalcIntJacobi(q, x, y, alpha, int_Jacobi_q);

                        //Add Basis Functions
//
//                       std::cout << "Bubble Value :" << int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_q[int_Jacobi_q.size()-1] << std::endl;
//                       std::cout << "int_legendre :" << int_Legendre_i[int_Legendre_i.size()-1] << std::endl;
//                       std::cout << "int_Jacobi_j :" << int_Jacobi_j[int_Jacobi_j.size()-1] << std::endl;
//                       std::cout << "int_Jacobi_l :" << int_Jacobi_l[int_Jacobi_l.size()-1] << std::endl;
//                       std::cout << "int_Jacobi_q :" << int_Jacobi_q[int_Jacobi_q.size()-1] << std::endl;
//
//                       std::cout << "Product Bary :" << La*Lb*Lc*Ld*Le << std::endl;
                       
                       u(o++) = (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_q[int_Jacobi_q.size()-1]);
                       //u(o++) = (2.0*La -1.0)*(2.0*Lb -1.0)*(2.0*Lc -1.0)*(2.0*Ld -1.0)*(2.0*Le-1.0);
                       //u(o++) = La*Lb*Lc*Ld*Le;

                   }
               }
           }
       }
   } //end of Interiors */
   //std::cout << "end of interiors" << std::endl;
    

   Ti.Mult(u, shape);
    
//    for (int i = 0; i<local_dof; i++) {
//        logfile << "shape " << i << "= " << shape(i) << std::endl;
//    }
//    logfile.close();

}

void H1_PentatopeElement_Barycentric::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_t(p+1),
           shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_t(p+1),
          dshape_l(p + 1);
   DenseMatrix du(Dof, Dim);
#endif
    int o = 0;

    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);
    std::vector<double> bary_vector{(1.0-L[0]-L[1]-L[2]-L[3]), L[0], L[1], L[2], L[3]};
    double La, Lb, Lc, Ld, Le;
    
    // compute grad of barycentric coordinates
    std::vector<std::vector<double>> grad_bary_vectors{{-1.0,-1.0,-1.0,-1.0},{1.0,0,0,0},{0,1.0,0,0},{0,0,1.0,0},{0,0,0,1.0}};
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    
   //Vertices
    for (int a=0; a<5; a++)
    {
        //compute vertex
        du(o,0) = grad_bary_vectors[a][0];
        du(o,1) = grad_bary_vectors[a][1];
        du(o,2) = grad_bary_vectors[a][2];
        du(o,3) = grad_bary_vectors[a][3];
        o++;
    }
    
   //Edges
    for(int i=2; i<=p; i++)
    {
        for(int a=0; a<5; a++)
        {
            for(int b=0; b<5; b++)
            {
                if(a<b)
                {
                    // Define Barycentric coordinates
                    La = bary_vector[a];
                    Lb = bary_vector[b];
                    
                    // compute polynomials
                    std::vector<double> Legendre_i; // Use i-1 polynomial
                    double x = Lb;
                    double y = La + Lb;
                
                    poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                    
                    std::vector<double> R_i;
                    x = Lb;
                    y = La + Lb;
                    poly1d.CalcRLegendre(i, x, y, R_i);
                    
                    
                    // Add Basis Funcitons
                    du(o,0) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][0] + R_i[R_i.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0]);
                    
                    du(o,1) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][1] + R_i[R_i.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1]);
                    
                    du(o,2) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][2] + R_i[R_i.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2]);
                    
                    du(o,3) = Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][3] + R_i[R_i.size()-2] * (grad_bary_vectors[a][3] + grad_bary_vectors[b][3]);
                
                    o++;
                }
            }
        }
    } // end of edges
    
   //Faces
   for(int i=2; i<=p;i++)
   {
       for(int j=1; j<=p;j++)
       {
           for(int a=0; a<5;a++)
           {
               for(int b=0; b<5;b++)
               {
                   for(int c=0; c<5;c++)
                   {
                       if((a<b)&&(b<c)&&((i+j)<=p))
                       {
                           // Define Barycentric coordinates
                           La = bary_vector[a];
                           Lb = bary_vector[b];
                           Lc = bary_vector[c];
                           
                           // compute polynomials
                           double x = Lb;
                           double y = La + Lb;
                           std::vector<double> int_Legendre_i;
                           poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                           
                           double alpha = 2.0*i;
                           x = Lc;
                           y = La + Lb + Lc;
                           std::vector<double> int_Jacobi_j;
                           poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                           
                           // compute polynomials for grad(L_i)
                           std::vector<double> Legendre_i; // Use i-1 polynomial
                           x = Lb;
                           y = La + Lb;
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           std::vector<double> R_i; // Use i-1 polynomial
                           x = Lb;
                           y = La + Lb;
                           poly1d.CalcRLegendre(i, x, y, R_i);
                           
                           // compute polynimials for grad(L_j^2i)
                           std::vector<double> Jacobi_j; // Use i-1 polynomial
                           x = Lc;
                           y = La + Lb + Lc;
                           alpha = 2.0*i;
                           poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                           
                           std::vector<double> R_j; // Use i-1 polynomial
                           x = Lc;
                           y = La + Lb + Lc;
                           poly1d.CalcRJacobi(j, x, y, alpha, R_j);

                           // Add Basis Funcitons
                           du(o,0) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][0] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][0] +grad_bary_vectors[b][0]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][0] + R_j[R_j.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0] + grad_bary_vectors[c][0]));
                           
                           du(o,1) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][1] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][1] +grad_bary_vectors[b][1]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][1] + R_j[R_j.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1] + grad_bary_vectors[c][1]));
                           
                           du(o,2) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][2] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][2] +grad_bary_vectors[b][2]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][2] + R_j[R_j.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2] + grad_bary_vectors[c][2]));
                           
                           du(o,3) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][3] +
                                      R_i[R_i.size()-2] * (grad_bary_vectors[a][3] +grad_bary_vectors[b][3]))
                                      * int_Jacobi_j[int_Jacobi_j.size()-1]
                           + int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][3] + R_j[R_j.size()-2] * (grad_bary_vectors[a][3] + grad_bary_vectors[b][3] + grad_bary_vectors[c][3]));
                           
                           o++;
                       }
                   }
               }
           }
       }
   } // end of faces
    
    
  //Facets
  int a;
  int b;
  int c;
  int d;
    
  for (int f=0; f<5; f++)
  {
      /*// Define each facet
      if (f==0)
      {
          a = 0;
          b = 1;
          c = 2;
          d = 3;
      }
      else if(f==1)
      {
          // Define Facet
          a = 0;
          b = 2;
          c = 1;
          d = 4;
      }
      // Define each facet
      else if (f==2)
      {
          // Define Facet
          a = 0;
          b = 1;
          c = 3;
          d = 4;
      }
      // Define each facet
      else if (f==3)
      {
          // Define Facet
          a = 0;
          b = 3;
          c = 2;
          d = 4;
      }
      // Define each facet
      else if (f==4)
      {
          // Define Facet
          a = 1;
          b = 2;
          c = 3;
          d = 4;
      }
      else
      {
          mfem_error("Invaild facet");
      }*/
      
      // Define each facet
      if (f==0)
      {
          a = 0;
          b = 1;
          c = 3;
          d = 2;
      }
      else if(f==1)
      {
          // Define Facet
          a = 0;
          b = 1;
          c = 4;
          d = 2;
      }
      // Define each facet
      else if (f==2)
      {
          // Define Facet
          a = 0;
          b = 1;
          c = 3;
          d = 4;
      }
      // Define each facet
      else if (f==3)
      {
          // Define Facet
          a = 0;
          b = 2;
          c = 3;
          d = 4;
      }
      // Define each facet
      else if (f==4)
      {
          // Define Facet
          a = 1;
          b = 2;
          c = 3;
          d = 4;
      }
      else
      {
          mfem_error("Invaild facet");
      }
      
      // Define Barycentric Coordinates
      La = bary_vector[a];
      Lb = bary_vector[b];
      Lc = bary_vector[c];
      Ld = bary_vector[d];
      
      for(int i=2; i<=p;i++)
      {
          for(int j=1; j<=p;j++)
          {
              for(int l=1; l<=p;l++)
              {
                  if((i+j+l)<=p)
                  {
                      
                      // compute polynomials
                      double x = Lb;
                      double y = La + Lb;
                      std::vector<double> int_Legendre_i;
                      poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                      
                      double alpha = 2.0*i;
                      x = Lc;
                      y = La + Lb + Lc;
                      std::vector<double> int_Jacobi_j;
                      poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                      
                      alpha = 2.0*(i+j);
                      x = Ld;
                      y = La + Lb + Lc + Ld;
                      std::vector<double> int_Jacobi_l;
                      poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                      
                      // compute polynomials for grad(L_i)
                      std::vector<double> Legendre_i; // Use i-1 polynomial
                      x = Lb;
                      y = La + Lb;
                      poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                      
                      std::vector<double> R_i; // Use i-1 polynomial
                      x = Lb;
                      y = La + Lb;
                      poly1d.CalcRLegendre(i, x, y, R_i);
                      
                      // compute polynimials for grad(L_j^2i)
                      std::vector<double> Jacobi_j; // Use i-1 polynomial
                      x = Lc;
                      y = La + Lb + Lc;
                      alpha = 2.0*i;
                      poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                      
                      std::vector<double> R_j; // Use i-1 polynomial
                      x = Lc;
                      y = La + Lb + Lc;
                      poly1d.CalcRJacobi(j, x, y, alpha, R_j);
                      
                      
                      // compute polynimials for grad(L_l^2(i+j))
                      std::vector<double> Jacobi_l; // Use i-1 polynomial
                      x = Ld;
                      y = La + Lb + Lc + Ld;
                      alpha = 2.0*(i+j);
                      poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                      
                      std::vector<double> R_l; // Use i-1 polynomial
                      x = Ld;
                      y = La + Lb + Lc + Ld;
                      poly1d.CalcRJacobi(l, x, y, alpha, R_l);
                      
                      // Add Basis Funcitons
                      du(o,0) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][0] +
                                 R_i[R_i.size()-2] * (grad_bary_vectors[a][0] +grad_bary_vectors[b][0]))
                                 * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                      + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][0] + R_j[R_j.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0] + grad_bary_vectors[c][0])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                      + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[d][0] + R_l[R_l.size()-2] * (grad_bary_vectors[a][0] + grad_bary_vectors[b][0] + grad_bary_vectors[c][0] + grad_bary_vectors[d][0])));
                      
                      du(o,1) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][1] +
                                 R_i[R_i.size()-2] * (grad_bary_vectors[a][1] +grad_bary_vectors[b][1]))
                                 * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                      + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][1] + R_j[R_j.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1] + grad_bary_vectors[c][1])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                      + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[d][1] + R_l[R_l.size()-2] * (grad_bary_vectors[a][1] + grad_bary_vectors[b][1] + grad_bary_vectors[c][1] + grad_bary_vectors[d][1])));
                      
                      du(o,2) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][2] +
                                 R_i[R_i.size()-2] * (grad_bary_vectors[a][2] +grad_bary_vectors[b][2]))
                                 * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                      + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][2] + R_j[R_j.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2] + grad_bary_vectors[c][2])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                      + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[d][2] + R_l[R_l.size()-2] * (grad_bary_vectors[a][2] + grad_bary_vectors[b][2] + grad_bary_vectors[c][2] + grad_bary_vectors[d][2])));
                      
                      du(o,3) = (Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[b][3] +
                                 R_i[R_i.size()-2] * (grad_bary_vectors[a][3] +grad_bary_vectors[b][3]))
                                 * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1]
                      + (int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[c][3] + R_j[R_j.size()-2] * (grad_bary_vectors[a][3] + grad_bary_vectors[b][3] + grad_bary_vectors[c][3])) * int_Jacobi_l[int_Jacobi_l.size()-1])
                      + (int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[d][3] + R_l[R_l.size()-2] * (grad_bary_vectors[a][3] + grad_bary_vectors[b][3] + grad_bary_vectors[c][3] + grad_bary_vectors[d][3])));
                      o++;
                          
                  }
              }
          }
      }
  }// end of Facets

           
    //Interiors
   for(int i=2; i<=p;i++)
   {
       for(int j=1; j<=p;j++)
       {
           for(int l=1; l<=p;l++)
           {
               for(int m=1; m<=p;m++)
               {
                   if ((i+j+l+m)<=p)
                   {
                       // define lamda
                       La = bary_vector[0];
                       Lb = bary_vector[1];
                       Lc = bary_vector[2];
                       Ld = bary_vector[3];
                       Le = bary_vector[4];
                       
                       // compute polynomials
                       double x = Lb;
                       double y = La + Lb;
                       std::vector<double> int_Legendre_i;
                       poly1d.CalcIntLegendre(i, x, y, int_Legendre_i);
                       
                       double alpha = 2.0*i;
                       x = Lc;
                       y = La + Lb + Lc;
                       std::vector<double> int_Jacobi_j;
                       poly1d.CalcIntJacobi(j, x, y, alpha, int_Jacobi_j);
                       
                       alpha = 2.0*(i+j);
                       x = Ld;
                       y = La + Lb + Lc + Ld;
                       std::vector<double> int_Jacobi_l;
                       poly1d.CalcIntJacobi(l, x, y, alpha, int_Jacobi_l);
                       
                       alpha = 2.0*(i + j + l);
                       x = Le;
                       y = 1.0;
                       std::vector<double> int_Jacobi_m;
                       poly1d.CalcIntJacobi(m, x, y, alpha, int_Jacobi_m);
                       
                       // compute polynomials for grad(L_i)
                       std::vector<double> Legendre_i; // Use i-1 polynomial
                       x = Lb;
                       y = La + Lb;
                       poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                       
                       std::vector<double> R_i; // Use i-1 polynomial
                       x = Lb;
                       y = La + Lb;
                       poly1d.CalcRLegendre(i, x, y, R_i);
                       
                       // compute polynimials for grad(L_j^2i)
                       std::vector<double> Jacobi_j; // Use i-1 polynomial
                       x = Lc;
                       y = La + Lb + Lc;
                       alpha = 2.0*i;
                       poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                       
                       std::vector<double> R_j; // Use i-1 polynomial
                       x = Lc;
                       y = La + Lb + Lc;
                       poly1d.CalcRJacobi(j, x, y, alpha, R_j);
                       
                       // compute polynimials for grad(L_l^2(i+j))
                       std::vector<double> Jacobi_l; // Use i-1 polynomial
                       x = Ld;
                       y = La + Lb + Lc + Ld;
                       alpha = 2.0*(i+j);
                       poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                       
                       std::vector<double> R_l; // Use i-1 polynomial
                       x = Ld;
                       y = La + Lb + Lc + Ld;
                       poly1d.CalcRJacobi(l, x, y, alpha, R_l);
                       
                       // compute polynimials for grad(L_m^2(i+j+l))
                       std::vector<double> Jacobi_m; // Use i-1 polynomial
                       x = Le;
                       y = 1.0;
                       alpha = 2.0*(i+j+l);
                       poly1d.CalcJacobi(m, x, y, alpha, Jacobi_m);
                       
//                       std::vector<double> R_m; // Use i-1 polynomial
//                       x = Le;
//                       y = 1;
//                       poly1d.CalcRJacobi(m, x, y, alpha, R_m);
                       
                       
                       // Add Basis Funcitons
                       //std::cout << "o = " << o << std::endl;
                       
                       double dA = ((Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][0] +
                                  R_i[R_i.size()-2] * (grad_bary_vectors[0][0] +grad_bary_vectors[1][0]))
                                  * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       +  int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][0] + R_j[R_j.size()-2] * (grad_bary_vectors[0][0] + grad_bary_vectors[1][0] + grad_bary_vectors[2][0])) * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][0] + R_l[R_l.size()-2] * (grad_bary_vectors[0][0] + grad_bary_vectors[1][0] + grad_bary_vectors[2][0] + grad_bary_vectors[3][0])) *  int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_m.size()-1] * (Jacobi_m[Jacobi_m.size()-2] * grad_bary_vectors[4][0]));
                       du(o, 0) = dA;
                       //std::cout << "dA = " << dA << std::endl;
                       
                       double dB = ((Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][1] +
                                  R_i[R_i.size()-2] * (grad_bary_vectors[0][1] +grad_bary_vectors[1][1]))
                                  * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       +  int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][1] + R_j[R_j.size()-2] * (grad_bary_vectors[0][1] + grad_bary_vectors[1][1] + grad_bary_vectors[2][1])) * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][1] + R_l[R_l.size()-2] * (grad_bary_vectors[0][1] + grad_bary_vectors[1][1] + grad_bary_vectors[2][1] + grad_bary_vectors[3][1])) *  int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_m.size()-1] * (Jacobi_m[Jacobi_m.size()-2] * grad_bary_vectors[4][1]));
                       du(o, 1) = dB;
                       //std::cout << "dB = " << dB << std::endl;

                       
                       double dC = ((Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][2] +
                                  R_i[R_i.size()-2] * (grad_bary_vectors[0][2] +grad_bary_vectors[1][2]))
                                  * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       +  int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][2] + R_j[R_j.size()-2] * (grad_bary_vectors[0][2] + grad_bary_vectors[1][2] + grad_bary_vectors[2][2])) * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][2] + R_l[R_l.size()-2] * (grad_bary_vectors[0][2] + grad_bary_vectors[1][2] + grad_bary_vectors[2][2] + grad_bary_vectors[3][2])) *  int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_m.size()-1] * (Jacobi_m[Jacobi_m.size()-2] * grad_bary_vectors[4][2]));
                       du(o, 2) = dC;
                       //std::cout << "dC = " << dC << std::endl;

                       
                       double dD = ((Legendre_i[Legendre_i.size()-2] * grad_bary_vectors[1][3] +
                                  R_i[R_i.size()-2] * (grad_bary_vectors[0][3] +grad_bary_vectors[1][3]))
                                  * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       +  int_Legendre_i[int_Legendre_i.size()-1] * (Jacobi_j[Jacobi_j.size()-2] * grad_bary_vectors[2][3] + R_j[R_j.size()-2] * (grad_bary_vectors[0][3] + grad_bary_vectors[1][3] + grad_bary_vectors[2][3])) * int_Jacobi_l[int_Jacobi_l.size()-1] * int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * (Jacobi_l[Jacobi_l.size()-2] * grad_bary_vectors[3][3] + R_l[R_l.size()-2] * (grad_bary_vectors[0][3] + grad_bary_vectors[1][3] + grad_bary_vectors[2][3] + grad_bary_vectors[3][3])) *  int_Jacobi_m[int_Jacobi_m.size()-1]
                       + int_Legendre_i[int_Legendre_i.size()-1] * int_Jacobi_j[int_Jacobi_j.size()-1] * int_Jacobi_l[int_Jacobi_m.size()-1] * (Jacobi_m[Jacobi_m.size()-2] * grad_bary_vectors[4][3]));
                       du(o, 3) = dD;
                       //std::cout << "dD = " << dD << std::endl;

                       /*
                       double testDA = (grad_bary_vectors[0][0]*Lb*Lc*Ld*Le + La*grad_bary_vectors[1][0]*Lc*Ld*Le
                                         + La*Lb*grad_bary_vectors[2][0]*Ld*Le + La*Lb*Lc*grad_bary_vectors[3][0]*Le
                                         + La*Lb*Lc*Ld*grad_bary_vectors[4][0]);
                       //std::cout << "testDA= " << testDA << std::endl;
                       du(o, 0) = testDA;

                       
                       double testDB = (grad_bary_vectors[0][1]*Lb*Lc*Ld*Le + La*grad_bary_vectors[1][1]*Lc*Ld*Le
                                         + La*Lb*grad_bary_vectors[2][1]*Ld*Le + La*Lb*Lc*grad_bary_vectors[3][1]*Le
                                         + La*Lb*Lc*Ld*grad_bary_vectors[4][1]);
                       //std::cout << "testDB= " << testDB << std::endl;
                       du(o, 1) = testDB;


                       
                       double testDC = (grad_bary_vectors[0][2]*Lb*Lc*Ld*Le + La*grad_bary_vectors[1][2]*Lc*Ld*Le
                                         + La*Lb*grad_bary_vectors[2][2]*Ld*Le + La*Lb*Lc*grad_bary_vectors[3][2]*Le
                                         + La*Lb*Lc*Ld*grad_bary_vectors[4][2]);
                       
                       //std::cout << "testDC= " << testDC << std::endl;
                       du(o, 2) = testDC;


                       
                       double testDD = (grad_bary_vectors[0][3]*Lb*Lc*Ld*Le + La*grad_bary_vectors[1][3]*Lc*Ld*Le
                                         + La*Lb*grad_bary_vectors[2][3]*Ld*Le + La*Lb*Lc*grad_bary_vectors[3][3]*Le
                                         + La*Lb*Lc*Ld*grad_bary_vectors[4][3]);
                       //std::cout << "testDD= " << testDD << std::endl;
                       du(o, 3) = testDD;


                       
                       o++;*/
                       
                   }
               }
           }
       }
   } //end of Interiors
    
   Ti.Mult(du, dshape);
}

void H1_PentatopeElement_Barycentric::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
    //Do Nothing
    mfem_error("H1_PentatopeElement_Barycentric::CalcHessian not implemented");
}


const double HSkwGrad_PentatopeElement::nk[20] =
{ 0,0,0,-1,  0,0,-1,0,  0,-1,0,0,  -1,0,0,0,  1,1,1,1};
// { .5,.5,.5, -.5,0,0, 0,-.5,0, 0,0,-.5}; // n_F |F|

const double HSkwGrad_PentatopeElement::c = 1./5.;

HSkwGrad_PentatopeElement::HSkwGrad_PentatopeElement(const int p)
   : VectorFiniteElement(4, Geometry::PENTATOPE, (p + 1)*(p + 2)*(p + 3)*(p + 5)/6,
                         p + 1, H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   const double *iop = (p > 0) ? poly1d.OpenPoints(p - 1) : NULL;
   const double *bop = poly1d.OpenPoints(p);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim);
   divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
#endif
   //std::cout << "p_V = " << p << std::endl;
   int o = 0;
   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp,
   //        the constructor of H1_PentatopeElement)
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,2,3)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, bop[j]/w, bop[k]/w, 0.);
            dof2nk[o++] = 0;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,2,1,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[j]/w, bop[i]/w, 0., bop[k]/w);
            dof2nk[o++] = 1;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, 0., bop[j]/w, bop[k]/w);
            dof2nk[o++] = 2;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,3,2,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(0., bop[j]/w, bop[i]/w, bop[k]/w);
            dof2nk[o++] = 3;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (1,2,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[p-i-j-k]/w, bop[i]/w, bop[j]/w, bop[k]/w);
            dof2nk[o++] = 4;
         }

   // interior
   for (int l = 0; l < p; l++)
      for (int k = 0; k + l < p; k++)
         for (int j = 0; j + k + l < p; j++)
            for (int i = 0; i + j + k + l < p; i++)
            {
                double w = iop[i] + iop[j] + iop[k] + iop[l] + iop[p-1-i-j-k-l];
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 1;
                dof2nk[o++] = 3;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 2;
                dof2nk[o++] = 2;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 3;
                dof2nk[o++] = 1;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 4;
                dof2nk[o++] = 0;
            }

   DenseMatrix T(dof);
   DenseMatrix B(dof, dim);
   

   mfem::out << "RT_PentatopeElement(" << p << ") : "; Ti.TestInversion();
}

void HSkwGrad_PentatopeElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   const int p = order - 1;
   //std::cout << "p_CalcShape = " << p << std::endl;
    int reset_o;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
   DenseMatrix u(Dof, Dim);
#endif
    
    /*//compute barycentric coordinates as function of ip

    std::vector<double> bary_vector{(1.0 -ip.x-ip.y-ip.z-ip.t), ip.x, ip.y, ip.z, ip.t};
    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{-1,-1,-1,-1};
    std::vector<double> gradL2{1,0,0,0};
    std::vector<double> gradL3{0,1,0,0};
    std::vector<double> gradL4{0,0,1,0};
    std::vector<double> gradL5{0,0,0,1};
    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;

   //Edges
    for(int i=0; i<=p-1; i++)
    {
        for(int a=0; a<5; a++)
        {
            for(int b=0; b<5; b++)
            {
                if(a<b)
                {
                    // set Lamda
                    La = bary_vector[a];
                    Lb = bary_vector[b];
                    // set grad Lamda
                    std::vector<double> grad_La = gradbary_vector[a];
                    std::vector<double> grad_Lb = gradbary_vector[b];
                    
                    // compute edge function
                    double x = Lb;
                    double y = (La+Lb);
                    std::vector<double> Legendre_i;
                    poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                    
                    for(int index=0; index<4; index++)
                    {
                        double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                        u(o++) = Legendre_i[Legendre_i.size()-1]  * scalar_val;
                    }
                    
                 }
               }
            }
        } // end of edges
    
   //Faces
   for(int i=0; i<=p-1;i++)
   {
       for(int j=1; j<=p-1;j++)
       {
           for(int a=0; a<5;a++)
           {
               for(int b=0; b<5;b++)
               {
                   for(int c=0; c<5;c++)
                   {
                       if((a<b)&&(b<c)&&((i+j)<=p-1))
                       {
                           
                           // Compute Family I:    <a,b,c>
                           
                           // set Lamda
                           La = bary_vector[a];
                           Lb = bary_vector[b];
                           Lc = bary_vector[c];
                           
                           // set grad Lamda
                           std::vector<double> grad_La = gradbary_vector[a];
                           std::vector<double> grad_Lb = gradbary_vector[b];
                           std::vector<double> grad_Lc = gradbary_vector[c];

                           // Compute Polynomials
                           double x = Lb;
                           double y = La + Lb;
                           std::vector<double> Legendre_i;
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           double alpha = (2*i+1);
                           std::vector<double> Int_Jacobi_j;
                           x = Lc;
                           y = La + Lb + Lc;
                           poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                           
                           for(int index=0; index<4; index++)
                           {
                               double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                               u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1]  * scalar_val;
                           }
                           
                           // Compute Family II:    <b,c,a>
                           
                           // set Lamda
                           La = bary_vector[b];
                           Lb = bary_vector[c];
                           Lc = bary_vector[a];
                           
                           // set grad Lamda
                           grad_La = gradbary_vector[b];
                           grad_Lb = gradbary_vector[c];
                           grad_Lc = gradbary_vector[a];
                           
                           // Compute Polynomials
                           x = Lb;
                           y = La + Lb;
                           Legendre_i.clear();
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           alpha = (2*i+1);
                           Int_Jacobi_j.clear();
                           x = Lc;
                           y = La + Lb + Lc;
                           poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                           
                           for(int index=0; index<4; index++)
                           {
                               double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                               u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * scalar_val;
                           }
                           
                       }
                   }
               }
           }
       }
   } // end of faces
           
    //Facets
    for(int i=0; i<=p-1;i++)
    {
        for(int j=1; j<=p-1;j++)
        {
            for(int l=1; 1<=p-1;l++)
            {
                for (int r=1; r<16; r++)
                {
                    if((i+j+l)<=p-1)
                    {
                        if (r==1)
                        {
                            // case 1: (a,b,c,d) = (0,1,2,3)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[1];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[2];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[2];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[1];
                            grad_Lb = gradbary_vector[2];
                            grad_Lc = gradbary_vector[0];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==2)
                        {
                            // case 2: (a,b,c,d) = (1,2,3,0)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[2];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[2];
                            grad_Lb = gradbary_vector[3];
                            grad_Lc = gradbary_vector[1];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==3)
                        {
                            // case 3: (a,b,c,d) = (2,3,0,1)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[0];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[3];
                            grad_Lb = gradbary_vector[0];
                            grad_Lc = gradbary_vector[2];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==4)
                        {
                            // case 4: (a,b,c,d) = (0,2,1,4)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[2];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[1];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[1];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[2];
                            grad_Lb = gradbary_vector[1];
                            grad_Lc = gradbary_vector[0];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==5)
                        {
                            // case 5: (a,b,c,d) = (2,1,4,0)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[1];
                            Lc = bary_vector[4];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[4];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[1];
                            grad_Lb = gradbary_vector[4];
                            grad_Lc = gradbary_vector[2];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==6)
                        {
                            // case 6: (a,b,c,d) = (1,4,0,2)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[4];
                            Lb = bary_vector[0];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[4];
                            grad_Lb = gradbary_vector[0];
                            grad_Lc = gradbary_vector[1];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==7)
                        {
                            // case 1: (a,b,c,d) = (0,1,3,4)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[1];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[3];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[1];
                            grad_Lb = gradbary_vector[3];
                            grad_Lc = gradbary_vector[0];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==8)
                        {
                            // case 8: (a,b,c,d) = (1,3,4,0)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[3];
                            Lc = bary_vector[4];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[3];
                            grad_Lb = gradbary_vector[4];
                            grad_Lc = gradbary_vector[1];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==9)
                        {
                            // case 9: (a,b,c,d) = (3,4,0,1)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[4];
                            Lb = bary_vector[0];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[4];
                            grad_Lb = gradbary_vector[0];
                            grad_Lc = gradbary_vector[3];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==10)
                        {
                            // case 10: (a,b,c,d) = (0,3,2,4)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[3];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[2];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[2];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[3];
                            grad_Lb = gradbary_vector[2];
                            grad_Lc = gradbary_vector[0];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==11)
                        {
                            // case 11: (a,b,c,d) = (3,2,4,0)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[2];
                            Lc = bary_vector[4];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[4];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[2];
                            grad_Lb = gradbary_vector[4];
                            grad_Lc = gradbary_vector[3];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==12)
                        {
                            // case 12: (a,b,c,d) = (2,4,0,3)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[4];
                            Lb = bary_vector[0];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[4];
                            grad_Lb = gradbary_vector[0];
                            grad_Lc = gradbary_vector[2];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==13)
                        {
                            // case 13: (a,b,c,d) = (1,2,3,4)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[2];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[2];
                            grad_Lb = gradbary_vector[3];
                            grad_Lc = gradbary_vector[1];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==14)
                        {
                            // case 1: (a,b,c,d) = (2,3,4,1)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[4];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[2];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[3];
                            grad_Lb = gradbary_vector[4];
                            grad_Lc = gradbary_vector[2];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                        if (r==15)
                        {
                            // case 15: (a,b,c,d) = (3,4,1,2)
                            
                            // Compute Family I:    <a,b,c>
                            
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[1];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[1];

                            // Compute Polynomials
                            double x = Lb;
                            double y = La + Lb;
                            std::vector<double> Legendre_i;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = (2*i+1);
                            std::vector<double> Int_Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            std::vector<double> Int_Jacobi_l;
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1]  * scalar_val;
                            }
                            
                            // Compute Family II:    <b,c,a>
                            
                            // set Lamda
                            La = bary_vector[4];
                            Lb = bary_vector[1];
                            Lc = bary_vector[3];
                            
                            // set grad Lamda
                            grad_La = gradbary_vector[4];
                            grad_Lb = gradbary_vector[1];
                            grad_Lc = gradbary_vector[3];
                            
                            // Compute Polynomials
                            x = Lb;
                            y = La + Lb;
                            Legendre_i.clear();
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            alpha = (2*i+1);
                            Int_Jacobi_j.clear();
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                            
                            alpha = 2*(i+j);
                            Int_Jacobi_l.clear();
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * scalar_val;
                            }
                        }
                    }
                }
            }
        }
    } //end of Facets
           
    //Interiors
    for (int i=0; i<=p-1; i++)
    {
        for (int j=1; j<=p-1; j++)
        {
            for (int l=1; l<=p-1; l++)
            {
                for (int m=1; m<=p-1; m++)
                {
                    if ((i+j+l+m)<=p-1)
                    {
                        for (int r=1; r<5; r++)
                        {
                            if (r==1)
                            {
                                //case where (a,b,c,d,e) = (0,1,2,3,4)
                                
                                // define lamda
                                La = bary_vector[0];
                                Lb = bary_vector[1];
                                Lc = bary_vector[2];
                                Ld = bary_vector[3];
                                Le = bary_vector[4];
                                
                                // define grad(lamda)
                                std::vector<double> grad_La = gradbary_vector[0];
                                std::vector<double> grad_Lb = gradbary_vector[1];
                                
                                // Compute Polynomials
                                double x = Lb;
                                double y = La + Lb;
                                std::vector<double> Legendre_i;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = (2*i+1);
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
                                alpha = 2*(i+j);
                                std::vector<double> Int_Jacobi_l;
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                
                                alpha = 2*(i+j+l);
                                std::vector<double> Int_Jacobi_m;
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                                
                            }
                            
                            if (r==2)
                            {
                                //case where (a,b,c,d,e) = (1,2,3,4,0)
                                
                                // define lamda
                                La = bary_vector[1];
                                Lb = bary_vector[2];
                                Lc = bary_vector[3];
                                Ld = bary_vector[4];
                                Le = bary_vector[0];
                                
                                // define grad(lamda)
                                std::vector<double> grad_La = gradbary_vector[1];
                                std::vector<double> grad_Lb = gradbary_vector[2];
                                
                                // Compute Polynomials
                                double x = Lb;
                                double y = La + Lb;
                                std::vector<double> Legendre_i;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = (2*i+1);
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
                                alpha = 2*(i+j);
                                std::vector<double> Int_Jacobi_l;
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                
                                alpha = 2*(i+j+l);
                                std::vector<double> Int_Jacobi_m;
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                            }
                            
                            if (r==3)
                            {
                                //case where (a,b,c,d,e) = (2,3,4,0,1)
                                
                                // define lamda
                                La = bary_vector[2];
                                Lb = bary_vector[3];
                                Lc = bary_vector[4];
                                Ld = bary_vector[0];
                                Le = bary_vector[1];
                                
                                // define grad(lamda)
                                std::vector<double> grad_La = gradbary_vector[2];
                                std::vector<double> grad_Lb = gradbary_vector[3];
                                
                                // Compute Polynomials
                                double x = Lb;
                                double y = La + Lb;
                                std::vector<double> Legendre_i;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = (2*i+1);
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
                                alpha = 2*(i+j);
                                std::vector<double> Int_Jacobi_l;
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                
                                alpha = 2*(i+j+l);
                                std::vector<double> Int_Jacobi_m;
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                                
                                
                            }
                            
                            if (r==4)
                            {
                                //case where (a,b,c,d,e) = (3,4,0,1,2)
                                
                                // define lamda
                                La = bary_vector[3];
                                Lb = bary_vector[4];
                                Lc = bary_vector[0];
                                Ld = bary_vector[1];
                                Le = bary_vector[2];
                                
                                // define grad(lamda)
                                std::vector<double> grad_La = gradbary_vector[3];
                                std::vector<double> grad_Lb = gradbary_vector[4];
                                
                                // Compute Polynomials
                                double x = Lb;
                                double y = La + Lb;
                                std::vector<double> Legendre_i;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = (2*i+1);
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
                                alpha = 2*(i+j);
                                std::vector<double> Int_Jacobi_l;
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                
                                alpha = 2*(i+j+l);
                                std::vector<double> Int_Jacobi_m;
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*grad_Lb[index] - Lb*grad_La[index];
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                                
                            }
                        }
                    }
                }
            }
        }
    } // end of interiors
   
    
    
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
            }

   Ti.Mult(u, shape);*/
    
}

void HSkwGrad_PentatopeElement::CalcSkwGradShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
    mfem_error("SkwGrad_PentatopeElement SkwGradShape not implemented");

}

void HSkwGrad_PentatopeElement::ProjectDivSkew(const FiniteElement& fe,
                                         ElementTransformation& Trans, DenseMatrix& DivSkew)
{
    mfem_error("SkwGrad_PentatopeElement_Fuentes::ProjectDivSkew not implemented");

}






Hcurl_PentatopeElement_Barycentric::Hcurl_PentatopeElement_Barycentric(const int p, const int type)
   : NodalFiniteElement(4, Geometry::PENTATOPE,
                        ((p + 1)*(p + 2)*(p + 3)*(p + 4))/24,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyClosed(type));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_t.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof,dim*(dim+1)/2 );
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
#endif

    int o = 0;

   // faces (see Mesh::GeneratePlanars in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,2)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[0], cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,1,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (0,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,3)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[j]/w, cp[0]);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,2,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[i]/w, cp[0], cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (1,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[p-i-j]/w, cp[0], cp[i]/w, cp[j]/w);
      }
   for (int j = 1; j < p; j++)
      for (int i=1; i + j < p; i++) // (2,3,4)
      {
         double w = cp[i] + cp[j] + cp[p-i-j];
         Nodes.IntPoint(o++).Set4(cp[0], cp[p-i-j]/w, cp[i]/w, cp[j]/w);
      }

   // facets (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,2,3)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[0]);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,2,1,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[j]/w, cp[i]/w, cp[0], cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,3,2,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[0], cp[j]/w, cp[i]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (1,2,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[p-i-j-k]/w, cp[i]/w, cp[j]/w, cp[k]/w);
         }

   // interior bubbles
   for (int l = 1; l < p; l++)
      for (int k = 1; k + l < p; k++)
         for (int j = 1; j + k + l < p; j++)
            for (int i = 1; i + j + k + l < p; i++)
            {
               double w = cp[i] + cp[j] + cp[k] + cp[l] + cp[p-i-j-k-l];
               Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[l]/w);
            }

   DenseMatrix T(dof);
    // Add code to assemble matrix

   Ti.Factor(T);
}

void Hcurl_PentatopeElement_Barycentric::CalcShape(const IntegrationPoint &ip,
                                                Vector &shape) const
{
const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
   Vector u(Dof);
#endif
    
    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);

    std::vector<double> bary_vector{L[0], L[1], L[2], L[3], 1-L[0]-L[1]-L[2]-L[3]};
    
    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};
    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> gradLa, gradLb, gradLc, gradLd;
    
    // Lamda Function for computing outter product with gradients of barycentric coordinates
    
    auto outer_prod_calc = [&](int a,int b) ->std::array<std::array<double,4>,4>{
        
        // set grad of Lamda
        gradLa = gradbary_vector[a];
        gradLb = gradbary_vector[b];
        std::array<std::array<double,4>,4> outer_prod = 
        {{{(gradLa[0]*gradLb[0]), (gradLa[0]*gradLb[1]), (gradLa[0]*gradLb[2]), (gradLa[0]*gradLb[3])},
            {(gradLa[1]*gradLb[0]), (gradLa[1]*gradLb[1]), (gradLa[1]*gradLb[2]), (gradLa[1]*gradLb[3])},
            {(gradLa[2]*gradLb[0]), (gradLa[2]*gradLb[1]), (gradLa[2]*gradLb[2]), (gradLa[2]*gradLb[3])},
            {(gradLa[3]*gradLb[0]), (gradLa[3]*gradLb[1]), (gradLa[3]*gradLb[2]), (gradLa[3]*gradLb[3])}}};
        return outer_prod;
    };
    
    
   //Faces
   for(int i=0; i<=p-1;i++)
   {
       for(int j=0;j<=p-1;j++)
       {
           for(int a=0; a<5;a++)
           {
               for(int b=0; b<5;b++)
               {
                   for(int c=0; c<5;c++)
                   {
                       if((a<b)&&(b<c)&&((i+j)<=p-1))
                       {
                           
                           // Define Barycentric coords
                           La = bary_vector[a];
                           Lb = bary_vector[b];
                           Lc = bary_vector[c];
                           
                           // Compute Polynomials
                           std::vector<double> Legendre_i;
                           double x = Lb;
                           double y = La + Lb;
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           std::vector<double> Jacobi_j;
                           double alpha = 2*i+1;
                           x = Lc;
                           y = La + Lb + Lc;
                           poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                           
                           // Compute outer produts
                           std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(b, c);
                           std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(c, b);
                           std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(c, a);
                           std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(a, c);
                           std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(a, b);
                           std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(b, a);

                           for (int row=0; row<5; row++) 
                           {
                               for (int col=0; col<5; col++) 
                               {
                                   u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                               }
                           }
                       }
                   }
               }
            }
        }
    } // end of faces
           
    //Facets
    for(int i=0; i<=p-1;i++)
    {
        for(int j=0; j<=p-1;j++)
        {
            for(int l=1; 1<=p-1;l++)
            {
                for (int r=1; r<16; r++)
                {
                    if((i+j+l)<=p-1)
                    {
                        if (r==1)
                        {
                            // case 1: (a,b,c,d) = (0,1,2,3)
                                                        
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[1];
                            Lc = bary_vector[2];
                            Ld = bary_vector[3];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[2];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(1, 2);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(2, 1);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(2, 0);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(0, 2);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(0, 1);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(1, 0);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                            
                        }
                        if (r==2)
                        {
                            // case 2: (a,b,c,d) = (1,2,3,0)
                                                        
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[2];
                            Lc = bary_vector[3];
                            Ld = bary_vector[0];
                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(2, 3);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(3, 2);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(3, 1);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(1, 3);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(1, 2);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(2, 1);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==3)
                        {
                            // case 3: (a,b,c,d) = (2,3,0,1)
                                                        
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[0];
                            Ld = bary_vector[1];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(3, 0);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(0, 3);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(0, 2);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(2, 0);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(2, 3);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(3, 2);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==4)
                        {
                            // case 4: (a,b,c,d) = (0,2,1,4)
                                                        
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[2];
                            Lc = bary_vector[1];
                            Ld = bary_vector[4];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[1];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(2, 1);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(1, 2);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(1, 0);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(0, 1);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(0, 2);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(2, 0);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==5)
                        {
                            // case 5: (a,b,c,d) = (2,1,4,0)
                                                        
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[1];
                            Lc = bary_vector[4];
                            Ld = bary_vector[0];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(1, 4);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(4, 1);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(4, 2);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(2, 4);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(2, 1);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(1, 2);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==6)
                        {
                            // case 6: (a,b,c,d) = (1,4,0,2)
                                                        
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            Ld = bary_vector[2];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(4, 0);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(0, 4);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(0, 1);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(1, 0);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(1, 4);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(4, 1);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==7)
                        {
                            // case 1: (a,b,c,d) = (0,1,3,4)
                                                        
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[1];
                            Lc = bary_vector[3];
                            Ld = bary_vector[4];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[1];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(1, 3);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(3, 1);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(3, 0);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(0, 3);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(0, 1);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(1, 0);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==8)
                        {
                            // case 8: (a,b,c,d) = (1,3,4,0)
                                                        
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[3];
                            Lc = bary_vector[4];
                            Ld = bary_vector[0];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(3, 4);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(4, 3);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(4, 1);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(1, 4);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(1, 3);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(3, 1);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==9)
                        {
                            // case 9: (a,b,c,d) = (3,4,0,1)
                                                        
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            Ld = bary_vector[1];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(4, 0);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(0, 4);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(0, 3);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(3, 0);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(3, 4);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(4, 3);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==10)
                        {
                            // case 10: (a,b,c,d) = (0,3,2,4)
                                                        
                            // set Lamda
                            La = bary_vector[0];
                            Lb = bary_vector[3];
                            Lc = bary_vector[2];
                            Ld = bary_vector[4];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[0];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[2];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(3, 2);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(2, 3);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(2, 0);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(0, 2);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(0, 3);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(3, 0);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==11)
                        {
                            // case 11: (a,b,c,d) = (3,2,4,0)
                                                        
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[2];
                            Lc = bary_vector[4];
                            Ld = bary_vector[0];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(2, 4);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(4, 2);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(4, 3);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(3, 4);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(3, 2);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(2, 3);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==12)
                        {
                            // case 12: (a,b,c,d) = (2,4,0,3)
                                                        
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[4];
                            Lc = bary_vector[0];
                            Ld = bary_vector[3];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[0];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(4, 0);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(0, 4);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(0, 2);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(2, 0);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(2, 4);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(4, 2);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==13)
                        {
                            // case 13: (a,b,c,d) = (1,2,3,4)
                                                        
                            // set Lamda
                            La = bary_vector[1];
                            Lb = bary_vector[2];
                            Lc = bary_vector[3];
                            Ld = bary_vector[4];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[1];
                            std::vector<double> grad_Lb = gradbary_vector[2];
                            std::vector<double> grad_Lc = gradbary_vector[3];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(2, 3);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(3, 2);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(3, 1);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(1, 3);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(1, 2);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(2, 1);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==14)
                        {
                            // case 1: (a,b,c,d) = (2,3,4,1)
                                                        
                            // set Lamda
                            La = bary_vector[2];
                            Lb = bary_vector[3];
                            Lc = bary_vector[4];
                            Ld = bary_vector[1];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[2];
                            std::vector<double> grad_Lb = gradbary_vector[3];
                            std::vector<double> grad_Lc = gradbary_vector[4];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(3, 4);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(4, 3);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(4, 2);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(2, 4);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(2, 3);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(3, 2);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                        if (r==15)
                        {
                            // case 15: (a,b,c,d) = (3,4,1,2)
                                                        
                            // set Lamda
                            La = bary_vector[3];
                            Lb = bary_vector[4];
                            Lc = bary_vector[1];
                            Ld = bary_vector[2];

                            
                            // set grad Lamda
                            std::vector<double> grad_La = gradbary_vector[3];
                            std::vector<double> grad_Lb = gradbary_vector[4];
                            std::vector<double> grad_Lc = gradbary_vector[1];

                            // Compute Polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            double alpha = 2*i+1;
                            x = Lc;
                            y = La + Lb + Lc;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            std::vector<double> Int_Jacobi_l;
                            alpha = 2*(i+j+1);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                            
                            // Compute outer produts
                            std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(4, 1);
                            std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(1, 4);
                            std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(1, 3);
                            std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(3, 1);
                            std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(3, 4);
                            std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(4, 3);

                            for (int row=0; row<5; row++)
                            {
                                for (int col=0; col<5; col++)
                                {
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                }
                            }
                        }
                       
                    }
                }
            }
        }
    } //end of Facets
           
    //Interiors
    for (int i=0; i<=p-1; i++)
    {
        for (int j=0; j<=p-1; j++)
        {
            for (int l=1; l<=p-1; l++)
            {
                for (int m=1; m<=p-1; m++)
                {
                    if ((i+j+l+m)<=p-1)
                    {
                        for (int r=1; r<7; r++)
                        {
                            if (r==1)
                            {
                                //case where (a,b,c,d,e) = (0,1,2,3,4)
                                
                                // set Lamda
                                La = bary_vector[0];
                                Lb = bary_vector[1];
                                Lc = bary_vector[2];
                                Ld = bary_vector[3];
                                Le = bary_vector[4];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4>  outer_prod_bc = outer_prod_calc(1, 2);
                                std::array<std::array<double,4>,4>  outer_prod_cb = outer_prod_calc(2, 1);
                                std::array<std::array<double,4>,4>  outer_prod_ca = outer_prod_calc(2, 0);
                                std::array<std::array<double,4>,4>  outer_prod_ac = outer_prod_calc(0, 2);
                                std::array<std::array<double,4>,4>  outer_prod_ab = outer_prod_calc(0, 1);
                                std::array<std::array<double,4>,4>  outer_prod_ba = outer_prod_calc(1, 0);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                                
                            }
                            
                            if (r==2)
                            {
                                //case where (a,b,c,d,e) = (1,2,3,4,0)
                                
                                // set Lamda
                                La = bary_vector[1];
                                Lb = bary_vector[2];
                                Lc = bary_vector[3];
                                Ld = bary_vector[4];
                                Le = bary_vector[0];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(2, 3);
                                std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(3, 2);
                                std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(3, 1);
                                std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(1, 3);
                                std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(1, 2);
                                std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(2, 1);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                            }
                            
                            if (r==3)
                            {
                                //case where (a,b,c,d,e) = (2,3,4,0,1)
                                
                                // set Lamda
                                La = bary_vector[2];
                                Lb = bary_vector[3];
                                Lc = bary_vector[4];
                                Ld = bary_vector[0];
                                Le = bary_vector[1];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(3, 4);
                                std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(4, 3);
                                std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(4, 2);
                                std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(2, 4);
                                std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(2, 3);
                                std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(3, 2);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                            }
                            
                            if (r==4)
                            {
                                //case where (a,b,c,d,e) = (3,4,0,1,2)
                                
                                // set Lamda
                                La = bary_vector[3];
                                Lb = bary_vector[4];
                                Lc = bary_vector[0];
                                Ld = bary_vector[1];
                                Le = bary_vector[2];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(4, 0);
                                std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(0, 4);
                                std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(0, 3);
                                std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(3, 0);
                                std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(3, 4);
                                std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(4, 3);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                            }
                            if (r==5) {
                                
                                // case where (a,b,c,d,e) = (4,0,1,2,3)
                                
                                // set Lamda
                                La = bary_vector[4];
                                Lb = bary_vector[0];
                                Lc = bary_vector[1];
                                Ld = bary_vector[2];
                                Le = bary_vector[3];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(0, 1);
                                std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(1, 0);
                                std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(1, 4);
                                std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(4, 1);
                                std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(4, 0);
                                std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(0, 4);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                            }
                            if (r==6) {
                                
                                // case where (a,b,c,d,e) = (0,1,3,2,4)
                                
                                // set Lamda
                                La = bary_vector[0];
                                Lb = bary_vector[1];
                                Lc = bary_vector[3];
                                Ld = bary_vector[2];
                                Le = bary_vector[4];
                                
                                // Compute Polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                std::vector<double> Int_Jacobi_l;
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                poly1d.CalcIntJacobi(l, x, y, alpha, Int_Jacobi_l);
                                
                                std::vector<double> Int_Jacobi_m;
                                alpha = 2*(i+j+l);
                                x = Le;
                                y = 1;
                                poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                // Compute outer produts
                                std::array<std::array<double,4>,4> outer_prod_bc = outer_prod_calc(1, 3);
                                std::array<std::array<double,4>,4> outer_prod_cb = outer_prod_calc(3, 1);
                                std::array<std::array<double,4>,4> outer_prod_ca = outer_prod_calc(3, 0);
                                std::array<std::array<double,4>,4> outer_prod_ac = outer_prod_calc(0, 3);
                                std::array<std::array<double,4>,4> outer_prod_ab = outer_prod_calc(0, 1);
                                std::array<std::array<double,4>,4> outer_prod_ba = outer_prod_calc(1, 0);

                                for (int row=0; row<5; row++)
                                {
                                    for (int col=0; col<5; col++)
                                    {
                                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Int_Jacobi_l[Int_Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] *  (La*(outer_prod_bc[row][col] - outer_prod_cb[row][col]) + Lb*(outer_prod_ca[row][col] - outer_prod_ac[row][col]) + Lc*(outer_prod_ab[row][col] - outer_prod_ba[row][col]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } // end of interiors
   
    
    
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
            }

   Ti.Mult(u, shape);
}

void Hcurl_PentatopeElement_Barycentric::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
    //Do Nothing
    mfem_error("Hcurl_PentatopeElement_Barycentric::CalcDShape not implemented");

}

void Hcurl_PentatopeElement_Barycentric::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
    //Do Nothing
    mfem_error("Hcurl_PentatopeElement_Barycentric::CalcHessian not implemented");
}

const double Hdiv_PentatopeElement::nk[20] =
{ 0,0,0,-1,  0,0,-1,0,  0,-1,0,0,  -1,0,0,0,  1,1,1,1};
// { .5,.5,.5, -.5,0,0, 0,-.5,0, 0,0,-.5}; // n_F |F|

const double Hdiv_PentatopeElement::c = 1./5.;

Hdiv_PentatopeElement::Hdiv_PentatopeElement(const int p)
   : VectorFiniteElement(4, Geometry::PENTATOPE, (p + 1)*(p + 2)*(p + 3)*(p + 5)/6,
                         p + 1, H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   const double *iop = (p > 0) ? poly1d.OpenPoints(p - 1) : NULL;
   const double *bop = poly1d.OpenPoints(p);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim);
   divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
#endif
   //std::cout << "p_V = " << p << std::endl;
   int o = 0;
   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp,
   //        the constructor of H1_PentatopeElement)
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,2,3)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, bop[j]/w, bop[k]/w, 0.);
            dof2nk[o++] = 0;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,2,1,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[j]/w, bop[i]/w, 0., bop[k]/w);
            dof2nk[o++] = 1;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, 0., bop[j]/w, bop[k]/w);
            dof2nk[o++] = 2;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,3,2,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(0., bop[j]/w, bop[i]/w, bop[k]/w);
            dof2nk[o++] = 3;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (1,2,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[p-i-j-k]/w, bop[i]/w, bop[j]/w, bop[k]/w);
            dof2nk[o++] = 4;
         }

   // interior
   for (int l = 0; l < p; l++)
      for (int k = 0; k + l < p; k++)
         for (int j = 0; j + k + l < p; j++)
            for (int i = 0; i + j + k + l < p; i++)
            {
                double w = iop[i] + iop[j] + iop[k] + iop[l] + iop[p-1-i-j-k-l];
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 1;
                dof2nk[o++] = 3;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 2;
                dof2nk[o++] = 2;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 3;
                dof2nk[o++] = 1;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                // old def dof2nk[o++] = 4;
                dof2nk[o++] = 0;
            }

   DenseMatrix T(dof);
   DenseMatrix B(dof, dim);
   int reset_o;
   //std::cout << "num_dof =" << num_dof << std::endl;

   for (int q = 0; q < dof; q++)
   {
       const IntegrationPoint &ip = Nodes.IntPoint(q);
       //std::cout << "Dof" << q << " = " << ip.x << ", " << ip.y << ", " << ip.z << ", " << ip.t << std::endl;

       //const double *nm = nk + 4*dof2nk[m];
       
       const Vector nm({nk[4*dof2nk[q]], nk[4*dof2nk[q]+1], nk[4*dof2nk[q]+2], nk[4*dof2nk[q]+3]});
       
       int o = 0;
       int num_dof = dof;

       //compute barycentric coordinates as function of ip
       std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};
       
       // compute the gradient of the barycentric coords
       std::vector<double> gradL1{1,0,0,0};
       std::vector<double> gradL2{0,1,0,0};
       std::vector<double> gradL3{0,0,1,0};
       std::vector<double> gradL4{0,0,0,1};
       std::vector<double> gradL5{-1,-1,-1,-1};

       std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
                  
       double La, Lb, Lc, Ld, Le;
       
       std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
       
       // Lamda Function for computing cross product with gradients of barycentric coordinates
       
       auto four_cross = [&](int A,int B,int C) -> std::vector<double>{
           
           // set grad of Lamda
           //std::cout << "vector info: " << A << ", " << B << ", " << C << std::endl;
           grad_La = gradbary_vector[A];
           //std::cout << grad_La[0] << ", " << grad_La[1] << ", " << grad_La[2] << ", " << grad_La[3] << std::endl;
           grad_Lb = gradbary_vector[B];
           //std::cout << grad_Lb[0] << ", " << grad_Lb[1] << ", " << grad_Lb[2] << ", " << grad_Lb[3] << std::endl;
           grad_Lc = gradbary_vector[C];
           //std::cout << grad_Lc[0] << ", " << grad_Lc[1] << ", " << grad_Lc[2] << ", " << grad_Lc[3] << std::endl;

                   
           // set variables
           double aa,bb,cc,dd,ee,ff,gg,hh,mm,nn,pp,qq;
           aa = grad_La[0];
           bb = grad_La[1];
           cc = grad_La[2];
           dd = grad_La[3];
           
           ee = grad_Lb[0];
           ff = grad_Lb[1];
           gg = grad_Lb[2];
           hh = grad_Lb[3];
           
           mm = grad_Lc[0];
           nn = grad_Lc[1];
           qq = grad_Lc[2];
           pp = grad_Lc[3];
           
           double i_hat, j_hat, k_hat, l_hat;
           
           i_hat = bb*((gg*pp)-(qq*hh)) - cc*((ff*pp)-(nn*hh)) + dd*((ff*qq)-(gg*nn));
           
           j_hat = aa*((gg*pp)-(qq*hh)) - cc*((ee*pp)-(mm*hh)) + dd*((ee*qq)-(gg*mm));
           
           k_hat = aa*((ff*pp)-(nn*hh)) - bb*((ee*pp)-(mm*hh)) + dd*((ee*nn)-(ff*mm));
           
           l_hat = aa*((ff*qq)-(nn*gg)) - bb*((ee*qq)-(mm*gg)) + cc*((ee*nn)-(mm*ff));
           
           //std::cout << i_hat << ", " << -1.0*j_hat << ", " << k_hat << ", " << -1.0*l_hat << std::endl;
           
           std::vector<double> four_cross = {i_hat, -1.0*j_hat, k_hat, -1.0*l_hat};
           return four_cross;

       }; // end of lamda function
       
              
       //Facets
       int a;
       int b;
       int c;
       int d;
       int e;
         
       for (int f=0; f<5; f++)
       {
           // Define each facet
           if (f==0)
           {
               a = 0;
               b = 1;
               c = 2;
               d = 3;
               
           }
           else if(f==1)
           {
               // Define Facet
               a = 0;
               b = 2;
               c = 1;
               d = 4;
               
           }
           // Define each facet
           else if (f==2)
           {
               // Define Facet
               a = 0;
               b = 1;
               c = 3;
               d = 4;
               
           }
           // Define each facet
           else if (f==3)
           {
               // Define Facet
               a = 0;
               b = 3;
               c = 2;
               d = 4;
               
           }
           // Define each facet
           else if (f==4)
           {
               // Define Facet
               a = 1;
               b = 2;
               c = 3;
               d = 4;
               
           }
           else
           {
               mfem_error("Invaild facet");
           }
           
           // Define Barycentric Coordinates
           La = bary_vector[a];
           Lb = bary_vector[b];
           Lc = bary_vector[c];
           Ld = bary_vector[d];
           
           // set grad of Lamda
           grad_La = gradbary_vector[a];
           grad_Lb = gradbary_vector[b];
           grad_Lc = gradbary_vector[c];
           grad_Ld = gradbary_vector[d];
           
           for(int i=0; i<=(p);i++)
           {
               for(int j=0; j<=(p);j++)
               {
                   for(int l=0; l<=(p);l++)
                   {
                       if((i+j+l)<=(p))
                       {
                           
                           // compute polynomials
                           
                           std::vector<double> Legendre_i;
                           double x = Lb;
                           double y = La + Lb;
                           poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                           
                           
                           double alpha = (2.0*i)+1.0;
                           x = Lc;
                           y = La + Lb + Lc;
                           std::vector<double> Jacobi_j;
                           poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                           
                           
                           alpha = 2.0*(i+j+1.0);
                           x = Ld;
                           y = La + Lb + Lc + Ld;
                           std::vector<double> Jacobi_l;
                           poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                           
                           
                           //Compute relevant four curl
                           std::vector<double> curl_bcd = four_cross(b,c,d);
                           std::vector<double> curl_cda = four_cross(c,d,a);
                           std::vector<double> curl_dab = four_cross(d,a,b);
                           std::vector<double> curl_abc = four_cross(a,b,c);
                           
                           // Compute Vector Components
                           double scalar_val_x = La*curl_bcd[0] - Lb*curl_cda[0] + Lc*curl_dab[0] - Ld*curl_abc[0];
                           double scalar_val_y = La*curl_bcd[1] - Lb*curl_cda[1] + Lc*curl_dab[1] - Ld*curl_abc[1];
                           double scalar_val_z = La*curl_bcd[2] - Lb*curl_cda[2] + Lc*curl_dab[2] - Ld*curl_abc[2];
                           double scalar_val_t = La*curl_bcd[3] - Lb*curl_cda[3] + Lc*curl_dab[3] - Ld*curl_abc[3];
                           
                           // Create Basis
                           
                           B(o ,0) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val_x;
                           
                           B(o ,1) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val_y;
                           
                           B(o ,2) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val_z;
                           
                           B(o ,3) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val_t;
                           
                           o++;

                       }
                   }
               }
           }
           
       }//end of Facets
       
              
       //Interiors
       for (int r = 0; r<4; r++)
       {
           if (r ==0)
           {
               //case where (a,b,c,d,e) = (0,1,2,3,4)
               // define lamda
               La = bary_vector[0];
               Lb = bary_vector[1];
               Lc = bary_vector[2];
               Ld = bary_vector[3];
               Le = bary_vector[4];
               
               // define a,b,c,d
               a = 0;
               b = 1;
               c = 2;
               d = 3;
               
           }
           else if (r==1)
           {
               //case where (a,b,c,d,e) = (1,2,3,4,0)
               // define lamda
               La = bary_vector[1];
               Lb = bary_vector[2];
               Lc = bary_vector[3];
               Ld = bary_vector[4];
               Le = bary_vector[0];
               
               // define a,b,c,d
               a = 1;
               b = 2;
               c = 3;
               d = 4;
               
               
           }
           else if (r==2)
           {
               //case where (a,b,c,d,e) = (2,3,4,0,1)
               // define lamda
               La = bary_vector[2];
               Lb = bary_vector[3];
               Lc = bary_vector[4];
               Ld = bary_vector[0];
               Le = bary_vector[1];
               
               // define a,b,c,d
               a = 2;
               b = 3;
               c = 4;
               d = 0;
               
           }
           else if (r==3)
           {
               //case where (a,b,c,d,e) = (3,4,0,1,2)
               // define lamda
               La = bary_vector[3];
               Lb = bary_vector[4];
               Lc = bary_vector[0];
               Ld = bary_vector[1];
               Le = bary_vector[2];
               
               // define a,b,c,d
               a = 3;
               b = 4;
               c = 0;
               d = 1;
               
           }
           else
           {
               mfem_error("Invaild Bubble");
           }
           
           for (int i=0; i<=(p+1); i++)
           {
               for (int j=0; j<=(p+1); j++)
               {
                   for (int l=0; l<=(p+1); l++)
                   {
                       for (int m=1; m<=(p+1); m++)
                       {
                           if ((i+j+l+m)<=((p+1)-1))
                           {
                               // compute polynomials
                               //std::cout << "inside bubble " << o << std::endl;
                               std::vector<double> Legendre_i;
                               double x = Lb;
                               double y = La + Lb;
                               poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                               //std::cout << "Legendre = " << Legendre_i[Legendre_i.size()-1] << std::endl;
                               
                               double alpha = 2*i+1;
                               x = Lc;
                               y = La + Lb + Lc;
                               std::vector<double> Jacobi_j;
                               poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                               //std::cout << "Jacobi_1 = " << Jacobi_j[Jacobi_j.size()-1] << std::endl;
                               
                               alpha = 2*(i+j+1);
                               x = Ld;
                               y = La + Lb + Lc + Ld;
                               std::vector<double> Jacobi_l;
                               poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                               //std::cout << "Jacobi_2 = " << Jacobi_l[Jacobi_l.size()-1] << std::endl;

                               
                               alpha = 2*(i+j+l)+3;
                               x = Le;
                               y = 1;
                               std::vector<double> Int_Jacobi_m;
                               poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                               //std::cout << "Jacobi_3 = " << Int_Jacobi_m[Int_Jacobi_m.size()-1] << std::endl;
                               //std::cout << "x = " << x << std::endl;

                               
                               //Compute relevant four curl
                               std::vector<double> curl_bcd = four_cross(b,c,d);
                               std::vector<double> curl_cda = four_cross(c,d,a);
                               std::vector<double> curl_dab = four_cross(d,a,b);
                               std::vector<double> curl_abc = four_cross(a,b,c);
                               
                               // Compute Vector Components
                               double scalar_val_x = La*curl_bcd[0] - Lb*curl_cda[0] + Lc*curl_dab[0] - Ld*curl_abc[0];
                               double scalar_val_y = La*curl_bcd[1] - Lb*curl_cda[1] + Lc*curl_dab[1] - Ld*curl_abc[1];
                               double scalar_val_z = La*curl_bcd[2] - Lb*curl_cda[2] + Lc*curl_dab[2] - Ld*curl_abc[2];
                               double scalar_val_t = La*curl_bcd[3] - Lb*curl_cda[3] + Lc*curl_dab[3] - Ld*curl_abc[3];
                               
                               // Create Basis
                               
                               B(o ,0) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1]  * scalar_val_x;
                               
                               B(o ,1) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1]  * scalar_val_y;
                               
                               B(o ,2) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val_z;
                               
                               B(o ,3) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val_t;
                               
                               o++;
                               
                           }
                       }
                   }
               }
           }
           
       }// End of Interiors
       
       
       B.Mult(nm, T.GetColumn(q));
//      std::ofstream Vander_file("Vander_matrix_Fuentes.txt");
//      T.PrintMatlab(Vander_file);
       
    

   }

//   std::ofstream Vander_file("Vander_matrix_Fuentes.txt");
//   T.PrintMatlab(Vander_file);
//    
//   std::ofstream basis_file("Basis_matrix.txt");
//   B.PrintMatlab(basis_file);
     
    


   Ti.Factor(T);
    
//   std::ofstream VanderInv_file("VanderInv_matrix.txt");
//   Ti.PrintMatlab(VanderInv_file);

   mfem::out << "RT_PentatopeElement(" << p << ") : "; Ti.TestInversion();
}

void Hdiv_PentatopeElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   const int p = order - 1;
   //std::cout << "p_CalcShape = " << p << std::endl;
    int reset_o;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
   DenseMatrix u(Dof, Dim);
#endif
    //ofstream logfile;
    //logfile.open ("logfile_shape.txt", "a");
    ofstream logfile("logfile_shape.txt", std::ios_base::app);
    //compute barycentric coordinates as function of ip
    std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};

    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};

    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
    int o = 0;
    int local_dof = 0;

    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    
    // Lamda Function for computing cross product with gradients of barycentric coordinates
    
    auto four_cross = [&](int A,int B,int C) -> std::vector<double>{
        
        // set grad of Lamda
        grad_La = gradbary_vector[A];
        grad_Lb = gradbary_vector[B];
        grad_Lc = gradbary_vector[C];
                
        // set variables
        double aa,bb,cc,dd,ee,ff,gg,hh,mm,nn,pp,qq;
        aa = grad_La[0];
        bb = grad_La[1];
        cc = grad_La[2];
        dd = grad_La[3];
        
        ee = grad_Lb[0];
        ff = grad_Lb[1];
        gg = grad_Lb[2];
        hh = grad_Lb[3];
        
        mm = grad_Lc[0];
        nn = grad_Lc[1];
        qq = grad_Lc[2];
        pp = grad_Lc[3];
        
        double i_hat, j_hat, k_hat, l_hat;
        
        i_hat = bb*((gg*pp)-(qq*hh)) - cc*((ff*pp)-(nn*hh)) + dd*((ff*qq)-(gg*nn));
        
        j_hat = aa*((gg*pp)-(qq*hh)) - cc*((ee*pp)-(mm*hh)) + dd*((ee*qq)-(gg*mm));
        
        k_hat = aa*((ff*pp)-(nn*hh)) - bb*((ee*pp)-(mm*hh)) + dd*((ee*nn)-(ff*mm));
        
        l_hat = aa*((ff*qq)-(nn*gg)) - bb*((ee*qq)-(mm*gg)) + cc*((ee*nn)-(mm*ff));
        
        std::vector<double> four_cross = {i_hat, -1.0*j_hat, k_hat, -1.0*l_hat};
        
//        double dotprod_La = four_cross[0]*grad_La[0] + four_cross[1]*grad_La[1] + four_cross[2]*grad_La[2] + four_cross[3]*grad_La[3];
//        double dotprod_Lb = four_cross[0]*grad_Lb[0] + four_cross[1]*grad_Lb[1] + four_cross[2]*grad_Lb[2] + four_cross[3]*grad_Lb[3];
//        double dotprod_Lc = four_cross[0]*grad_Lc[0] + four_cross[1]*grad_Lc[1] + four_cross[2]*grad_Lc[2] + four_cross[3]*grad_Lc[3];
//        
//        std::cout << "dotprod_La = " << dotprod_La << std::endl;
//        std::cout << "dotprod_Lb = " << dotprod_Lb << std::endl;
//        std::cout << "dotprod_Lc = " << dotprod_Lc << std::endl;

        
        return four_cross;

    }; // end of lamda function
    
    //Facets
    int a;
    int b;
    int c;
    int d;
    int e;
      
    for (int f=0; f<5; f++)
    {
        // Define each facet
        if (f==0)
        {
            a = 0;
            b = 1;
            c = 2;
            d = 3;
        }
        else if(f==1)
        {
            // Define Facet
            a = 0;
            b = 2;
            c = 1;
            d = 4;
        }
        // Define each facet
        else if (f==2)
        {
            // Define Facet
            a = 0;
            b = 1;
            c = 3;
            d = 4;
        }
        // Define each facet
        else if (f==3)
        {
            // Define Facet
            a = 0;
            b = 3;
            c = 2;
            d = 4;
        }
        // Define each facet
        else if (f==4)
        {
            // Define Facet
            a = 1;
            b = 2;
            c = 3;
            d = 4;
        }
        else
        {
            mfem_error("Invaild facet");
        }
        
        // Define Barycentric Coordinates
        La = bary_vector[a];
        Lb = bary_vector[b];
        Lc = bary_vector[c];
        Ld = bary_vector[d];
        
        // set grad of Lamda
        grad_La = gradbary_vector[a];
        grad_Lb = gradbary_vector[b];
        grad_Lc = gradbary_vector[c];
        grad_Ld = gradbary_vector[d];
        
        //std::cout << "Section" << std::endl;
        
        for(int i=0; i<=(p);i++)
        {
            for(int j=0; j<=(p);j++)
            {
                for(int l=0; l<=(p);l++)
                {
                    if((i+j+l)<=(p))
                    {

                        // compute polynomials
                        std::vector<double> Legendre_i;
                        double x = Lb;
                        double y = La + Lb;
                        poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                        
                        double alpha = 2.0*i+1.0;
                        x = Lc;
                        y = La + Lb + Lc;
                        std::vector<double> Jacobi_j;
                        poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                        
                        alpha = 2.0*(i+j+1.0);
                        x = Ld;
                        y = La + Lb + Lc + Ld;
                        std::vector<double> Jacobi_l;
                        poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                        
                        //Compute relevant four curl
                        std::vector<double> curl_bcd = four_cross(b,c,d);
                        std::vector<double> curl_cda = four_cross(c,d,a);
                        std::vector<double> curl_dab = four_cross(d,a,b);
                        std::vector<double> curl_abc = four_cross(a,b,c);
                        
                        //std::cout << "Before Lamda" << std::endl;
                        for(int index=0; index<4; index++)
                        {
                            double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                            // Compute Facet Funtion
                            u(o, index) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val;
                        }
                        o++;
                    }
                }
            }
        }
                        

    }//end of Facets

    
    
    //Interiors
    for (int r = 0; r<4; r++)
    {
        if (r ==0)
        {
            //case where (a,b,c,d,e) = (0,1,2,3,4)
            // define lamda
            La = bary_vector[0];
            Lb = bary_vector[1];
            Lc = bary_vector[2];
            Ld = bary_vector[3];
            Le = bary_vector[4];
            
            // define a,b,c,d
            a = 0;
            b = 1;
            c = 2;
            d = 3;
            
        }
        else if (r==1)
        {
            //case where (a,b,c,d,e) = (1,2,3,4,0)
            // define lamda
            La = bary_vector[1];
            Lb = bary_vector[2];
            Lc = bary_vector[3];
            Ld = bary_vector[4];
            Le = bary_vector[0];
            
            // define a,b,c,d
            a = 1;
            b = 2;
            c = 3;
            d = 4;
            
            
        }
        else if (r==2)
        {
            //case where (a,b,c,d,e) = (2,3,4,0,1)
            // define lamda
            La = bary_vector[2];
            Lb = bary_vector[3];
            Lc = bary_vector[4];
            Ld = bary_vector[0];
            Le = bary_vector[1];
            
            // define a,b,c,d
            a = 2;
            b = 3;
            c = 4;
            d = 0;
            
        }
        else if (r==3)
        {
            //case where (a,b,c,d,e) = (3,4,0,1,2)
            // define lamda
            La = bary_vector[3];
            Lb = bary_vector[4];
            Lc = bary_vector[0];
            Ld = bary_vector[1];
            Le = bary_vector[2];
            
            // define a,b,c,d
            a = 3;
            b = 4;
            c = 0;
            d = 1;
            
        }
        else
        {
            mfem_error("Invaild Bubble");
        }
        
        for (int i=0; i<=p; i++)
        {
            for (int j=0; j<=p; j++)
            {
                for (int l=0; l<=p; l++)
                {
                    for (int m=1; m<=p; m++)
                    {
                        if ((i+j+l+m)<=p)
                        {
                            // compute polynomials
                            
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            double alpha = 2.0*i + 1.0;
                            x = Lc;
                            y = La + Lb + Lc;
                            std::vector<double> Jacobi_j;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            alpha = 2.0*(i+j+1.0);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            std::vector<double> Jacobi_l;
                            poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                            
                            alpha = 2.0*(i+j+l)+3.0;
                            x = Le;
                            y = 1.0;
                            std::vector<double> Int_Jacobi_m;
                            poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);

                            //Compute relevant four curl
                            std::vector<double> curl_bcd = four_cross(b,c,d);
                            std::vector<double> curl_cda = four_cross(c,d,a);
                            std::vector<double> curl_dab = four_cross(d,a,b);
                            std::vector<double> curl_abc = four_cross(a,b,c);
                            
                            for(int index=0; index<4; index++)
                            {
                                double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                                // Compute Facet Funtion
                                u(o, index) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                            }
                            o++;
                            local_dof = o;
                        }
                    }
                }
            }
        }
        
    }// End of Interiors
    

    


   Ti.Mult(u, shape);
    
//   std::ofstream shape_file("shape_matrix_Fuentes.txt");
//   shape.PrintMatlab(shape_file);
//    
//   std::ofstream u_file("u_matrix.txt");
//   u.PrintMatlab(u_file);

   //std::cout << "Made it to Shape" << std::endl;
}

void Hdiv_PentatopeElement::CalcDivShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
   const int p = order - 1;
   //std::cout << "Value of p = " << p << std::endl;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   Vector divu(Dof);
#endif
    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);

    std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};

    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};

    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    
    // Lamda Function for computing cross product with gradients of barycentric coordinates
    
    auto four_cross = [&](int A,int B,int C) -> std::vector<double>{
        
        // set grad of Lamda
        grad_La = gradbary_vector[A];
        grad_Lb = gradbary_vector[B];
        grad_Lc = gradbary_vector[C];
                
        // set variables
        double aa,bb,cc,dd,ee,ff,gg,hh,mm,nn,pp,qq;
        aa = grad_La[0];
        bb = grad_La[1];
        cc = grad_La[2];
        dd = grad_La[3];
        
        ee = grad_Lb[0];
        ff = grad_Lb[1];
        gg = grad_Lb[2];
        hh = grad_Lb[3];
        
        mm = grad_Lc[0];
        nn = grad_Lc[1];
        qq = grad_Lc[2];
        pp = grad_Lc[3];
        
        double i_hat, j_hat, k_hat, l_hat;
        
        i_hat = bb*((gg*pp)-(qq*hh)) - cc*((ff*pp)-(nn*hh)) + dd*((ff*qq)-(gg*nn));
        
        j_hat = aa*((gg*pp)-(qq*hh)) - cc*((ee*pp)-(mm*hh)) + dd*((ee*qq)-(gg*mm));
        
        k_hat = aa*((ff*pp)-(nn*hh)) - bb*((ee*pp)-(mm*hh)) + dd*((ee*nn)-(ff*mm));
        
        l_hat = aa*((ff*qq)-(nn*gg)) - bb*((ee*qq)-(mm*gg)) + cc*((ee*nn)-(mm*ff));
        
        std::vector<double> four_cross = {i_hat, -1.0*j_hat, k_hat, -1.0*l_hat};
        return four_cross;

    }; // end of lamda function
    
           
    //Facets
    int a;
    int b;
    int c;
    int d;
    int e;
      
    for (int f=0; f<5; f++)
    {
        // Define each facet
        if (f==0)
        {
            a = 0;
            b = 1;
            c = 2;
            d = 3;
        }
        else if(f==1)
        {
            // Define Facet
            a = 0;
            b = 2;
            c = 1;
            d = 4;
        }
        // Define each facet
        else if (f==2)
        {
            // Define Facet
            a = 0;
            b = 1;
            c = 3;
            d = 4;
        }
        // Define each facet
        else if (f==3)
        {
            // Define Facet
            a = 0;
            b = 3;
            c = 2;
            d = 4;
        }
        // Define each facet
        else if (f==4)
        {
            // Define Facet
            a = 1;
            b = 2;
            c = 3;
            d = 4;
        }
        else
        {
            mfem_error("Invaild facet");
        }
        
        // Define Barycentric Coordinates
        La = bary_vector[a];
        Lb = bary_vector[b];
        Lc = bary_vector[c];
        Ld = bary_vector[d];
        
        // set grad of Lamda
        grad_La = gradbary_vector[a];
        grad_Lb = gradbary_vector[b];
        grad_Lc = gradbary_vector[c];
        grad_Ld = gradbary_vector[d];
        
        for(int i=0; i<=p;i++)
        {
            for(int j=0; j<=p;j++)
            {
                for(int l=0; l<=p;l++)
                {
                    //std::cout << "i+j+l = " << i+j+l << std::endl;
                    if((i+j+l)<=p)
                    {
                        
                        // compute polynomials
                        //std::cout << "Inside Shape" << std::endl;

                        
                        // Legendre Poly
                        std::vector<double> Legendre_i;
                        std::vector<double> Legendre_i_ref;
                        double x = Lb;
                        double y = La + Lb;
                        poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                        // Derivative Legendre Poly
                        std::vector<double> Legendre_i_dx;
                        std::vector<double> Legendre_i_dt;
                        poly1d.CalcScaledLegendreDerivative(i, x, y, Legendre_i_ref, Legendre_i_dx, Legendre_i_dt);
                        
                        // Jacobi Poly
                        double alpha = 2.0*i+1.0;
                        x = Lc;
                        y = La + Lb + Lc;
                        std::vector<double> Jacobi_j;
                        std::vector<double> Jacobi_j_ref;
                        poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                        // Derivative Jacobi Poly
                        std::vector<double> Jacobi_j_dx;
                        std::vector<double> Jacobi_j_dt;
                        poly1d.CalcScaledJacobiDerivative(j, alpha, x, y, Jacobi_j_ref, Jacobi_j_dx, Jacobi_j_dt);
                        
                        // Jacobi Poly
                        alpha = 2.0*(i+j+1.0);
                        x = Ld;
                        y = La + Lb + Lc + Ld;
                        std::vector<double> Jacobi_l;
                        std::vector<double> Jacobi_l_ref;
                        poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                        // Derivative Jacobi Poly
                        std::vector<double> Jacobi_l_dx;
                        std::vector<double> Jacobi_l_dt;
                        poly1d.CalcScaledJacobiDerivative(l, alpha, x, y, Jacobi_l_ref, Jacobi_l_dx, Jacobi_l_dt);
                        
                        //Compute relevant four curl
                        std::vector<double> curl_bcd = four_cross(b,c,d);
                        std::vector<double> curl_cda = four_cross(c,d,a);
                        std::vector<double> curl_dab = four_cross(d,a,b);
                        std::vector<double> curl_abc = four_cross(a,b,c);
                        
                        // Basis Functions can be decomposed into vector and scalar parts
                        
                        // Scalar part
                        double Scalar = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1];
                        
                        // d Scalar / dx
                        double dscalar_x = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[0] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[0] + grad_Lb[0]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[0] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[0] + grad_Lb[0] + grad_Lc[0]))*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[0] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[0] + grad_Lb[0] + grad_Lc[0] + grad_Ld[0]));
                        
                        //std::cout << "dScalar / dx = " << dscalar_x << std::endl;
                        
                        
                        // d Scalar / dy
                        double dscalar_y = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[1] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[1] + grad_Lb[1]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[1] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[1] + grad_Lb[1] + grad_Lc[1]))*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[1] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[1] + grad_Lb[1] + grad_Lc[1] + grad_Ld[1]));
                        
                        //std::cout << "dScalar / dy = " << dscalar_y << std::endl;

                        
                        
                        // d Scalar / dz
                        double dscalar_z = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[2] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[2] + grad_Lb[2]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[2] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[2] + grad_Lb[2] + grad_Lc[2]))*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[2] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[2] + grad_Lb[2] + grad_Lc[2] + grad_Ld[2]));
                        
                        //std::cout << "dScalar / dz = " << dscalar_z << std::endl;

                        
                        
                        // d Scalar / dt
                        double dscalar_t = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[3] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[3] + grad_Lb[3]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[3] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[3] + grad_Lb[3] + grad_Lc[3]))*Jacobi_l[Jacobi_l.size()-1] +
                        
                            Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[3] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[3] + grad_Lb[3] + grad_Lc[3] + grad_Ld[3]));
                        
                        //std::cout << "dScalar / dt = " << dscalar_t << std::endl;

                        
                        
                        // components of vectors part
                        double vec_val_x = La*curl_bcd[0] - Lb*curl_cda[0] + Lc*curl_dab[0] - Ld*curl_abc[0];
                        
                        double vec_val_y = La*curl_bcd[1] - Lb*curl_cda[1] + Lc*curl_dab[1] - Ld*curl_abc[1];

                        double vec_val_z = La*curl_bcd[2] - Lb*curl_cda[2] + Lc*curl_dab[2] - Ld*curl_abc[2];

                        double vec_val_t = La*curl_bcd[3] - Lb*curl_cda[3] + Lc*curl_dab[3] - Ld*curl_abc[3];
                        
                        // Each term of div(vector)
                        double div_val_x = grad_La[0]*curl_bcd[0] - grad_Lb[0]*curl_cda[0] + grad_Lc[0]*curl_dab[0] - grad_Ld[0]*curl_abc[0];
                        
                        double div_val_y = grad_La[1]*curl_bcd[1] - grad_Lb[1]*curl_cda[1] + grad_Lc[1]*curl_dab[1] - grad_Ld[1]*curl_abc[1];

                        double div_val_z = grad_La[2]*curl_bcd[2] - grad_Lb[2]*curl_cda[2] + grad_Lc[2]*curl_dab[2] - grad_Ld[2]*curl_abc[2];

                        double div_val_t = grad_La[3]*curl_bcd[3] - grad_Lb[3]*curl_cda[3] + grad_Lc[3]*curl_dab[3] - grad_Ld[3]*curl_abc[3];
                        
                        // Div value
                        double div_vec = div_val_x + div_val_y + div_val_z + div_val_t;
                        
                        
                        // Add Divergence of Funtion
                        divu(o++) = (dscalar_x * vec_val_x + dscalar_y * vec_val_y + dscalar_z * vec_val_z + dscalar_t * vec_val_t) + (div_vec * Scalar);

                    }
                }
            }
        }
    }//end of Facets
    
           
    //Interiors
    for (int r = 0; r<4; r++)
    {
        if (r ==0)
        {
            //case where (a,b,c,d,e) = (0,1,2,3,4)
            // define lamda
            La = bary_vector[0];
            Lb = bary_vector[1];
            Lc = bary_vector[2];
            Ld = bary_vector[3];
            Le = bary_vector[4];
            
            // define grad(lamda)
            grad_La = gradbary_vector[0];
            grad_Lb = gradbary_vector[1];
            grad_Lc = gradbary_vector[2];
            grad_Ld = gradbary_vector[3];
            grad_Le = gradbary_vector[4];
            
            a = 0;
            b = 1;
            c = 2;
            d = 3;
            e = 4;
            
        }
        else if (r==1)
        {
            //case where (a,b,c,d,e) = (1,2,3,4,0)
            // define lamda
            La = bary_vector[1];
            Lb = bary_vector[2];
            Lc = bary_vector[3];
            Ld = bary_vector[4];
            Le = bary_vector[0];
            
            // define grad(lamda)
            grad_La = gradbary_vector[1];
            grad_Lb = gradbary_vector[2];
            grad_Lc = gradbary_vector[3];
            grad_Ld = gradbary_vector[4];
            grad_Le = gradbary_vector[0];
            
            a = 1;
            b = 2;
            c = 3;
            d = 4;
            e = 0;
            
        }
        else if (r==2)
        {
            //case where (a,b,c,d,e) = (2,3,4,0,1)
            // define lamda
            La = bary_vector[2];
            Lb = bary_vector[3];
            Lc = bary_vector[4];
            Ld = bary_vector[0];
            Le = bary_vector[1];
            
            // define grad(lamda)
            grad_La = gradbary_vector[2];
            grad_Lb = gradbary_vector[3];
            grad_Lc = gradbary_vector[4];
            grad_Ld = gradbary_vector[0];
            grad_Le = gradbary_vector[1];
            
            a = 2;
            b = 3;
            c = 4;
            d = 0;
            e = 1;
            
        }
        else if (r==3)
        {
            //case where (a,b,c,d,e) = (3,4,0,1,2)
            // define lamda
            La = bary_vector[3];
            Lb = bary_vector[4];
            Lc = bary_vector[0];
            Ld = bary_vector[1];
            Le = bary_vector[2];
            
            // define grad(lamda)
            grad_La = gradbary_vector[3];
            grad_Lb = gradbary_vector[4];
            grad_Lc = gradbary_vector[0];
            grad_Ld = gradbary_vector[1];
            grad_Le = gradbary_vector[2];
            
            a = 3;
            b = 4;
            c = 0;
            d = 1;
            e = 2;
            
        }
        else
        {
            mfem_error("Invaild Bubble");
        }
        
        for (int i=0; i<=p; i++)
        {
            for (int j=0; j<=p; j++)
            {
                for (int l=0; l<=p; l++)
                {
                    for (int m=1; m<=p; m++)
                    {
                        if ((i+j+l+m)<=p)
                        {
                            // compute polynomials
                            
                            // Legendre Poly
                            std::vector<double> Legendre_i;
                            std::vector<double> Legendre_i_ref;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            // Derivative Legendre Poly
                            std::vector<double> Legendre_i_dx;
                            std::vector<double> Legendre_i_dt;
                            poly1d.CalcScaledLegendreDerivative(i, x, y, Legendre_i_ref, Legendre_i_dx, Legendre_i_dt);
                            
                            // Jacobi Poly
                            double alpha = 2.0*i+1.0;
                            x = Lc;
                            y = La + Lb + Lc;
                            std::vector<double> Jacobi_j;
                            std::vector<double> Jacobi_j_ref;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            // Derivative Jacobi Poly
                            std::vector<double> Jacobi_j_dx;
                            std::vector<double> Jacobi_j_dt;
                            poly1d.CalcScaledJacobiDerivative(j, alpha, x, y, Jacobi_j_ref, Jacobi_j_dx, Jacobi_j_dt);
                            
                            // Jacobi Poly
                            alpha = 2.0*(i+j+1.0);
                            x = Ld;
                            y = La + Lb + Lc + Ld;
                            std::vector<double> Jacobi_l;
                            std::vector<double> Jacobi_l_ref;
                            poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                            // Derivative Jacobi Poly
                            std::vector<double> Jacobi_l_dx;
                            std::vector<double> Jacobi_l_dt;
                            poly1d.CalcScaledJacobiDerivative(l, alpha, x, y, Jacobi_l_ref, Jacobi_l_dx, Jacobi_l_dt);
                            
                            //  Int Jacobi Poly
                            alpha = 2.0*(i+j+l)+3.0;
                            x = Le;
                            y = 1.0;
                            std::vector<double> Int_Jacobi_m;
                            poly1d.CalcIntJacobi(m, x, y, alpha, Int_Jacobi_m);
                            // Derivative Int Jacobi Poly
                            std::vector<double> Int_Jacobi_m_dx; // Use i-1 polynomial
                            x = Le;
                            y = 1.0;
                            alpha = 2.0*(i+j+l)+3.0; // added + 3.0
                            poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m_dx);
                            
                            //Compute relevant four curl
                            std::vector<double> curl_bcd = four_cross(b,c,d);
                            std::vector<double> curl_cda = four_cross(c,d,a);
                            std::vector<double> curl_dab = four_cross(d,a,b);
                            std::vector<double> curl_abc = four_cross(a,b,c);
                            
                            
                            // Basis Functions can be decomposed into vector and scalar parts

                            // Scalar Part
                            double Scalar = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1];
                            
                            // d scalar / dx
                            double dscalar_x = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[0] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[0] + grad_Lb[0]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[0] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[0] + grad_Lb[0] + grad_Lc[0]))*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[0] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[0] + grad_Lb[0] + grad_Lc[0] + grad_Ld[0]))*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*(Int_Jacobi_m_dx[Int_Jacobi_m_dx.size()-2]*grad_Le[0]);
                                                        
                            // d scalar / dy
                            double dscalar_y = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[1] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[1] + grad_Lb[1]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[1] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[1] + grad_Lb[1] + grad_Lc[1]))*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[1] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[1] + grad_Lb[1] + grad_Lc[1] + grad_Ld[1]))*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*(Int_Jacobi_m_dx[Int_Jacobi_m_dx.size()-2]*grad_Le[1]);
                            
                            
                            // d scalar / dz
                            double dscalar_z = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[2] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[2] + grad_Lb[2]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[2] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[2] + grad_Lb[2] + grad_Lc[2]))*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[2] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[2] + grad_Lb[2] + grad_Lc[2] + grad_Ld[2]))*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*(Int_Jacobi_m_dx[Int_Jacobi_m_dx.size()-2]*grad_Le[2]);
                            
                            
                            // d scalar / dt
                            double dscalar_t = (Legendre_i_dx[Legendre_i_dx.size()-1]*grad_Lb[3] + Legendre_i_dt[Legendre_i_dt.size()-1]*(grad_La[3] + grad_Lb[3]))*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*(Jacobi_j_dx[Jacobi_j_dx.size()-1]*grad_Lc[3] + Jacobi_j_dt[Jacobi_j_dt.size()-1]*(grad_La[3] + grad_Lb[3] + grad_Lc[3]))*Jacobi_l[Jacobi_l.size()-1]*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*(Jacobi_l_dx[Jacobi_l_dx.size()-1]*grad_Ld[3] + Jacobi_l_dt[Jacobi_l_dt.size()-1]*(grad_La[3] + grad_Lb[3] + grad_Lc[3] + grad_Ld[3]))*Int_Jacobi_m[Int_Jacobi_m.size()-1] +
                            
                                Legendre_i[Legendre_i.size()-1]*Jacobi_j[Jacobi_j.size()-1]*Jacobi_l[Jacobi_l.size()-1]*(Int_Jacobi_m_dx[Int_Jacobi_m_dx.size()-2]*grad_Le[3]);
                            
                            // components of vectors part
                            double vec_val_x = La*curl_bcd[0] - Lb*curl_cda[0] + Lc*curl_dab[0] - Ld*curl_abc[0];
                            
                            double vec_val_y = La*curl_bcd[1] - Lb*curl_cda[1] + Lc*curl_dab[1] - Ld*curl_abc[1];

                            double vec_val_z = La*curl_bcd[2] - Lb*curl_cda[2] + Lc*curl_dab[2] - Ld*curl_abc[2];

                            double vec_val_t = La*curl_bcd[3] - Lb*curl_cda[3] + Lc*curl_dab[3] - Ld*curl_abc[3];
                            
                            // Each term of div(vector)
                            double div_val_x = grad_La[0]*curl_bcd[0] - grad_Lb[0]*curl_cda[0] + grad_Lc[0]*curl_dab[0] - grad_Ld[0]*curl_abc[0];
                            
                            double div_val_y = grad_La[1]*curl_bcd[1] - grad_Lb[1]*curl_cda[1] + grad_Lc[1]*curl_dab[1] - grad_Ld[1]*curl_abc[1];

                            double div_val_z = grad_La[2]*curl_bcd[2] - grad_Lb[2]*curl_cda[2] + grad_Lc[2]*curl_dab[2] - grad_Ld[2]*curl_abc[2];

                            double div_val_t = grad_La[3]*curl_bcd[3] - grad_Lb[3]*curl_cda[3] + grad_Lc[3]*curl_dab[3] - grad_Ld[3]*curl_abc[3];
                            
                            // Div value
                            double div_vec = div_val_x + div_val_y + div_val_z + div_val_t;
                            
                            
                            // Add Divergence of Funtion
                            divu(o++) = (dscalar_x * vec_val_x + dscalar_y * vec_val_y + dscalar_z * vec_val_z + dscalar_t * vec_val_t) + (div_vec * Scalar);
                            
                            

                        }
                    }
                }
            }
        }
    }// End of Interiors
    

   Ti.Mult(divu, divshape);
    //std::cout << "Made it divshape" << std::endl;
}

void Hdiv_PentatopeElement::ProjectDivSkew(const FiniteElement& fe,
                                         ElementTransformation& Trans, DenseMatrix& DivSkew)
{
    mfem_error("HDiv_PentatopeElement_Fuentes::ProjectDivSkew not implemented");

}

//void Hdiv_PentatopeElement::GetLocalInterpolation(ElementTransformation &Trans,
//                                   DenseMatrix &I) const
//{
//    mfem_error("HDiv_PentatopeElement_Fuentes::GetLocalInterpolation not implemented");
//
//}
//
//void Hdiv_PentatopeElement::GetLocalRestriction(ElementTransformation &Trans,
//                                 DenseMatrix &R) const
//{
//    mfem_error("HDiv_PentatopeElement_Fuentes::GetLocalRestriction not implemented");
//
//}
//
//void Hdiv_PentatopeElement::GetTransferMatrix(const FiniteElement &fe,
//                               ElementTransformation &Trans,
//                               DenseMatrix &I) const
//{
//    mfem_error("HDiv_PentatopeElement_Fuentes::GetTransferMatrix not implemented");
//
//}
//
//void Hdiv_PentatopeElement::Project(VectorCoefficient &vc,
//                     ElementTransformation &Trans, Vector &dofs) const
//{
//    //mfem_error("HDiv_PentatopeElement_Fuentes::Project not implemented");
//    Project_RT(nk, dof2nk, vc, Trans, dofs);
//
//}
//
//void Hdiv_PentatopeElement::ProjectMatrixCoefficient(
//                                      MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const
//{
//    mfem_error("HDiv_PentatopeElement_Fuentes::ProjectMatrixCoefficient not implemented");
//
//}
//
//void Hdiv_PentatopeElement::Project(const FiniteElement &fe, ElementTransformation &Trans,
//                     DenseMatrix &I) const
//{
//    //mfem_error("HDiv_PentatopeElement_Fuentes::Project not implemented");
//    Project_RT(nk, dof2nk, fe, Trans, I);
//
//}

Hdiv_PentatopeElement_Barycentric::Hdiv_PentatopeElement_Barycentric(const int p, const int type)
   : NodalFiniteElement(4, Geometry::PENTATOPE,
                        ((p + 1)*(p + 2)*(p + 3)*(p + 4))/24,
                        p, FunctionSpace::Pk)
{
   const double *cp = poly1d.ClosedPoints(p, VerifyClosed(type));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   ddshape_x.SetSize(p + 1);
   ddshape_y.SetSize(p + 1);
   ddshape_z.SetSize(p + 1);
   ddshape_t.SetSize(p + 1);
   ddshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
   ddu.SetSize(dof,dim*(dim+1)/2 );
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
#endif

    int o = 0;

   // facets (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,2,3)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[0]);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,2,1,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[j]/w, cp[i]/w, cp[0], cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,1,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[i]/w, cp[0], cp[j]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (0,3,2,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[0], cp[j]/w, cp[i]/w, cp[k]/w);
         }
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)  // (1,2,3,4)
         {
            double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
            Nodes.IntPoint(o++).Set4(cp[p-i-j-k]/w, cp[i]/w, cp[j]/w, cp[k]/w);
         }

   // interior bubbles
   for (int l = 1; l < p; l++)
      for (int k = 1; k + l < p; k++)
         for (int j = 1; j + k + l < p; j++)
            for (int i = 1; i + j + k + l < p; i++)
            {
               double w = cp[i] + cp[j] + cp[k] + cp[l] + cp[p-i-j-k-l];
               Nodes.IntPoint(o++).Set4(cp[i]/w, cp[j]/w, cp[k]/w, cp[l]/w);
            }

   DenseMatrix T(dof);
    // Add code to assemble matrix
   Ti.Factor(T);
}

void Hdiv_PentatopeElement_Barycentric::CalcShape(const IntegrationPoint &ip,
                                                Vector &shape) const
{
const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
   Vector u(Dof);
#endif
    
    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);

    std::vector<double> bary_vector{L[0], L[1], L[2], L[3], 1-L[0]-L[1]-L[2]-L[3]};
    
    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};
    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    
    // Lamda Function for computing cross product with gradients of barycentric coordinates
    
    auto four_cross = [&](int A,int B,int C) -> std::vector<double>{
        
        // set grad of Lamda
        grad_La = gradbary_vector[A];
        grad_Lb = gradbary_vector[B];
        grad_Lc = gradbary_vector[C];
                
        // set variables
        double a,b,c,d,e,f,g,h,m,n,p,q;
        a = grad_La[0];
        b = grad_La[1];
        c = grad_La[2];
        d = grad_La[3];
        
        e = grad_Lb[0];
        f = grad_Lb[1];
        g = grad_Lb[2];
        h = grad_Lb[3];
        
        m = grad_Lc[0];
        n = grad_Lc[1];
        q = grad_Lc[2];
        p = grad_Lc[3];
        
        double i_hat, j_hat, k_hat, l_hat;
        
        i_hat = b*((g*p)-(q*h)) - c*((f*p)-(n*h)) + d*((f*q)-(g*n));
        
        j_hat = a*((g*p)-(q*h)) - c*((e*p)-(m*h)) + d*((e*q)-(g*m));
        
        k_hat = a*((f*p)-(n*h)) - b*((e*p)-(m*h)) + d*((e*n)-(f*m));
        
        l_hat = a*((f*q)-(n*g)) - b*((e*q)-(m*g)) + c*((e*n)-(m*f));
        
        std::vector<double> four_curl = {i_hat, j_hat, k_hat, l_hat};
        return four_curl;

    }; // end of lamda function
    
           
    //Facets
    for(int i=0; i<=p-1;i++)
    {
        for(int j=0; j<=p-1;j++)
        {
            for(int l=1; 1<=p-1;l++)
            {
                if((i+j+l)<=p-1)
                {
                    
                    // set Lamda
                    La = bary_vector[0];
                    Lb = bary_vector[1];
                    Lc = bary_vector[2];
                    Ld = bary_vector[3];
                    
                    // set grad of Lamda
                    grad_La = gradbary_vector[0];
                    grad_Lb = gradbary_vector[1];
                    grad_Lc = gradbary_vector[2];
                    grad_Ld = gradbary_vector[3];
                    
                    // compute polynomials
                    
                    std::vector<double> Legendre_i;
                    double x = Lb;
                    double y = La + Lb;
                    poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                    
                    double alpha = 2*i+1;
                    x = Lc;
                    y = La + Lb + Lc;
                    std::vector<double> Jacobi_j;
                    poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                    
                    alpha = 2*(i+j+1);
                    x = Ld;
                    y = La + Lb + Lc + Ld;
                    std::vector<double> Jacobi_l;
                    poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                    
                    //Compute relevant four curl
                    std::vector<double> curl_bcd = four_cross(1,2,3);
                    std::vector<double> curl_cda = four_cross(2,3,0);
                    std::vector<double> curl_dab = four_cross(1,0,3);
                    std::vector<double> curl_abc = four_cross(0,1,2);

                    
                    for(int index=0; index<4; index++)
                    {
                        double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                        // Compute Facet Funtion
                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * scalar_val;
                    }
                }
            }
        }
    } //end of Facets
           
    /*//Interiors
    for (int i=0; i<=p-1; i++)
    {
        for (int j=1; j<=p-1; j++)
        {
            for (int l=1; l<=p-1; l++)
            {
                for (int m=1; m<=p-1; m++)
                {
                    if ((i+j+l+m)<=p-1)
                    {
                        for (int r=1; r<5; r++)
                        {
                            if (r==1)
                            {
                                //case where (a,b,c,d,e) = (0,1,2,3,4)
                                
                                // define lamda
                                La = bary_vector[0];
                                Lb = bary_vector[1];
                                Lc = bary_vector[2];
                                Ld = bary_vector[3];
                                Le = bary_vector[4];
                                
                                // define grad(lamda)
                                grad_La = gradbary_vector[0];
                                grad_Lb = gradbary_vector[1];
                                grad_Lc = gradbary_vector[2];
                                grad_Ld = gradbary_vector[3];
                                grad_Le = gradbary_vector[4];

                                
                                // compute polynomials
                                
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                std::vector<double> Jacobi_j;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                std::vector<double> Jacobi_l;
                                poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                                
                                alpha = 2*(i+j+l)+3;
                                x = Le;
                                y = 1;
                                std::vector<double> Int_Jacobi_m;
                                poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                //Compute relevant four curl
                                std::vector<double> curl_bcd = four_cross(1,2,3);
                                std::vector<double> curl_cda = four_cross(2,3,0);
                                std::vector<double> curl_dab = four_cross(1,0,3);
                                std::vector<double> curl_abc = four_cross(0,1,2);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                                    // Compute Facet Funtion
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                            }
                            
                            if (r==2)
                            {
                                //case where (a,b,c,d,e) = (1,2,3,4,0)
                                
                                // define lamda
                                La = bary_vector[1];
                                Lb = bary_vector[2];
                                Lc = bary_vector[3];
                                Ld = bary_vector[4];
                                Le = bary_vector[0];
                                
                                // define grad(lamda)
                                grad_La = gradbary_vector[1];
                                grad_Lb = gradbary_vector[2];
                                grad_Lc = gradbary_vector[3];
                                grad_Ld = gradbary_vector[4];
                                grad_Le = gradbary_vector[0];
                                
                                // compute polynomials
                                
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                std::vector<double> Jacobi_j;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                std::vector<double> Jacobi_l;
                                poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                                
                                alpha = 2*(i+j+l)+3;
                                x = Le;
                                y = 1;
                                std::vector<double> Int_Jacobi_m;
                                poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                //Compute relevant four curl
                                std::vector<double> curl_bcd = four_curl(2,3,4);
                                std::vector<double> curl_cda = four_curl(3,4,1);
                                std::vector<double> curl_dab = four_curl(2,1,4);
                                std::vector<double> curl_abc = four_curl(1,2,3);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                                    // Compute Facet Funtion
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                            }
                            
                            if (r==3)
                            {
                                //case where (a,b,c,d,e) = (2,3,4,0,1)
                                
                                // define lamda
                                La = bary_vector[2];
                                Lb = bary_vector[3];
                                Lc = bary_vector[4];
                                Ld = bary_vector[0];
                                Le = bary_vector[1];
                                
                                // define grad(lamda)
                                grad_La = gradbary_vector[2];
                                grad_Lb = gradbary_vector[3];
                                grad_Lc = gradbary_vector[4];
                                grad_Ld = gradbary_vector[0];
                                grad_Le = gradbary_vector[1];
                                
                                // compute polynomials
                                
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                std::vector<double> Jacobi_j;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                std::vector<double> Jacobi_l;
                                poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                                
                                alpha = 2*(i+j+l)+3;
                                x = Le;
                                y = 1;
                                std::vector<double> Int_Jacobi_m;
                                poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                //Compute relevant four curl
                                std::vector<double> curl_bcd = four_curl(3,4,0);
                                std::vector<double> curl_cda = four_curl(4,0,2);
                                std::vector<double> curl_dab = four_curl(3,2,0);
                                std::vector<double> curl_abc = four_curl(2,3,4);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                                    // Compute Facet Funtion
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                            }
                            
                            if (r==4)
                            {
                                //case where (a,b,c,d,e) = (3,4,0,1,2)
                                
                                // define lamda
                                La = bary_vector[3];
                                Lb = bary_vector[4];
                                Lc = bary_vector[0];
                                Ld = bary_vector[1];
                                Le = bary_vector[2];
                                
                                // define grad(lamda)
                                grad_La = gradbary_vector[3];
                                grad_Lb = gradbary_vector[4];
                                grad_Lc = gradbary_vector[0];
                                grad_Ld = gradbary_vector[1];
                                grad_Le = gradbary_vector[2];



                                
                                // compute polynomials
                                
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                double alpha = 2*i+1;
                                x = Lc;
                                y = La + Lb + Lc;
                                std::vector<double> Jacobi_j;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                alpha = 2*(i+j+1);
                                x = Ld;
                                y = La + Lb + Lc + Ld;
                                std::vector<double> Jacobi_l;
                                poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                                
                                alpha = 2*(i+j+l)+3;
                                x = Le;
                                y = 1;
                                std::vector<double> Int_Jacobi_m;
                                poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m);
                                
                                //Compute relevant four curl
                                std::vector<double> curl_bcd = four_curl(4,0,1);
                                std::vector<double> curl_cda = four_curl(0,1,3);
                                std::vector<double> curl_dab = four_curl(4,3,1);
                                std::vector<double> curl_abc = four_curl(3,4,0);
                                
                                
                                for(int index=0; index<4; index++)
                                {
                                    double scalar_val = La*curl_bcd[index] - Lb*curl_cda[index] + Lc*curl_dab[index] - Ld*curl_abc[index];
                                    // Compute Facet Funtion
                                    u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1] * scalar_val;
                                }
                            }
                        }
                    }
                }
            }
        }
    } // end of interiors*/
   
    
    
//   poly1d.CalcBasis(p, ip.x, shape_x);
//   poly1d.CalcBasis(p, ip.y, shape_y);
//   poly1d.CalcBasis(p, ip.z, shape_z);
//   poly1d.CalcBasis(p, ip.t, shape_t);
//   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);
//
//   for (int o = 0, l = 0; l <= p; l++)
//      for (int k = 0; k + l <= p; k++)
//         for (int j = 0; j + k + l <= p; j++)
//            for (int i = 0; i + j + k + l <= p; i++)
//            {
//               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
//            }

   Ti.Mult(u, shape);
}

void Hdiv_PentatopeElement_Barycentric::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
    //Do Nothing
    mfem_error("Hdiv_PentatopeElement_Barycentric::CalcDShape not implemented");

}

void Hdiv_PentatopeElement_Barycentric::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
    //Do Nothing
    mfem_error("Hdiv_PentatopeElement_Barycentric::CalcHessian not implemented");
}

// No Constructor for L2 Yet!!

void L2_PentatopeElement_Barycentric::CalcShape(const IntegrationPoint &ip,
                                                Vector &shape) const
{
const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p+1),
          shape_l(p + 1);
   Vector u(Dof);
#endif
    
    double L1, L2, L3, L4, L5;
    //compute barycentric coordinates as function of ip
    double L[4];
    ip.Get(L,4);

    std::vector<double> bary_vector{L[0], L[1], L[2], L[3], 1-L[0]-L[1]-L[2]-L[3]};
    
    int o = 0;
    
    double La, Lb, Lc, Ld, Le;
    
    //Interiors
    for (int i=0; i<=p-1; i++)
    {
        for (int j=0; j<=p-1; j++)
        {
            for (int l=0; l<=p-1; l++)
            {
                for (int m=0; m<=p-1; m++)
                {
                    if ((i+j+l+m)<=p-1)
                    {
                        double La = bary_vector[0];
                        double Lb = bary_vector[1];
                        double Lc = bary_vector[2];
                        double Ld = bary_vector[3];
                        double Le = bary_vector[4];
                        
                        // compute polynomials
                        
                        std::vector<double> Legendre_i;
                        double x = Lb;
                        double y = La + Lb;
                        poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                        
                        double alpha = 2*i+1;
                        x = Lc;
                        y = La + Lb + Lc;
                        std::vector<double> Jacobi_j;
                        poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                        
                        alpha = 2*(i+j+1);
                        x = Ld;
                        y = La + Lb + Lc + Ld;
                        std::vector<double> Jacobi_l;
                        poly1d.CalcJacobi(l, x, y, alpha, Jacobi_l);
                        
                        alpha = 2*(i+j+l)+3;
                        x = Le;
                        y = 1;
                        std::vector<double> Int_Jacobi_m;
                        poly1d.CalcJacobi(m, x, y, alpha, Int_Jacobi_m);

                        u(o++) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * Jacobi_l[Jacobi_l.size()-1] * Int_Jacobi_m[Int_Jacobi_m.size()-1];
                    }
                }
            }
        }
    } // end of interiors
   
    
    
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
            }

   Ti.Mult(u, shape);
}

void L2_PentatopeElement_Barycentric::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
    //Do Nothing
    mfem_error("L2_PentatopeElement_Barycentric::CalcDShape not implemented");

}

void L2_PentatopeElement_Barycentric::CalcHessian(const IntegrationPoint &ip,
                                      DenseMatrix &ddshape) const
{
    //Do Nothing
    mfem_error("L2_PentatopeElement_Barycentric::CalcHessian not implemented");
}

H1Pos_TriangleElement::H1Pos_TriangleElement(const int p)
   : PositiveFiniteElement(2, Geometry::TRIANGLE, ((p + 1)*(p + 2))/2, p,
                           FunctionSpace::Pk)
{
#ifndef MFEM_THREAD_SAFE
   m_shape.SetSize(dof);
   dshape_1d.SetSize(p + 1);
   m_dshape.SetSize(dof, dim);
#endif
   dof_map.SetSize(dof);

   struct Index
   {
      int p2p3;
      Index(int p) { p2p3 = 2*p + 3; }
      int operator()(int i, int j) { return ((p2p3-j)*j)/2+i; }
   };
   Index idx(p);

   // vertices
   dof_map[idx(0,0)] = 0;
   Nodes.IntPoint(0).Set2(0., 0.);
   dof_map[idx(p,0)] = 1;
   Nodes.IntPoint(1).Set2(1., 0.);
   dof_map[idx(0,p)] = 2;
   Nodes.IntPoint(2).Set2(0., 1.);

   // edges
   int o = 3;
   for (int i = 1; i < p; i++)
   {
      dof_map[idx(i,0)] = o;
      Nodes.IntPoint(o++).Set2(double(i)/p, 0.);
   }
   for (int i = 1; i < p; i++)
   {
      dof_map[idx(p-i,i)] = o;
      Nodes.IntPoint(o++).Set2(double(p-i)/p, double(i)/p);
   }
   for (int i = 1; i < p; i++)
   {
      dof_map[idx(0,p-i)] = o;
      Nodes.IntPoint(o++).Set2(0., double(p-i)/p);
   }

   // interior
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)
      {
         dof_map[idx(i,j)] = o;
         Nodes.IntPoint(o++).Set2(double(i)/p, double(j)/p);
      }
}

// static method
void H1Pos_TriangleElement::CalcShape(
   const int p, const double l1, const double l2, double *shape)
{
   const double l3 = 1. - l1 - l2;

   // The (i,j) basis function is given by: T(i,j,p-i-j) l1^i l2^j l3^{p-i-j},
   // where T(i,j,k) = (i+j+k)! / (i! j! k!)
   // Another expression is given by the terms of the expansion:
   //    (l1 + l2 + l3)^p =
   //       \sum_{j=0}^p \binom{p}{j} l2^j
   //          \sum_{i=0}^{p-j} \binom{p-j}{i} l1^i l3^{p-j-i}
   const int *bp = Poly_1D::Binom(p);
   double z = 1.;
   for (int o = 0, j = 0; j <= p; j++)
   {
      Poly_1D::CalcBinomTerms(p - j, l1, l3, &shape[o]);
      double s = bp[j]*z;
      for (int i = 0; i <= p - j; i++)
      {
         shape[o++] *= s;
      }
      z *= l2;
   }
}

// static method
void H1Pos_TriangleElement::CalcDShape(
   const int p, const double l1, const double l2,
   double *dshape_1d, double *dshape)
{
   const int dof = ((p + 1)*(p + 2))/2;
   const double l3 = 1. - l1 - l2;

   const int *bp = Poly_1D::Binom(p);
   double z = 1.;
   for (int o = 0, j = 0; j <= p; j++)
   {
      Poly_1D::CalcDBinomTerms(p - j, l1, l3, dshape_1d);
      double s = bp[j]*z;
      for (int i = 0; i <= p - j; i++)
      {
         dshape[o++] = s*dshape_1d[i];
      }
      z *= l2;
   }
   z = 1.;
   for (int i = 0; i <= p; i++)
   {
      Poly_1D::CalcDBinomTerms(p - i, l2, l3, dshape_1d);
      double s = bp[i]*z;
      for (int o = i, j = 0; j <= p - i; j++)
      {
         dshape[dof + o] = s*dshape_1d[j];
         o += p + 1 - j;
      }
      z *= l1;
   }
}

void H1Pos_TriangleElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector m_shape(dof);
#endif
   CalcShape(order, ip.x, ip.y, m_shape.GetData());
   for (int i = 0; i < dof; i++)
   {
      shape(dof_map[i]) = m_shape(i);
   }
}

void H1Pos_TriangleElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector dshape_1d(order + 1);
   DenseMatrix m_dshape(dof, dim);
#endif
   CalcDShape(order, ip.x, ip.y, dshape_1d.GetData(), m_dshape.Data());
   for (int d = 0; d < 2; d++)
   {
      for (int i = 0; i < dof; i++)
      {
         dshape(dof_map[i],d) = m_dshape(i,d);
      }
   }
}


H1Pos_TetrahedronElement::H1Pos_TetrahedronElement(const int p)
   : PositiveFiniteElement(3, Geometry::TETRAHEDRON,
                           ((p + 1)*(p + 2)*(p + 3))/6, p, FunctionSpace::Pk)
{
#ifndef MFEM_THREAD_SAFE
   m_shape.SetSize(dof);
   dshape_1d.SetSize(p + 1);
   m_dshape.SetSize(dof, dim);
#endif
   dof_map.SetSize(dof);

   struct Index
   {
      int p, dof;
      int tri(int k) { return (k*(k + 1))/2; }
      int tet(int k) { return (k*(k + 1)*(k + 2))/6; }
      Index(int p_) { p = p_; dof = tet(p + 1); }
      int operator()(int i, int j, int k)
      { return dof - tet(p - k) - tri(p + 1 - k - j) + i; }
   };
   Index idx(p);

   // vertices
   dof_map[idx(0,0,0)] = 0;
   Nodes.IntPoint(0).Set3(0., 0., 0.);
   dof_map[idx(p,0,0)] = 1;
   Nodes.IntPoint(1).Set3(1., 0., 0.);
   dof_map[idx(0,p,0)] = 2;
   Nodes.IntPoint(2).Set3(0., 1., 0.);
   dof_map[idx(0,0,p)] = 3;
   Nodes.IntPoint(3).Set3(0., 0., 1.);

   // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
   int o = 4;
   for (int i = 1; i < p; i++)  // (0,1)
   {
      dof_map[idx(i,0,0)] = o;
      Nodes.IntPoint(o++).Set3(double(i)/p, 0., 0.);
   }
   for (int i = 1; i < p; i++)  // (0,2)
   {
      dof_map[idx(0,i,0)] = o;
      Nodes.IntPoint(o++).Set3(0., double(i)/p, 0.);
   }
   for (int i = 1; i < p; i++)  // (0,3)
   {
      dof_map[idx(0,0,i)] = o;
      Nodes.IntPoint(o++).Set3(0., 0., double(i)/p);
   }
   for (int i = 1; i < p; i++)  // (1,2)
   {
      dof_map[idx(p-i,i,0)] = o;
      Nodes.IntPoint(o++).Set3(double(p-i)/p, double(i)/p, 0.);
   }
   for (int i = 1; i < p; i++)  // (1,3)
   {
      dof_map[idx(p-i,0,i)] = o;
      Nodes.IntPoint(o++).Set3(double(p-i)/p, 0., double(i)/p);
   }
   for (int i = 1; i < p; i++)  // (2,3)
   {
      dof_map[idx(0,p-i,i)] = o;
      Nodes.IntPoint(o++).Set3(0., double(p-i)/p, double(i)/p);
   }

   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp)
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (1,2,3)
      {
         dof_map[idx(p-i-j,i,j)] = o;
         Nodes.IntPoint(o++).Set3(double(p-i-j)/p, double(i)/p, double(j)/p);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,3,2)
      {
         dof_map[idx(0,j,i)] = o;
         Nodes.IntPoint(o++).Set3(0., double(j)/p, double(i)/p);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,1,3)
      {
         dof_map[idx(i,0,j)] = o;
         Nodes.IntPoint(o++).Set3(double(i)/p, 0., double(j)/p);
      }
   for (int j = 1; j < p; j++)
      for (int i = 1; i + j < p; i++)  // (0,2,1)
      {
         dof_map[idx(j,i,0)] = o;
         Nodes.IntPoint(o++).Set3(double(j)/p, double(i)/p, 0.);
      }

   // interior
   for (int k = 1; k < p; k++)
      for (int j = 1; j + k < p; j++)
         for (int i = 1; i + j + k < p; i++)
         {
            dof_map[idx(i,j,k)] = o;
            Nodes.IntPoint(o++).Set3(double(i)/p, double(j)/p, double(k)/p);
         }
}

// static method
void H1Pos_TetrahedronElement::CalcShape(
   const int p, const double l1, const double l2, const double l3,
   double *shape)
{
   const double l4 = 1. - l1 - l2 - l3;

   // The basis functions are the terms in the expansion:
   //   (l1 + l2 + l3 + l4)^p =
   //      \sum_{k=0}^p \binom{p}{k} l3^k
   //         \sum_{j=0}^{p-k} \binom{p-k}{j} l2^j
   //            \sum_{i=0}^{p-k-j} \binom{p-k-j}{i} l1^i l4^{p-k-j-i}
   const int *bp = Poly_1D::Binom(p);
   double l3k = 1.;
   for (int o = 0, k = 0; k <= p; k++)
   {
      const int *bpk = Poly_1D::Binom(p - k);
      const double ek = bp[k]*l3k;
      double l2j = 1.;
      for (int j = 0; j <= p - k; j++)
      {
         Poly_1D::CalcBinomTerms(p - k - j, l1, l4, &shape[o]);
         double ekj = ek*bpk[j]*l2j;
         for (int i = 0; i <= p - k - j; i++)
         {
            shape[o++] *= ekj;
         }
         l2j *= l2;
      }
      l3k *= l3;
   }
}

// static method
void H1Pos_TetrahedronElement::CalcDShape(
   const int p, const double l1, const double l2, const double l3,
   double *dshape_1d, double *dshape)
{
   const int dof = ((p + 1)*(p + 2)*(p + 3))/6;
   const double l4 = 1. - l1 - l2 - l3;

   // For the x derivatives, differentiate the terms of the expression:
   //   \sum_{k=0}^p \binom{p}{k} l3^k
   //      \sum_{j=0}^{p-k} \binom{p-k}{j} l2^j
   //         \sum_{i=0}^{p-k-j} \binom{p-k-j}{i} l1^i l4^{p-k-j-i}
   const int *bp = Poly_1D::Binom(p);
   double l3k = 1.;
   for (int o = 0, k = 0; k <= p; k++)
   {
      const int *bpk = Poly_1D::Binom(p - k);
      const double ek = bp[k]*l3k;
      double l2j = 1.;
      for (int j = 0; j <= p - k; j++)
      {
         Poly_1D::CalcDBinomTerms(p - k - j, l1, l4, dshape_1d);
         double ekj = ek*bpk[j]*l2j;
         for (int i = 0; i <= p - k - j; i++)
         {
            dshape[o++] = dshape_1d[i]*ekj;
         }
         l2j *= l2;
      }
      l3k *= l3;
   }
   // For the y derivatives, differentiate the terms of the expression:
   //   \sum_{k=0}^p \binom{p}{k} l3^k
   //      \sum_{i=0}^{p-k} \binom{p-k}{i} l1^i
   //         \sum_{j=0}^{p-k-i} \binom{p-k-i}{j} l2^j l4^{p-k-j-i}
   l3k = 1.;
   for (int ok = 0, k = 0; k <= p; k++)
   {
      const int *bpk = Poly_1D::Binom(p - k);
      const double ek = bp[k]*l3k;
      double l1i = 1.;
      for (int i = 0; i <= p - k; i++)
      {
         Poly_1D::CalcDBinomTerms(p - k - i, l2, l4, dshape_1d);
         double eki = ek*bpk[i]*l1i;
         int o = ok + i;
         for (int j = 0; j <= p - k - i; j++)
         {
            dshape[dof + o] = dshape_1d[j]*eki;
            o += p - k - j + 1;
         }
         l1i *= l1;
      }
      l3k *= l3;
      ok += ((p - k + 2)*(p - k + 1))/2;
   }
   // For the z derivatives, differentiate the terms of the expression:
   //   \sum_{j=0}^p \binom{p}{j} l2^j
   //      \sum_{i=0}^{p-j} \binom{p-j}{i} l1^i
   //         \sum_{k=0}^{p-j-i} \binom{p-j-i}{k} l3^k l4^{p-k-j-i}
   double l2j = 1.;
   for (int j = 0; j <= p; j++)
   {
      const int *bpj = Poly_1D::Binom(p - j);
      const double ej = bp[j]*l2j;
      double l1i = 1.;
      for (int i = 0; i <= p - j; i++)
      {
         Poly_1D::CalcDBinomTerms(p - j - i, l3, l4, dshape_1d);
         double eji = ej*bpj[i]*l1i;
         int m = ((p + 2)*(p + 1))/2;
         int n = ((p - j + 2)*(p - j + 1))/2;
         for (int o = i, k = 0; k <= p - j - i; k++)
         {
            // m = ((p - k + 2)*(p - k + 1))/2;
            // n = ((p - k - j + 2)*(p - k - j + 1))/2;
            o += m;
            dshape[2*dof + o - n] = dshape_1d[k]*eji;
            m -= p - k + 1;
            n -= p - k - j + 1;
         }
         l1i *= l1;
      }
      l2j *= l2;
   }
}

void H1Pos_TetrahedronElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector m_shape(dof);
#endif
   CalcShape(order, ip.x, ip.y, ip.z, m_shape.GetData());
   for (int i = 0; i < dof; i++)
   {
      shape(dof_map[i]) = m_shape(i);
   }
}

void H1Pos_TetrahedronElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector dshape_1d(order + 1);
   DenseMatrix m_dshape(dof, dim);
#endif
   CalcDShape(order, ip.x, ip.y, ip.z, dshape_1d.GetData(), m_dshape.Data());
   for (int d = 0; d < 3; d++)
   {
      for (int i = 0; i < dof; i++)
      {
         dshape(dof_map[i],d) = m_dshape(i,d);
      }
   }
}


H1_WedgeElement::H1_WedgeElement(const int p,
                                 const int btype)
   : NodalFiniteElement(3, Geometry::PRISM, ((p + 1)*(p + 1)*(p + 2))/2,
                        p, FunctionSpace::Qk),
     TriangleFE(p, btype),
     SegmentFE(p, btype)
{
#ifndef MFEM_THREAD_SAFE
   t_shape.SetSize(TriangleFE.GetDof());
   s_shape.SetSize(SegmentFE.GetDof());
   t_dshape.SetSize(TriangleFE.GetDof(), 2);
   s_dshape.SetSize(SegmentFE.GetDof(), 1);
#endif

   t_dof.SetSize(dof);
   s_dof.SetSize(dof);

   int p2p3 = 2*p + 3, ntri = ((p + 1)*(p + 2))/2;
   auto idx = [p2p3,ntri](int i, int j, int k)
   {
      return k*ntri + ((p2p3-j)*j)/2+i;
   };

   lex_ordering.SetSize(dof);
   int o = 0;

   // Nodal DoFs
   lex_ordering[idx(0,0,0)] = o++;
   lex_ordering[idx(p,0,0)] = o++;
   lex_ordering[idx(0,p,0)] = o++;
   lex_ordering[idx(0,0,p)] = o++;
   lex_ordering[idx(p,0,p)] = o++;
   lex_ordering[idx(0,p,p)] = o++;
   t_dof[0] = 0; s_dof[0] = 0;
   t_dof[1] = 1; s_dof[1] = 0;
   t_dof[2] = 2; s_dof[2] = 0;
   t_dof[3] = 0; s_dof[3] = 1;
   t_dof[4] = 1; s_dof[4] = 1;
   t_dof[5] = 2; s_dof[5] = 1;

   // Edge DoFs
   int k = 0;
   int ne = p-1;
   for (int i=1; i<p; i++)
   {
      lex_ordering[idx(i,0,0)] = o + 0*ne + k;
      lex_ordering[idx(p-i,i,0)] = o + 1*ne + k;
      lex_ordering[idx(0,p-i,0)] = o + 2*ne + k;
      lex_ordering[idx(i,0,p)] = o + 3*ne + k;
      lex_ordering[idx(p-i,i,p)] = o + 4*ne + k;
      lex_ordering[idx(0,p-i,p)] = o + 5*ne + k;
      lex_ordering[idx(0,0,i)] = o + 6*ne + k;
      lex_ordering[idx(p,0,i)] = o + 7*ne + k;
      lex_ordering[idx(0,p,i)] = o + 8*ne + k;
      t_dof[5 + 0 * ne + i] = 2 + 0 * ne + i; s_dof[5 + 0 * ne + i] = 0;
      t_dof[5 + 1 * ne + i] = 2 + 1 * ne + i; s_dof[5 + 1 * ne + i] = 0;
      t_dof[5 + 2 * ne + i] = 2 + 2 * ne + i; s_dof[5 + 2 * ne + i] = 0;
      t_dof[5 + 3 * ne + i] = 2 + 0 * ne + i; s_dof[5 + 3 * ne + i] = 1;
      t_dof[5 + 4 * ne + i] = 2 + 1 * ne + i; s_dof[5 + 4 * ne + i] = 1;
      t_dof[5 + 5 * ne + i] = 2 + 2 * ne + i; s_dof[5 + 5 * ne + i] = 1;
      t_dof[5 + 6 * ne + i] = 0;              s_dof[5 + 6 * ne + i] = i + 1;
      t_dof[5 + 7 * ne + i] = 1;              s_dof[5 + 7 * ne + i] = i + 1;
      t_dof[5 + 8 * ne + i] = 2;              s_dof[5 + 8 * ne + i] = i + 1;
      ++k;
   }
   o += 9*ne;

   // Triangular Face DoFs
   k=0;
   int nt = (p-1)*(p-2)/2;
   for (int j=1; j<p; j++)
   {
      for (int i=1; i<p-j; i++)
      {
         int l = j - p + (((2 * p - 1) - i) * i) / 2;
         lex_ordering[idx(i,j,0)] = o+l;
         lex_ordering[idx(i,j,p)] = o+nt+k;
         t_dof[6 + 9 * ne + k]      = 3 * p + l; s_dof[6 + 9 * ne + k]      = 0;
         t_dof[6 + 9 * ne + nt + k] = 3 * p + k; s_dof[6 + 9 * ne + nt + k] = 1;
         k++;
      }
   }
   o += 2*nt;

   // Quadrilateral Face DoFs
   k=0;
   int nq = (p-1)*(p-1);
   for (int j=1; j<p; j++)
   {
      for (int i=1; i<p; i++)
      {
         lex_ordering[idx(i,0,j)] = o+k;
         lex_ordering[idx(p-i,i,j)] = o+nq+k;
         lex_ordering[idx(0,p-i,j)] = o+2*nq+k;

         t_dof[6 + 9 * ne + 2 * nt + 0 * nq + k] = 2 + 0 * ne + i;
         t_dof[6 + 9 * ne + 2 * nt + 1 * nq + k] = 2 + 1 * ne + i;
         t_dof[6 + 9 * ne + 2 * nt + 2 * nq + k] = 2 + 2 * ne + i;

         s_dof[6 + 9 * ne + 2 * nt + 0 * nq + k] = 1 + j;
         s_dof[6 + 9 * ne + 2 * nt + 1 * nq + k] = 1 + j;
         s_dof[6 + 9 * ne + 2 * nt + 2 * nq + k] = 1 + j;

         k++;
      }
   }
   o += 3*nq;

   // Interior DoFs
   int m=0;
   for (int k=1; k<p; k++)
   {
      int l=0;
      for (int j=1; j<p; j++)
      {
         for (int i=1; i+j<p; i++)
         {
            lex_ordering[idx(i,j,k)] = o++;
            t_dof[6 + 9 * ne + 2 * nt + 3 * nq + m] = 3 * p + l;
            s_dof[6 + 9 * ne + 2 * nt + 3 * nq + m] = 1 + k;
            l++; m++;
         }
      }
   }

   // Define Nodes
   const IntegrationRule & t_Nodes = TriangleFE.GetNodes();
   const IntegrationRule & s_Nodes = SegmentFE.GetNodes();
   for (int i=0; i<dof; i++)
   {
      Nodes.IntPoint(i).x = t_Nodes.IntPoint(t_dof[i]).x;
      Nodes.IntPoint(i).y = t_Nodes.IntPoint(t_dof[i]).y;
      Nodes.IntPoint(i).z = s_Nodes.IntPoint(s_dof[i]).x;
   }
}

void H1_WedgeElement::CalcShape(const IntegrationPoint &ip,
                                Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector t_shape(TriangleFE.GetDof());
   Vector s_shape(SegmentFE.GetDof());
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   SegmentFE.CalcShape(ipz, s_shape);

   for (int i=0; i<dof; i++)
   {
      shape[i] = t_shape[t_dof[i]] * s_shape[s_dof[i]];
   }
}

void H1_WedgeElement::CalcDShape(const IntegrationPoint &ip,
                                 DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector      t_shape(TriangleFE.GetDof());
   DenseMatrix t_dshape(TriangleFE.GetDof(), 2);
   Vector      s_shape(SegmentFE.GetDof());
   DenseMatrix s_dshape(SegmentFE.GetDof(), 1);
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   TriangleFE.CalcDShape(ip, t_dshape);
   SegmentFE.CalcShape(ipz, s_shape);
   SegmentFE.CalcDShape(ipz, s_dshape);

   for (int i=0; i<dof; i++)
   {
      dshape(i, 0) = t_dshape(t_dof[i],0) * s_shape[s_dof[i]];
      dshape(i, 1) = t_dshape(t_dof[i],1) * s_shape[s_dof[i]];
      dshape(i, 2) = t_shape[t_dof[i]] * s_dshape(s_dof[i],0);
   }
}


H1Pos_WedgeElement::H1Pos_WedgeElement(const int p)
   : PositiveFiniteElement(3, Geometry::PRISM,
                           ((p + 1)*(p + 1)*(p + 2))/2, p, FunctionSpace::Qk),
     TriangleFE(p),
     SegmentFE(p)
{
#ifndef MFEM_THREAD_SAFE
   t_shape.SetSize(TriangleFE.GetDof());
   s_shape.SetSize(SegmentFE.GetDof());
   t_dshape.SetSize(TriangleFE.GetDof(), 2);
   s_dshape.SetSize(SegmentFE.GetDof(), 1);
#endif

   t_dof.SetSize(dof);
   s_dof.SetSize(dof);

   // Nodal DoFs
   t_dof[0] = 0; s_dof[0] = 0;
   t_dof[1] = 1; s_dof[1] = 0;
   t_dof[2] = 2; s_dof[2] = 0;
   t_dof[3] = 0; s_dof[3] = 1;
   t_dof[4] = 1; s_dof[4] = 1;
   t_dof[5] = 2; s_dof[5] = 1;

   // Edge DoFs
   int ne = p-1;
   for (int i=1; i<p; i++)
   {
      t_dof[5 + 0 * ne + i] = 2 + 0 * ne + i; s_dof[5 + 0 * ne + i] = 0;
      t_dof[5 + 1 * ne + i] = 2 + 1 * ne + i; s_dof[5 + 1 * ne + i] = 0;
      t_dof[5 + 2 * ne + i] = 2 + 2 * ne + i; s_dof[5 + 2 * ne + i] = 0;
      t_dof[5 + 3 * ne + i] = 2 + 0 * ne + i; s_dof[5 + 3 * ne + i] = 1;
      t_dof[5 + 4 * ne + i] = 2 + 1 * ne + i; s_dof[5 + 4 * ne + i] = 1;
      t_dof[5 + 5 * ne + i] = 2 + 2 * ne + i; s_dof[5 + 5 * ne + i] = 1;
      t_dof[5 + 6 * ne + i] = 0;              s_dof[5 + 6 * ne + i] = i + 1;
      t_dof[5 + 7 * ne + i] = 1;              s_dof[5 + 7 * ne + i] = i + 1;
      t_dof[5 + 8 * ne + i] = 2;              s_dof[5 + 8 * ne + i] = i + 1;
   }

   // Triangular Face DoFs
   int k=0;
   int nt = (p-1)*(p-2)/2;
   for (int j=1; j<p; j++)
   {
      for (int i=1; i<j; i++)
      {
         t_dof[6 + 9 * ne + k]      = 3 * p + k; s_dof[6 + 9 * ne + k]      = 0;
         t_dof[6 + 9 * ne + nt + k] = 3 * p + k; s_dof[6 + 9 * ne + nt + k] = 1;
         k++;
      }
   }

   // Quadrilateral Face DoFs
   k=0;
   int nq = (p-1)*(p-1);
   for (int j=1; j<p; j++)
   {
      for (int i=1; i<p; i++)
      {
         t_dof[6 + 9 * ne + 2 * nt + 0 * nq + k] = 2 + 0 * ne + i;
         t_dof[6 + 9 * ne + 2 * nt + 1 * nq + k] = 2 + 1 * ne + i;
         t_dof[6 + 9 * ne + 2 * nt + 2 * nq + k] = 2 + 2 * ne + i;

         s_dof[6 + 9 * ne + 2 * nt + 0 * nq + k] = 1 + j;
         s_dof[6 + 9 * ne + 2 * nt + 1 * nq + k] = 1 + j;
         s_dof[6 + 9 * ne + 2 * nt + 2 * nq + k] = 1 + j;

         k++;
      }
   }

   // Interior DoFs
   int m=0;
   for (int k=1; k<p; k++)
   {
      int l=0;
      for (int j=1; j<p; j++)
      {
         for (int i=1; i<j; i++)
         {
            t_dof[6 + 9 * ne + 2 * nt + 3 * nq + m] = 3 * p + l;
            s_dof[6 + 9 * ne + 2 * nt + 3 * nq + m] = 1 + k;
            l++; m++;
         }
      }
   }

   // Define Nodes
   const IntegrationRule & t_Nodes = TriangleFE.GetNodes();
   const IntegrationRule & s_Nodes = SegmentFE.GetNodes();
   for (int i=0; i<dof; i++)
   {
      Nodes.IntPoint(i).x = t_Nodes.IntPoint(t_dof[i]).x;
      Nodes.IntPoint(i).y = t_Nodes.IntPoint(t_dof[i]).y;
      Nodes.IntPoint(i).z = s_Nodes.IntPoint(s_dof[i]).x;
   }
}

void H1Pos_WedgeElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector t_shape(TriangleFE.GetDof());
   Vector s_shape(SegmentFE.GetDof());
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   SegmentFE.CalcShape(ipz, s_shape);

   for (int i=0; i<dof; i++)
   {
      shape[i] = t_shape[t_dof[i]] * s_shape[s_dof[i]];
   }
}

void H1Pos_WedgeElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector      t_shape(TriangleFE.GetDof());
   DenseMatrix t_dshape(TriangleFE.GetDof(), 2);
   Vector      s_shape(SegmentFE.GetDof());
   DenseMatrix s_dshape(SegmentFE.GetDof(), 1);
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   TriangleFE.CalcDShape(ip, t_dshape);
   SegmentFE.CalcShape(ipz, s_shape);
   SegmentFE.CalcDShape(ipz, s_dshape);

   for (int i=0; i<dof; i++)
   {
      dshape(i, 0) = t_dshape(t_dof[i],0) * s_shape[s_dof[i]];
      dshape(i, 1) = t_dshape(t_dof[i],1) * s_shape[s_dof[i]];
      dshape(i, 2) = t_shape[t_dof[i]] * s_dshape(s_dof[i],0);
   }
}


L2_SegmentElement::L2_SegmentElement(const int p, const int btype)
   : NodalTensorFiniteElement(1, p, VerifyOpen(btype), L2_DOF_MAP)
{
   const double *op = poly1d.OpenPoints(p, btype);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   dshape_x.SetDataAndSize(NULL, p + 1);
#endif

   for (int i = 0; i <= p; i++)
   {
      Nodes.IntPoint(i).x = op[i];
   }
}

void L2_SegmentElement::CalcShape(const IntegrationPoint &ip,
                                  Vector &shape) const
{
   basis1d.Eval(ip.x, shape);
}

void L2_SegmentElement::CalcDShape(const IntegrationPoint &ip,
                                   DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(dof), dshape_x(dshape.Data(), dof);
#else
   dshape_x.SetData(dshape.Data());
#endif
   basis1d.Eval(ip.x, shape_x, dshape_x);
}

void L2_SegmentElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *op = poly1d.OpenPoints(p, b_type);

   switch (vertex)
   {
      case 0:
         for (int i = 0; i <= p; i++)
         {
            dofs(i) = poly1d.CalcDelta(p,(1.0 - op[i]));
         }
         break;

      case 1:
         for (int i = 0; i <= p; i++)
         {
            dofs(i) = poly1d.CalcDelta(p,op[i]);
         }
         break;
   }
}


L2Pos_SegmentElement::L2Pos_SegmentElement(const int p)
   : PositiveTensorFiniteElement(1, p, L2_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   dshape_x.SetDataAndSize(NULL, p + 1);
#endif

   if (p == 0)
   {
      Nodes.IntPoint(0).x = 0.5;
   }
   else
   {
      for (int i = 0; i <= p; i++)
      {
         Nodes.IntPoint(i).x = double(i)/p;
      }
   }
}

void L2Pos_SegmentElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   Poly_1D::CalcBernstein(order, ip.x, shape);
}

void L2Pos_SegmentElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(dof), dshape_x(dshape.Data(), dof);
#else
   dshape_x.SetData(dshape.Data());
#endif
   Poly_1D::CalcBernstein(order, ip.x, shape_x, dshape_x);
}

void L2Pos_SegmentElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   dofs[vertex*order] = 1.0;
}


L2_QuadrilateralElement::L2_QuadrilateralElement(const int p, const int btype)
   : NodalTensorFiniteElement(2, p, VerifyOpen(btype), L2_DOF_MAP)
{
   const double *op = poly1d.OpenPoints(p, b_type);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
#endif

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         Nodes.IntPoint(o++).Set2(op[i], op[j]);
      }
}

void L2_QuadrilateralElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   basis1d.Eval(ip.x, shape_x);
   basis1d.Eval(ip.y, shape_y);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         shape(o++) = shape_x(i)*shape_y(j);
      }
}

void L2_QuadrilateralElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), dshape_x(p+1), dshape_y(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         dshape(o,0) = dshape_x(i)* shape_y(j);
         dshape(o,1) =  shape_x(i)*dshape_y(j);  o++;
      }
}

void L2_QuadrilateralElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *op = poly1d.OpenPoints(p, b_type);

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   for (int i = 0; i <= p; i++)
   {
      shape_x(i) = poly1d.CalcDelta(p,(1.0 - op[i]));
      shape_y(i) = poly1d.CalcDelta(p,op[i]);
   }

   switch (vertex)
   {
      case 0:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs[o++] = shape_x(i)*shape_x(j);
            }
         break;
      case 1:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs[o++] = shape_y(i)*shape_x(j);
            }
         break;
      case 2:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs[o++] = shape_y(i)*shape_y(j);
            }
         break;
      case 3:
         for (int o = 0, j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               dofs[o++] = shape_x(i)*shape_y(j);
            }
         break;
   }
}


L2Pos_QuadrilateralElement::L2Pos_QuadrilateralElement(const int p)
   : PositiveTensorFiniteElement(2, p, L2_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
#endif

   if (p == 0)
   {
      Nodes.IntPoint(0).Set2(0.5, 0.5);
   }
   else
   {
      for (int o = 0, j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            Nodes.IntPoint(o++).Set2(double(i)/p, double(j)/p);
         }
   }
}

void L2Pos_QuadrilateralElement::CalcShape(const IntegrationPoint &ip,
                                           Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x);
   Poly_1D::CalcBernstein(p, ip.y, shape_y);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         shape(o++) = shape_x(i)*shape_y(j);
      }
}

void L2Pos_QuadrilateralElement::CalcDShape(const IntegrationPoint &ip,
                                            DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), dshape_x(p+1), dshape_y(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x, dshape_x);
   Poly_1D::CalcBernstein(p, ip.y, shape_y, dshape_y);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i <= p; i++)
      {
         dshape(o,0) = dshape_x(i)* shape_y(j);
         dshape(o,1) =  shape_x(i)*dshape_y(j);  o++;
      }
}

void L2Pos_QuadrilateralElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;

   dofs = 0.0;
   switch (vertex)
   {
      case 0: dofs[0] = 1.0; break;
      case 1: dofs[p] = 1.0; break;
      case 2: dofs[p*(p + 2)] = 1.0; break;
      case 3: dofs[p*(p + 1)] = 1.0; break;
   }
}


L2_HexahedronElement::L2_HexahedronElement(const int p, const int btype)
   : NodalTensorFiniteElement(3, p, VerifyOpen(btype), L2_DOF_MAP)
{
   const double *op = poly1d.OpenPoints(p, btype);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
#endif

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            Nodes.IntPoint(o++).Set3(op[i], op[j], op[k]);
         }
}

void L2_HexahedronElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), shape_z(p+1);
#endif

   basis1d.Eval(ip.x, shape_x);
   basis1d.Eval(ip.y, shape_y);
   basis1d.Eval(ip.z, shape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            shape(o++) = shape_x(i)*shape_y(j)*shape_z(k);
         }
}

void L2_HexahedronElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1),  shape_y(p+1),  shape_z(p+1);
   Vector dshape_x(p+1), dshape_y(p+1), dshape_z(p+1);
#endif

   basis1d.Eval(ip.x, shape_x, dshape_x);
   basis1d.Eval(ip.y, shape_y, dshape_y);
   basis1d.Eval(ip.z, shape_z, dshape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dshape(o,0) = dshape_x(i)* shape_y(j)* shape_z(k);
            dshape(o,1) =  shape_x(i)*dshape_y(j)* shape_z(k);
            dshape(o,2) =  shape_x(i)* shape_y(j)*dshape_z(k);  o++;
         }
}

void L2_HexahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;
   const double *op = poly1d.OpenPoints(p, b_type);

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1);
#endif

   for (int i = 0; i <= p; i++)
   {
      shape_x(i) = poly1d.CalcDelta(p,(1.0 - op[i]));
      shape_y(i) = poly1d.CalcDelta(p,op[i]);
   }

   switch (vertex)
   {
      case 0:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_x(i)*shape_x(j)*shape_x(k);
               }
         break;
      case 1:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_y(i)*shape_x(j)*shape_x(k);
               }
         break;
      case 2:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_y(i)*shape_y(j)*shape_x(k);
               }
         break;
      case 3:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_x(i)*shape_y(j)*shape_x(k);
               }
         break;
      case 4:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_x(i)*shape_x(j)*shape_y(k);
               }
         break;
      case 5:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_y(i)*shape_x(j)*shape_y(k);
               }
         break;
      case 6:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_y(i)*shape_y(j)*shape_y(k);
               }
         break;
      case 7:
         for (int o = 0, k = 0; k <= p; k++)
            for (int j = 0; j <= p; j++)
               for (int i = 0; i <= p; i++)
               {
                  dofs[o++] = shape_x(i)*shape_y(j)*shape_y(k);
               }
         break;
   }
}


L2Pos_HexahedronElement::L2Pos_HexahedronElement(const int p)
   : PositiveTensorFiniteElement(3, p, L2_DOF_MAP)
{
#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
#endif

   if (p == 0)
   {
      Nodes.IntPoint(0).Set3(0.5, 0.5, 0.5);
   }
   else
   {
      for (int o = 0, k = 0; k <= p; k++)
         for (int j = 0; j <= p; j++)
            for (int i = 0; i <= p; i++)
            {
               Nodes.IntPoint(o++).Set3(double(i)/p, double(j)/p, double(k)/p);
            }
   }
}

void L2Pos_HexahedronElement::CalcShape(const IntegrationPoint &ip,
                                        Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1), shape_y(p+1), shape_z(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x);
   Poly_1D::CalcBernstein(p, ip.y, shape_y);
   Poly_1D::CalcBernstein(p, ip.z, shape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            shape(o++) = shape_x(i)*shape_y(j)*shape_z(k);
         }
}

void L2Pos_HexahedronElement::CalcDShape(const IntegrationPoint &ip,
                                         DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p+1),  shape_y(p+1),  shape_z(p+1);
   Vector dshape_x(p+1), dshape_y(p+1), dshape_z(p+1);
#endif

   Poly_1D::CalcBernstein(p, ip.x, shape_x, dshape_x);
   Poly_1D::CalcBernstein(p, ip.y, shape_y, dshape_y);
   Poly_1D::CalcBernstein(p, ip.z, shape_z, dshape_z);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dshape(o,0) = dshape_x(i)* shape_y(j)* shape_z(k);
            dshape(o,1) =  shape_x(i)*dshape_y(j)* shape_z(k);
            dshape(o,2) =  shape_x(i)* shape_y(j)*dshape_z(k);  o++;
         }
}

void L2Pos_HexahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   const int p = order;

   dofs = 0.0;
   switch (vertex)
   {
      case 0: dofs[0] = 1.0; break;
      case 1: dofs[p] = 1.0; break;
      case 2: dofs[p*(p + 2)] = 1.0; break;
      case 3: dofs[p*(p + 1)] = 1.0; break;
      case 4: dofs[p*(p + 1)*(p + 1)] = 1.0; break;
      case 5: dofs[p + p*(p + 1)*(p + 1)] = 1.0; break;
      case 6: dofs[dof - 1] = 1.0; break;
      case 7: dofs[dof - p - 1] = 1.0; break;
   }
}


L2_TriangleElement::L2_TriangleElement(const int p, const int btype)
   : NodalFiniteElement(2, Geometry::TRIANGLE, ((p + 1)*(p + 2))/2, p,
                        FunctionSpace::Pk)
{
   const double *op = poly1d.OpenPoints(p, VerifyNodal(VerifyOpen(btype)));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
#endif

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         double w = op[i] + op[j] + op[p-i-j];
         Nodes.IntPoint(o++).Set2(op[i]/w, op[j]/w);
      }

   DenseMatrix T(dof);
   for (int k = 0; k < dof; k++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(k);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);

      for (int o = 0, j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            T(o++, k) = shape_x(i)*shape_y(j)*shape_l(p-i-j);
         }
   }

   Ti.Factor(T);
   // mfem::out << "L2_TriangleElement(" << p << ") : "; Ti.TestInversion();
}

void L2_TriangleElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1), u(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         u(o++) = shape_x(i)*shape_y(j)*shape_l(p-i-j);
      }

   Ti.Mult(u, shape);
}

void L2_TriangleElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_l(p + 1);
   DenseMatrix du(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l, dshape_l);

   for (int o = 0, j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         int k = p - i - j;
         du(o,0) = ((dshape_x(i)* shape_l(k)) -
                    ( shape_x(i)*dshape_l(k)))*shape_y(j);
         du(o,1) = ((dshape_y(j)* shape_l(k)) -
                    ( shape_y(j)*dshape_l(k)))*shape_x(i);
         o++;
      }

   Ti.Mult(du, dshape);
}

void L2_TriangleElement::ProjectDelta(int vertex, Vector &dofs) const
{
   switch (vertex)
   {
      case 0:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(1.0 - ip.x - ip.y, order);
         }
         break;
      case 1:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.x, order);
         }
         break;
      case 2:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.y, order);
         }
         break;
   }
}


L2Pos_TriangleElement::L2Pos_TriangleElement(const int p)
   : PositiveFiniteElement(2, Geometry::TRIANGLE, ((p + 1)*(p + 2))/2, p,
                           FunctionSpace::Pk)
{
#ifndef MFEM_THREAD_SAFE
   dshape_1d.SetSize(p + 1);
#endif

   if (p == 0)
   {
      Nodes.IntPoint(0).Set2(1./3, 1./3);
   }
   else
   {
      for (int o = 0, j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            Nodes.IntPoint(o++).Set2(double(i)/p, double(j)/p);
         }
   }
}

void L2Pos_TriangleElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   H1Pos_TriangleElement::CalcShape(order, ip.x, ip.y, shape.GetData());
}

void L2Pos_TriangleElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector dshape_1d(order + 1);
#endif

   H1Pos_TriangleElement::CalcDShape(order, ip.x, ip.y, dshape_1d.GetData(),
                                     dshape.Data());
}

void L2Pos_TriangleElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   switch (vertex)
   {
      case 0: dofs[0] = 1.0; break;
      case 1: dofs[order] = 1.0; break;
      case 2: dofs[dof-1] = 1.0; break;
   }
}


L2_TetrahedronElement::L2_TetrahedronElement(const int p, const int btype)
   : NodalFiniteElement(3, Geometry::TETRAHEDRON, ((p + 1)*(p + 2)*(p + 3))/6,
                        p, FunctionSpace::Pk)
{
   const double *op = poly1d.OpenPoints(p, VerifyNodal(VerifyOpen(btype)));

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
#endif

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            double w = op[i] + op[j] + op[k] + op[p-i-j-k];
            Nodes.IntPoint(o++).Set3(op[i]/w, op[j]/w, op[k]/w);
         }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);

      for (int o = 0, k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               T(o++, m) = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
            }
   }

   Ti.Factor(T);
   // mfem::out << "L2_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void L2_TetrahedronElement::CalcShape(const IntegrationPoint &ip,
                                      Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   Vector u(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
         }

   Ti.Mult(u, shape);
}

void L2_TetrahedronElement::CalcDShape(const IntegrationPoint &ip,
                                       DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   DenseMatrix du(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l);

   for (int o = 0, k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            int l = p - i - j - k;
            du(o,0) = ((dshape_x(i)* shape_l(l)) -
                       ( shape_x(i)*dshape_l(l)))*shape_y(j)*shape_z(k);
            du(o,1) = ((dshape_y(j)* shape_l(l)) -
                       ( shape_y(j)*dshape_l(l)))*shape_x(i)*shape_z(k);
            du(o,2) = ((dshape_z(k)* shape_l(l)) -
                       ( shape_z(k)*dshape_l(l)))*shape_x(i)*shape_y(j);
            o++;
         }

   Ti.Mult(du, dshape);
}

void L2_TetrahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   switch (vertex)
   {
      case 0:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(1.0 - ip.x - ip.y - ip.z, order);
         }
         break;
      case 1:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.x, order);
         }
         break;
      case 2:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.y, order);
         }
         break;
      case 3:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.z, order);
         }
         break;
   }
}

L2_PentatopeElement::L2_PentatopeElement(const int p, const int _type)
   : NodalFiniteElement(4, Geometry::PENTATOPE,
                        ((p + 1)*(p + 2)*(p + 3)*(p + 4))/24,
                        p, FunctionSpace::Pk), T(dof)
{
   const double *op;

   type = _type;
   switch (type)
   {
      case 0: op = poly1d.OpenPoints(p); break;
      case 1:
      default: op = poly1d.ClosedPoints(p);
   }

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof);
   du.SetSize(dof, dim);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
#endif

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; l + k <= p; k++)
         for (int j = 0; j + l + k <= p; j++)
            for (int i = 0; i + j + l + k <= p; i++)
            {
               double w = op[i] + op[j] + op[k] + op[l] + op[p-i-j-k-l];
               Nodes.IntPoint(o++).Set4(op[i]/w, op[j]/w, op[k]/w, op[l]/w);
            }

   for (int m = 0; m < dof; m++)
   {
      IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, ip.t, shape_t);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

      for (int o = 0, l = 0; l <= p; l++)
         for (int k = 0; l + k <= p; k++)
            for (int j = 0; j + l + k <= p; j++)
               for (int i = 0; i + j + l + k <= p; i++)
               {
                  T(o++, m) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
               }
   }

   T.Invert();
}

void L2_PentatopeElement::CalcShape(const IntegrationPoint &ip,
                                    Vector &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   Vector u(Dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

   for (int o = 0, l = 0; l <= p; l++)
      for (int k = 0; l + k <= p; k++)
         for (int j = 0; j + l + k <= p; j++)
            for (int i = 0; i + j + l + k <= p; i++)
            {
               u(o++) = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
            }

   T.Mult(u, shape);
}

void L2_PentatopeElement::CalcDShape(const IntegrationPoint &ip,
                                     DenseMatrix &dshape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector  shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1), shape_t(p + 1),
           shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_t(p + 1),
          dshape_l(p + 1);
   DenseMatrix du(Dof, Dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, ip.t, shape_t, dshape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l, dshape_l);

   for (int o = 0, m = 0; m <= p; m++)
      for (int k = 0; k + m <= p; k++)
         for (int j = 0; j + k + m <= p; j++)
            for (int i = 0; i + j + k + m <= p; i++)
            {
               int l = p - i - j - k - m;
               du(o,0) = ((dshape_x(i)* shape_l(l)) -
                          ( shape_x(i)*dshape_l(l)))*shape_y(j)*shape_z(k)*shape_t(m);
               du(o,1) = ((dshape_y(j)* shape_l(l)) -
                          ( shape_y(j)*dshape_l(l)))*shape_x(i)*shape_z(k)*shape_t(m);
               du(o,2) = ((dshape_z(k)* shape_l(l)) -
                          ( shape_z(k)*dshape_l(l)))*shape_x(i)*shape_y(j)*shape_t(m);
               du(o,3) = ((dshape_t(m)* shape_l(l)) -
                          ( shape_t(m)*dshape_l(l)))*shape_x(i)*shape_y(j)*shape_z(k);
               o++;
            }

   Mult(T, du, dshape);
}

void L2_PentatopeElement::ProjectDelta(int vertex, Vector &dofs) const
{
   switch (vertex)
   {
      case 0:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(1.0 - ip.x - ip.y - ip.z - ip.t, order);
         }
         break;
      case 1:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.x, order);
         }
         break;
      case 2:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.y, order);
         }
         break;
      case 3:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.z, order);
         }
         break;
      case 4:
         for (int i = 0; i < dof; i++)
         {
            const IntegrationPoint &ip = Nodes.IntPoint(i);
            dofs[i] = pow(ip.t, order);
         }
         break;
   }
}


L2Pos_TetrahedronElement::L2Pos_TetrahedronElement(const int p)
   : PositiveFiniteElement(3, Geometry::TETRAHEDRON,
                           ((p + 1)*(p + 2)*(p + 3))/6, p, FunctionSpace::Pk)
{
#ifndef MFEM_THREAD_SAFE
   dshape_1d.SetSize(p + 1);
#endif

   if (p == 0)
   {
      Nodes.IntPoint(0).Set3(0.25, 0.25, 0.25);
   }
   else
   {
      for (int o = 0, k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               Nodes.IntPoint(o++).Set3(double(i)/p, double(j)/p, double(k)/p);
            }
   }
}

void L2Pos_TetrahedronElement::CalcShape(const IntegrationPoint &ip,
                                         Vector &shape) const
{
   H1Pos_TetrahedronElement::CalcShape(order, ip.x, ip.y, ip.z,
                                       shape.GetData());
}

void L2Pos_TetrahedronElement::CalcDShape(const IntegrationPoint &ip,
                                          DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector dshape_1d(order + 1);
#endif

   H1Pos_TetrahedronElement::CalcDShape(order, ip.x, ip.y, ip.z,
                                        dshape_1d.GetData(), dshape.Data());
}

void L2Pos_TetrahedronElement::ProjectDelta(int vertex, Vector &dofs) const
{
   dofs = 0.0;
   switch (vertex)
   {
      case 0: dofs[0] = 1.0; break;
      case 1: dofs[order] = 1.0; break;
      case 2: dofs[(order*(order+3))/2] = 1.0; break;
      case 3: dofs[dof-1] = 1.0; break;
   }
}


L2_WedgeElement::L2_WedgeElement(const int p, const int btype)
   : NodalFiniteElement(3, Geometry::PRISM, ((p + 1)*(p + 1)*(p + 2))/2,
                        p, FunctionSpace::Qk),
     TriangleFE(p, btype),
     SegmentFE(p, btype)
{
#ifndef MFEM_THREAD_SAFE
   t_shape.SetSize(TriangleFE.GetDof());
   s_shape.SetSize(SegmentFE.GetDof());
   t_dshape.SetSize(TriangleFE.GetDof(), 2);
   s_dshape.SetSize(SegmentFE.GetDof(), 1);
#endif

   t_dof.SetSize(dof);
   s_dof.SetSize(dof);

   // Interior DoFs
   int m=0;
   for (int k=0; k<=p; k++)
   {
      int l=0;
      for (int j=0; j<=p; j++)
      {
         for (int i=0; i<=j; i++)
         {
            t_dof[m] = l;
            s_dof[m] = k;
            l++; m++;
         }
      }
   }

   // Define Nodes
   const IntegrationRule & t_Nodes = TriangleFE.GetNodes();
   const IntegrationRule & s_Nodes = SegmentFE.GetNodes();
   for (int i=0; i<dof; i++)
   {
      Nodes.IntPoint(i).x = t_Nodes.IntPoint(t_dof[i]).x;
      Nodes.IntPoint(i).y = t_Nodes.IntPoint(t_dof[i]).y;
      Nodes.IntPoint(i).z = s_Nodes.IntPoint(s_dof[i]).x;
   }
}

void L2_WedgeElement::CalcShape(const IntegrationPoint &ip,
                                Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector t_shape(TriangleFE.GetDof());
   Vector s_shape(SegmentFE.GetDof());
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   SegmentFE.CalcShape(ipz, s_shape);

   for (int i=0; i<dof; i++)
   {
      shape[i] = t_shape[t_dof[i]] * s_shape[s_dof[i]];
   }
}

void L2_WedgeElement::CalcDShape(const IntegrationPoint &ip,
                                 DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector      t_shape(TriangleFE.GetDof());
   DenseMatrix t_dshape(TriangleFE.GetDof(), 2);
   Vector      s_shape(SegmentFE.GetDof());
   DenseMatrix s_dshape(SegmentFE.GetDof(), 1);
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   TriangleFE.CalcDShape(ip, t_dshape);
   SegmentFE.CalcShape(ipz, s_shape);
   SegmentFE.CalcDShape(ipz, s_dshape);

   for (int i=0; i<dof; i++)
   {
      dshape(i, 0) = t_dshape(t_dof[i],0) * s_shape[s_dof[i]];
      dshape(i, 1) = t_dshape(t_dof[i],1) * s_shape[s_dof[i]];
      dshape(i, 2) = t_shape[t_dof[i]] * s_dshape(s_dof[i],0);
   }
}


L2Pos_WedgeElement::L2Pos_WedgeElement(const int p)
   : PositiveFiniteElement(3, Geometry::PRISM,
                           ((p + 1)*(p + 1)*(p + 2))/2, p, FunctionSpace::Qk),
     TriangleFE(p),
     SegmentFE(p)
{
#ifndef MFEM_THREAD_SAFE
   t_shape.SetSize(TriangleFE.GetDof());
   s_shape.SetSize(SegmentFE.GetDof());
   t_dshape.SetSize(TriangleFE.GetDof(), 2);
   s_dshape.SetSize(SegmentFE.GetDof(), 1);
#endif

   t_dof.SetSize(dof);
   s_dof.SetSize(dof);

   // Interior DoFs
   int m=0;
   for (int k=0; k<=p; k++)
   {
      int l=0;
      for (int j=0; j<=p; j++)
      {
         for (int i=0; i<=j; i++)
         {
            t_dof[m] = l;
            s_dof[m] = k;
            l++; m++;
         }
      }
   }

   // Define Nodes
   const IntegrationRule & t_Nodes = TriangleFE.GetNodes();
   const IntegrationRule & s_Nodes = SegmentFE.GetNodes();
   for (int i=0; i<dof; i++)
   {
      Nodes.IntPoint(i).x = t_Nodes.IntPoint(t_dof[i]).x;
      Nodes.IntPoint(i).y = t_Nodes.IntPoint(t_dof[i]).y;
      Nodes.IntPoint(i).z = s_Nodes.IntPoint(s_dof[i]).x;
   }
}

void L2Pos_WedgeElement::CalcShape(const IntegrationPoint &ip,
                                   Vector &shape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector t_shape(TriangleFE.GetDof());
   Vector s_shape(SegmentFE.GetDof());
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   SegmentFE.CalcShape(ipz, s_shape);

   for (int i=0; i<dof; i++)
   {
      shape[i] = t_shape[t_dof[i]] * s_shape[s_dof[i]];
   }
}

void L2Pos_WedgeElement::CalcDShape(const IntegrationPoint &ip,
                                    DenseMatrix &dshape) const
{
#ifdef MFEM_THREAD_SAFE
   Vector      t_shape(TriangleFE.GetDof());
   DenseMatrix t_dshape(TriangleFE.GetDof(), 2);
   Vector      s_shape(SegmentFE.GetDof());
   DenseMatrix s_dshape(SegmentFE.GetDof(), 1);
#endif

   IntegrationPoint ipz; ipz.x = ip.z; ipz.y = 0.0; ipz.z = 0.0;

   TriangleFE.CalcShape(ip, t_shape);
   TriangleFE.CalcDShape(ip, t_dshape);
   SegmentFE.CalcShape(ipz, s_shape);
   SegmentFE.CalcDShape(ipz, s_dshape);

   for (int i=0; i<dof; i++)
   {
      dshape(i, 0) = t_dshape(t_dof[i],0) * s_shape[s_dof[i]];
      dshape(i, 1) = t_dshape(t_dof[i],1) * s_shape[s_dof[i]];
      dshape(i, 2) = t_shape[t_dof[i]] * s_dshape(s_dof[i],0);
   }
}


const double RT_QuadrilateralElement::nk[8] =
{ 0., -1.,  1., 0.,  0., 1.,  -1., 0. };

RT_QuadrilateralElement::RT_QuadrilateralElement(const int p,
                                                 const int cb_type,
                                                 const int ob_type)
   : VectorTensorFiniteElement(2, 2*(p + 1)*(p + 2), p + 1, cb_type, ob_type,
                               H_DIV, DofMapType::L2_DOF_MAP),
     dof2nk(dof)
{
   dof_map.SetSize(dof);

   const double *cp = poly1d.ClosedPoints(p + 1, cb_type);
   const double *op = poly1d.OpenPoints(p, ob_type);
   const int dof2 = dof/2;

#ifndef MFEM_THREAD_SAFE
   shape_cx.SetSize(p + 2);
   shape_ox.SetSize(p + 1);
   shape_cy.SetSize(p + 2);
   shape_oy.SetSize(p + 1);
   dshape_cx.SetSize(p + 2);
   dshape_cy.SetSize(p + 2);
#endif

   // edges
   int o = 0;
   for (int i = 0; i <= p; i++)  // (0,1)
   {
      dof_map[1*dof2 + i + 0*(p + 1)] = o++;
   }
   for (int i = 0; i <= p; i++)  // (1,2)
   {
      dof_map[0*dof2 + (p + 1) + i*(p + 2)] = o++;
   }
   for (int i = 0; i <= p; i++)  // (2,3)
   {
      dof_map[1*dof2 + (p - i) + (p + 1)*(p + 1)] = o++;
   }
   for (int i = 0; i <= p; i++)  // (3,0)
   {
      dof_map[0*dof2 + 0 + (p - i)*(p + 2)] = o++;
   }

   // interior
   for (int j = 0; j <= p; j++)  // x-components
      for (int i = 1; i <= p; i++)
      {
         dof_map[0*dof2 + i + j*(p + 2)] = o++;
      }
   for (int j = 1; j <= p; j++)  // y-components
      for (int i = 0; i <= p; i++)
      {
         dof_map[1*dof2 + i + j*(p + 1)] = o++;
      }

   // dof orientations
   // x-components
   for (int j = 0; j <= p; j++)
      for (int i = 0; i <= p/2; i++)
      {
         int idx = 0*dof2 + i + j*(p + 2);
         dof_map[idx] = -1 - dof_map[idx];
      }
   if (p%2 == 1)
      for (int j = p/2 + 1; j <= p; j++)
      {
         int idx = 0*dof2 + (p/2 + 1) + j*(p + 2);
         dof_map[idx] = -1 - dof_map[idx];
      }
   // y-components
   for (int j = 0; j <= p/2; j++)
      for (int i = 0; i <= p; i++)
      {
         int idx = 1*dof2 + i + j*(p + 1);
         dof_map[idx] = -1 - dof_map[idx];
      }
   if (p%2 == 1)
      for (int i = 0; i <= p/2; i++)
      {
         int idx = 1*dof2 + i + (p/2 + 1)*(p + 1);
         dof_map[idx] = -1 - dof_map[idx];
      }

   o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i <= p + 1; i++)
      {
         int idx;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx;
            dof2nk[idx] = 3;
         }
         else
         {
            dof2nk[idx] = 1;
         }
         Nodes.IntPoint(idx).Set2(cp[i], op[j]);
      }
   for (int j = 0; j <= p + 1; j++)
      for (int i = 0; i <= p; i++)
      {
         int idx;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx;
            dof2nk[idx] = 0;
         }
         else
         {
            dof2nk[idx] = 2;
         }
         Nodes.IntPoint(idx).Set2(op[i], cp[j]);
      }
}

void RT_QuadrilateralElement::CalcVShape(const IntegrationPoint &ip,
                                         DenseMatrix &shape) const
{
   const int pp1 = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(pp1 + 1), shape_ox(pp1), shape_cy(pp1 + 1), shape_oy(pp1);
#endif

   cbasis1d.Eval(ip.x, shape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy);
   obasis1d.Eval(ip.y, shape_oy);

   int o = 0;
   for (int j = 0; j < pp1; j++)
      for (int i = 0; i <= pp1; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         shape(idx,0) = s*shape_cx(i)*shape_oy(j);
         shape(idx,1) = 0.;
      }
   for (int j = 0; j <= pp1; j++)
      for (int i = 0; i < pp1; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         shape(idx,0) = 0.;
         shape(idx,1) = s*shape_ox(i)*shape_cy(j);
      }
}

void RT_QuadrilateralElement::CalcDivShape(const IntegrationPoint &ip,
                                           Vector &divshape) const
{
   const int pp1 = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(pp1 + 1), shape_ox(pp1), shape_cy(pp1 + 1), shape_oy(pp1);
   Vector dshape_cx(pp1 + 1), dshape_cy(pp1 + 1);
#endif

   cbasis1d.Eval(ip.x, shape_cx, dshape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy, dshape_cy);
   obasis1d.Eval(ip.y, shape_oy);

   int o = 0;
   for (int j = 0; j < pp1; j++)
      for (int i = 0; i <= pp1; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         divshape(idx) = s*dshape_cx(i)*shape_oy(j);
      }
   for (int j = 0; j <= pp1; j++)
      for (int i = 0; i < pp1; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         divshape(idx) = s*shape_ox(i)*dshape_cy(j);
      }
}


const double RT_HexahedronElement::nk[18] =
{ 0.,0.,-1.,  0.,-1.,0.,  1.,0.,0.,  0.,1.,0.,  -1.,0.,0.,  0.,0.,1. };

RT_HexahedronElement::RT_HexahedronElement(const int p,
                                           const int cb_type,
                                           const int ob_type)
   : VectorTensorFiniteElement(3, 3*(p + 1)*(p + 1)*(p + 2), p + 1, cb_type,
                               ob_type, H_DIV, DofMapType::L2_DOF_MAP),
     dof2nk(dof)
{
   dof_map.SetSize(dof);

   const double *cp = poly1d.ClosedPoints(p + 1, cb_type);
   const double *op = poly1d.OpenPoints(p, ob_type);
   const int dof3 = dof/3;

#ifndef MFEM_THREAD_SAFE
   shape_cx.SetSize(p + 2);
   shape_ox.SetSize(p + 1);
   shape_cy.SetSize(p + 2);
   shape_oy.SetSize(p + 1);
   shape_cz.SetSize(p + 2);
   shape_oz.SetSize(p + 1);
   dshape_cx.SetSize(p + 2);
   dshape_cy.SetSize(p + 2);
   dshape_cz.SetSize(p + 2);
#endif

   // faces
   int o = 0;
   for (int j = 0; j <= p; j++)  // (3,2,1,0) -- bottom
      for (int i = 0; i <= p; i++)
      {
         dof_map[2*dof3 + i + ((p - j) + 0*(p + 1))*(p + 1)] = o++;
      }
   for (int j = 0; j <= p; j++)  // (0,1,5,4) -- front
      for (int i = 0; i <= p; i++)
      {
         dof_map[1*dof3 + i + (0 + j*(p + 2))*(p + 1)] = o++;
      }
   for (int j = 0; j <= p; j++)  // (1,2,6,5) -- right
      for (int i = 0; i <= p; i++)
      {
         dof_map[0*dof3 + (p + 1) + (i + j*(p + 1))*(p + 2)] = o++;
      }
   for (int j = 0; j <= p; j++)  // (2,3,7,6) -- back
      for (int i = 0; i <= p; i++)
      {
         dof_map[1*dof3 + (p - i) + ((p + 1) + j*(p + 2))*(p + 1)] = o++;
      }
   for (int j = 0; j <= p; j++)  // (3,0,4,7) -- left
      for (int i = 0; i <= p; i++)
      {
         dof_map[0*dof3 + 0 + ((p - i) + j*(p + 1))*(p + 2)] = o++;
      }
   for (int j = 0; j <= p; j++)  // (4,5,6,7) -- top
      for (int i = 0; i <= p; i++)
      {
         dof_map[2*dof3 + i + (j + (p + 1)*(p + 1))*(p + 1)] = o++;
      }

   // interior
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 1; i <= p; i++)
         {
            dof_map[0*dof3 + i + (j + k*(p + 1))*(p + 2)] = o++;
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 1; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dof_map[1*dof3 + i + (j + k*(p + 2))*(p + 1)] = o++;
         }
   // z-components
   for (int k = 1; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            dof_map[2*dof3 + i + (j + k*(p + 1))*(p + 1)] = o++;
         }

   // dof orientations
   // for odd p, do not change the orientations in the mid-planes
   // {i = p/2 + 1}, {j = p/2 + 1}, {k = p/2 + 1} in the x, y, z-components
   // respectively.
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p/2; i++)
         {
            int idx = 0*dof3 + i + (j + k*(p + 1))*(p + 2);
            dof_map[idx] = -1 - dof_map[idx];
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p/2; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx = 1*dof3 + i + (j + k*(p + 2))*(p + 1);
            dof_map[idx] = -1 - dof_map[idx];
         }
   // z-components
   for (int k = 0; k <= p/2; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx = 2*dof3 + i + (j + k*(p + 1))*(p + 1);
            dof_map[idx] = -1 - dof_map[idx];
         }

   o = 0;
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p + 1; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx;
               dof2nk[idx] = 4;
            }
            else
            {
               dof2nk[idx] = 2;
            }
            Nodes.IntPoint(idx).Set3(cp[i], op[j], op[k]);
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p + 1; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx;
               dof2nk[idx] = 1;
            }
            else
            {
               dof2nk[idx] = 3;
            }
            Nodes.IntPoint(idx).Set3(op[i], cp[j], op[k]);
         }
   // z-components
   for (int k = 0; k <= p + 1; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx;
               dof2nk[idx] = 0;
            }
            else
            {
               dof2nk[idx] = 5;
            }
            Nodes.IntPoint(idx).Set3(op[i], op[j], cp[k]);
         }
}

void RT_HexahedronElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   const int pp1 = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(pp1 + 1), shape_ox(pp1), shape_cy(pp1 + 1), shape_oy(pp1);
   Vector shape_cz(pp1 + 1), shape_oz(pp1);
#endif

   cbasis1d.Eval(ip.x, shape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy);
   obasis1d.Eval(ip.y, shape_oy);
   cbasis1d.Eval(ip.z, shape_cz);
   obasis1d.Eval(ip.z, shape_oz);

   int o = 0;
   // x-components
   for (int k = 0; k < pp1; k++)
      for (int j = 0; j < pp1; j++)
         for (int i = 0; i <= pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = s*shape_cx(i)*shape_oy(j)*shape_oz(k);
            shape(idx,1) = 0.;
            shape(idx,2) = 0.;
         }
   // y-components
   for (int k = 0; k < pp1; k++)
      for (int j = 0; j <= pp1; j++)
         for (int i = 0; i < pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = 0.;
            shape(idx,1) = s*shape_ox(i)*shape_cy(j)*shape_oz(k);
            shape(idx,2) = 0.;
         }
   // z-components
   for (int k = 0; k <= pp1; k++)
      for (int j = 0; j < pp1; j++)
         for (int i = 0; i < pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = 0.;
            shape(idx,1) = 0.;
            shape(idx,2) = s*shape_ox(i)*shape_oy(j)*shape_cz(k);
         }
}

void RT_HexahedronElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   const int pp1 = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(pp1 + 1), shape_ox(pp1), shape_cy(pp1 + 1), shape_oy(pp1);
   Vector shape_cz(pp1 + 1), shape_oz(pp1);
   Vector dshape_cx(pp1 + 1), dshape_cy(pp1 + 1), dshape_cz(pp1 + 1);
#endif

   cbasis1d.Eval(ip.x, shape_cx, dshape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy, dshape_cy);
   obasis1d.Eval(ip.y, shape_oy);
   cbasis1d.Eval(ip.z, shape_cz, dshape_cz);
   obasis1d.Eval(ip.z, shape_oz);

   int o = 0;
   // x-components
   for (int k = 0; k < pp1; k++)
      for (int j = 0; j < pp1; j++)
         for (int i = 0; i <= pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            divshape(idx) = s*dshape_cx(i)*shape_oy(j)*shape_oz(k);
         }
   // y-components
   for (int k = 0; k < pp1; k++)
      for (int j = 0; j <= pp1; j++)
         for (int i = 0; i < pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            divshape(idx) = s*shape_ox(i)*dshape_cy(j)*shape_oz(k);
         }
   // z-components
   for (int k = 0; k <= pp1; k++)
      for (int j = 0; j < pp1; j++)
         for (int i = 0; i < pp1; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            divshape(idx) = s*shape_ox(i)*shape_oy(j)*dshape_cz(k);
         }
}


const double RT_TriangleElement::nk[6] =
{ 0., -1., 1., 1., -1., 0. };

const double RT_TriangleElement::c = 1./3.;

RT_TriangleElement::RT_TriangleElement(const int p)
   : VectorFiniteElement(2, Geometry::TRIANGLE, (p + 1)*(p + 3), p + 1,
                         H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   const double *iop = (p > 0) ? poly1d.OpenPoints(p - 1) : NULL;
   const double *bop = poly1d.OpenPoints(p);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim);
   divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
#endif

   // edges
   int o = 0;
   for (int i = 0; i <= p; i++)  // (0,1)
   {
      Nodes.IntPoint(o).Set2(bop[i], 0.);
      dof2nk[o++] = 0;
   }
   for (int i = 0; i <= p; i++)  // (1,2)
   {
      Nodes.IntPoint(o).Set2(bop[p-i], bop[i]);
      dof2nk[o++] = 1;
   }
   for (int i = 0; i <= p; i++)  // (2,0)
   {
      Nodes.IntPoint(o).Set2(0., bop[p-i]);
      dof2nk[o++] = 2;
   }

   // interior
   for (int j = 0; j < p; j++)
      for (int i = 0; i + j < p; i++)
      {
         double w = iop[i] + iop[j] + iop[p-1-i-j];
         Nodes.IntPoint(o).Set2(iop[i]/w, iop[j]/w);
         dof2nk[o++] = 0;
         Nodes.IntPoint(o).Set2(iop[i]/w, iop[j]/w);
         dof2nk[o++] = 2;
      }

   DenseMatrix T(dof);
   for (int k = 0; k < dof; k++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(k);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);
      const double *n_k = nk + 2*dof2nk[k];

      o = 0;
      for (int j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_l(p-i-j);
            T(o++, k) = s*n_k[0];
            T(o++, k) = s*n_k[1];
         }
      for (int i = 0; i <= p; i++)
      {
         double s = shape_x(i)*shape_y(p-i);
         T(o++, k) = s*((ip.x - c)*n_k[0] + (ip.y - c)*n_k[1]);
      }
   }

   Ti.Factor(T);
   // mfem::out << "RT_TriangleElement(" << p << ") : "; Ti.TestInversion();
}

void RT_TriangleElement::CalcVShape(const IntegrationPoint &ip,
                                    DenseMatrix &shape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
   DenseMatrix u(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);

   int o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         double s = shape_x(i)*shape_y(j)*shape_l(p-i-j);
         u(o,0) = s;  u(o,1) = 0;  o++;
         u(o,0) = 0;  u(o,1) = s;  o++;
      }
   for (int i = 0; i <= p; i++)
   {
      double s = shape_x(i)*shape_y(p-i);
      u(o,0) = (ip.x - c)*s;
      u(o,1) = (ip.y - c)*s;
      o++;
   }

   Ti.Mult(u, shape);
}

void RT_TriangleElement::CalcDivShape(const IntegrationPoint &ip,
                                      Vector &divshape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1),  shape_y(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_l(p + 1);
   Vector divu(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l, dshape_l);

   int o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         int k = p - i - j;
         divu(o++) = (dshape_x(i)*shape_l(k) -
                      shape_x(i)*dshape_l(k))*shape_y(j);
         divu(o++) = (dshape_y(j)*shape_l(k) -
                      shape_y(j)*dshape_l(k))*shape_x(i);
      }
   for (int i = 0; i <= p; i++)
   {
      int j = p - i;
      divu(o++) = ((shape_x(i) + (ip.x - c)*dshape_x(i))*shape_y(j) +
                   (shape_y(j) + (ip.y - c)*dshape_y(j))*shape_x(i));
   }

   Ti.Mult(divu, divshape);
}


const double RT_TetrahedronElement::nk[12] =
{ 1,1,1,  -1,0,0,  0,-1,0,  0,0,-1 };
// { .5,.5,.5, -.5,0,0, 0,-.5,0, 0,0,-.5}; // n_F |F|

const double RT_TetrahedronElement::c = 1./4.;

RT_TetrahedronElement::RT_TetrahedronElement(const int p)
   : VectorFiniteElement(3, Geometry::TETRAHEDRON, (p + 1)*(p + 2)*(p + 4)/2,
                         p + 1, H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   const double *iop = (p > 0) ? poly1d.OpenPoints(p - 1) : NULL;
   const double *bop = poly1d.OpenPoints(p);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim);
   divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
#endif

   int o = 0;
   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp,
   //        the constructor of H1_TetrahedronElement)
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)  // (1,2,3)
      {
         double w = bop[i] + bop[j] + bop[p-i-j];
         Nodes.IntPoint(o).Set3(bop[p-i-j]/w, bop[i]/w, bop[j]/w);
         dof2nk[o++] = 0;
      }
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)  // (0,3,2)
      {
         double w = bop[i] + bop[j] + bop[p-i-j];
         Nodes.IntPoint(o).Set3(0., bop[j]/w, bop[i]/w);
         dof2nk[o++] = 1;
      }
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)  // (0,1,3)
      {
         double w = bop[i] + bop[j] + bop[p-i-j];
         Nodes.IntPoint(o).Set3(bop[i]/w, 0., bop[j]/w);
         dof2nk[o++] = 2;
      }
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)  // (0,2,1)
      {
         double w = bop[i] + bop[j] + bop[p-i-j];
         Nodes.IntPoint(o).Set3(bop[j]/w, bop[i]/w, 0.);
         dof2nk[o++] = 3;
      }

   // interior
   for (int k = 0; k < p; k++)
      for (int j = 0; j + k < p; j++)
         for (int i = 0; i + j + k < p; i++)
         {
            double w = iop[i] + iop[j] + iop[k] + iop[p-1-i-j-k];
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2nk[o++] = 1;
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2nk[o++] = 2;
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2nk[o++] = 3;
         }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);
      const double *nm = nk + 3*dof2nk[m];

      o = 0;
      for (int k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
               T(o++, m) = s * nm[0];
               T(o++, m) = s * nm[1];
               T(o++, m) = s * nm[2];
            }
      for (int j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_z(p-i-j);
            T(o++, m) = s*((ip.x - c)*nm[0] + (ip.y - c)*nm[1] +
                           (ip.z - c)*nm[2]);
         }
   }

   Ti.Factor(T);
   // mfem::out << "RT_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void RT_TetrahedronElement::CalcVShape(const IntegrationPoint &ip,
                                       DenseMatrix &shape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   DenseMatrix u(dof, dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);

   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
            u(o,0) = s;  u(o,1) = 0;  u(o,2) = 0;  o++;
            u(o,0) = 0;  u(o,1) = s;  u(o,2) = 0;  o++;
            u(o,0) = 0;  u(o,1) = 0;  u(o,2) = s;  o++;
         }
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         double s = shape_x(i)*shape_y(j)*shape_z(p-i-j);
         u(o,0) = (ip.x - c)*s;  u(o,1) = (ip.y - c)*s;  u(o,2) = (ip.z - c)*s;
         o++;
      }

   Ti.Mult(u, shape);
}

void RT_TetrahedronElement::CalcDivShape(const IntegrationPoint &ip,
                                         Vector &divshape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   Vector divu(dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l);

   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            int l = p - i - j - k;
            divu(o++) = (dshape_x(i)*shape_l(l) -
                         shape_x(i)*dshape_l(l))*shape_y(j)*shape_z(k);
            divu(o++) = (dshape_y(j)*shape_l(l) -
                         shape_y(j)*dshape_l(l))*shape_x(i)*shape_z(k);
            divu(o++) = (dshape_z(k)*shape_l(l) -
                         shape_z(k)*dshape_l(l))*shape_x(i)*shape_y(j);
         }
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         int k = p - i - j;
         divu(o++) =
            (shape_x(i) + (ip.x - c)*dshape_x(i))*shape_y(j)*shape_z(k) +
            (shape_y(j) + (ip.y - c)*dshape_y(j))*shape_x(i)*shape_z(k) +
            (shape_z(k) + (ip.z - c)*dshape_z(k))*shape_x(i)*shape_y(j);
      }

   Ti.Mult(divu, divshape);
}

const double RT_PentatopeElement::nk[20] =
{ 0,0,0,-1,  0,0,-1,0,  0,-1,0,0,  -1,0,0,0,  1,1,1,1};
// { .5,.5,.5, -.5,0,0, 0,-.5,0, 0,0,-.5}; // n_F |F|

const double RT_PentatopeElement::c = 1./5.;

RT_PentatopeElement::RT_PentatopeElement(const int p)
   : VectorFiniteElement(4, Geometry::PENTATOPE, (p + 1)*(p + 2)*(p + 3)*(p + 5)/6,
                         p + 1, H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   const double *iop = (p > 0) ? poly1d.OpenPoints(p - 1) : NULL;
   const double *bop = poly1d.OpenPoints(p);

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1);
   shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1);
   shape_t.SetSize(p + 1);
   shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1);
   dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1);
   dshape_t.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim);
   divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
#endif

   int o = 0;
   // faces (see Mesh::GenerateFaces in mesh/mesh.cpp,
   //        the constructor of H1_PentatopeElement)
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,2,3)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, bop[j]/w, bop[k]/w, 0.);
            dof2nk[o++] = 0;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,2,1,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[j]/w, bop[i]/w, 0., bop[k]/w);
            dof2nk[o++] = 1;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,1,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[i]/w, 0., bop[j]/w, bop[k]/w);
            dof2nk[o++] = 2;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (0,3,2,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(0., bop[j]/w, bop[i]/w, bop[k]/w);
            dof2nk[o++] = 3;
         }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)  // (1,2,3,4)
         {
            double w = bop[i] + bop[j] + bop[k] + bop[p-i-j-k];
            Nodes.IntPoint(o).Set4(bop[p-i-j-k]/w, bop[i]/w, bop[j]/w, bop[k]/w);
            dof2nk[o++] = 4;
         }

   // interior
   for (int l = 0; l < p; l++)
      for (int k = 0; k + l < p; k++)
         for (int j = 0; j + k + l < p; j++)
            for (int i = 0; i + j + k + l < p; i++)
            {
               double w = iop[i] + iop[j] + iop[k] + iop[l] + iop[p-1-i-j-k-l];
               Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
               dof2nk[o++] = 1;
               Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
               dof2nk[o++] = 2;
               Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
               dof2nk[o++] = 3;
               Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
               dof2nk[o++] = 4;
            }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, ip.t, shape_t);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);
      const double *nm = nk + 4*dof2nk[m];

      o = 0;
      for (int l = 0; l<= p; l++)
         for (int k = 0; k + l <= p; k++)
            for (int j = 0; j + k + l <= p; j++)
               for (int i = 0; i + j + k + l <= p; i++)
               {
                  double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
                  T(o++, m) = s * nm[0];
                  T(o++, m) = s * nm[1];
                  T(o++, m) = s * nm[2];
                  T(o++, m) = s * nm[3];
               }
      for (int k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(p-i-j-k);
               T(o++, m) = s*((ip.x - c)*nm[0] + (ip.y - c)*nm[1] +
                              (ip.z - c)*nm[2] + (ip.t - c)*nm[3]);
            }
   }
    std::ofstream Vander_file("Vander_matrix_Baseline.txt");
    T.PrintMatlab(Vander_file);

   Ti.Factor(T);
   // mfem::out << "RT_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void RT_PentatopeElement::CalcVShape(const IntegrationPoint &ip,
                                     DenseMatrix &shape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
          shape_l(p + 1);
   DenseMatrix u(Dof, Dim);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, ip.t, shape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l);

   int o = 0;
   for (int l = 0; l <= p; l++)
      for (int k = 0; k + l <= p; k++)
         for (int j = 0; j + k + l <= p; j++)
            for (int i = 0; i + j + k + l <= p; i++)
            {
               double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(l)*shape_l(p-i-j-k-l);
               u(o,0) = s;  u(o,1) = 0;  u(o,2) = 0;  u(o,3) = 0;  o++;
               u(o,0) = 0;  u(o,1) = s;  u(o,2) = 0;  u(o,3) = 0;  o++;
               u(o,0) = 0;  u(o,1) = 0;  u(o,2) = s;  u(o,3) = 0;  o++;
               u(o,0) = 0;  u(o,1) = 0;  u(o,2) = 0;  u(o,3) = s;  o++;
            }
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_t(p-i-j-k);
            u(o,0) = (ip.x - c)*s;  u(o,1) = (ip.y - c)*s;  u(o,2) = (ip.z - c)*s;
            u(o,3) = (ip.t - c)*s;
            o++;
         }
//    std::ofstream u_file("u.txt");
//    u.PrintMatlab(u_file);

   Ti.Mult(u, shape);
    
//    std::ofstream shape_file("shape_baseline.txt");
//    shape.PrintMatlab(shape_file);
}

void RT_PentatopeElement::CalcDivShape(const IntegrationPoint &ip,
                                       Vector &divshape) const
{
   const int p = order - 1;

#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1),  shape_y(p + 1),  shape_z(p + 1),  shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   Vector divu(Dof);
#endif

   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, ip.t, shape_t, dshape_t);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z - ip.t, shape_l, dshape_l);

   int o = 0;
   for (int m = 0; m <= p; m++)
      for (int k = 0; k + m <= p; k++)
         for (int j = 0; j + k + m <= p; j++)
            for (int i = 0; i + j + k + m <= p; i++)
            {
               int l = p - i - j - k - m;
               divu(o++) = (dshape_x(i)*shape_l(l) -
                            shape_x(i)*dshape_l(l))*shape_y(j)*shape_z(k)*shape_t(m);
               divu(o++) = (dshape_y(j)*shape_l(l) -
                            shape_y(j)*dshape_l(l))*shape_x(i)*shape_z(k)*shape_t(m);
               divu(o++) = (dshape_z(k)*shape_l(l) -
                            shape_z(k)*dshape_l(l))*shape_x(i)*shape_y(j)*shape_t(m);
               divu(o++) = (dshape_t(m)*shape_l(l) -
                            shape_t(m)*dshape_l(l))*shape_x(i)*shape_y(j)*shape_z(k);
            }
   for (int l = 0; l <= p; l++)
      for (int j = 0; j + l<= p; j++)
         for (int i = 0; i + j + l <= p; i++)
         {
            int k = p - i - j - l;
            divu(o++) =
               (shape_x(i) + (ip.x - c)*dshape_x(i))*shape_y(j)*shape_z(l)*shape_t(k) +
               (shape_y(j) + (ip.y - c)*dshape_y(j))*shape_x(i)*shape_z(l)*shape_t(k) +
               (shape_z(l) + (ip.z - c)*dshape_z(l))*shape_x(i)*shape_y(j)*shape_t(k) +
               (shape_t(k) + (ip.t - c)*dshape_t(k))*shape_x(i)*shape_y(j)*shape_z(l);
         }

   Ti.Mult(divu, divshape);
}

void RT_PentatopeElement::ProjectDivSkew(const FiniteElement& fe,
                                         ElementTransformation& Trans, DenseMatrix& DivSkew)
{
   int dof = fe.GetDof();

   mfem_warning("RT_PentatopeElement::ProjectDivSkew(...) Implementation not tested!"); // TODO

#ifdef MFEM_THREAD_SAFE
   DenseMatrix Jinv(dim, dim);
#endif

   DivSkew.SetSize(dof,dof);
   DivSkew = 0.0;

   double n[4];
   Vector ni(n, 4);
   Vector vecF(4);

   DenseMatrix DivSkewshape(dof,4);
   DenseMatrix DivSkew_dFt(dof,4);
   for (int k = 0; k < dof; k++)
   {
      Trans.SetIntPoint(&Nodes.IntPoint(k));
      const DenseMatrix &J = Trans.Jacobian();
      CalcAdjugateTranspose(J, Jinv);

      fe.CalcDivSkewShape(Nodes.IntPoint(k), DivSkewshape);
      MultABt(DivSkewshape, J, DivSkew_dFt);
      DivSkew_dFt *= (1.0 / Trans.Weight());

      Jinv.Mult(nk + dof2nk[k] * dim,n);

      for (int j=0; j<dof; j++)
      {
         vecF(0) = DivSkew_dFt(j,0);
         vecF(1) = DivSkew_dFt(j,1);
         vecF(2) = DivSkew_dFt(j,2);
         vecF(3) = DivSkew_dFt(j,3);

         DivSkew(k, j) = vecF * ni;
      }
   }
}

const double ND_HexahedronElement::tk[18] =
{ 1.,0.,0.,  0.,1.,0.,  0.,0.,1., -1.,0.,0.,  0.,-1.,0.,  0.,0.,-1. };

ND_HexahedronElement::ND_HexahedronElement(const int p,
                                           const int cb_type, const int ob_type)
   : VectorTensorFiniteElement(3, 3*p*(p + 1)*(p + 1), p, cb_type, ob_type,
                               H_CURL, DofMapType::L2_DOF_MAP),
     dof2tk(dof)
{
   dof_map.SetSize(dof);

   const double *cp = poly1d.ClosedPoints(p, cb_type);
   const double *op = poly1d.OpenPoints(p - 1, ob_type);
   const int dof3 = dof/3;

#ifndef MFEM_THREAD_SAFE
   shape_cx.SetSize(p + 1);
   shape_ox.SetSize(p);
   shape_cy.SetSize(p + 1);
   shape_oy.SetSize(p);
   shape_cz.SetSize(p + 1);
   shape_oz.SetSize(p);
   dshape_cx.SetSize(p + 1);
   dshape_cy.SetSize(p + 1);
   dshape_cz.SetSize(p + 1);
#endif

   // edges
   int o = 0;
   for (int i = 0; i < p; i++)  // (0,1)
   {
      dof_map[0*dof3 + i + (0 + 0*(p + 1))*p] = o++;
   }
   for (int i = 0; i < p; i++)  // (1,2)
   {
      dof_map[1*dof3 + p + (i + 0*p)*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (3,2)
   {
      dof_map[0*dof3 + i + (p + 0*(p + 1))*p] = o++;
   }
   for (int i = 0; i < p; i++)  // (0,3)
   {
      dof_map[1*dof3 + 0 + (i + 0*p)*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (4,5)
   {
      dof_map[0*dof3 + i + (0 + p*(p + 1))*p] = o++;
   }
   for (int i = 0; i < p; i++)  // (5,6)
   {
      dof_map[1*dof3 + p + (i + p*p)*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (7,6)
   {
      dof_map[0*dof3 + i + (p + p*(p + 1))*p] = o++;
   }
   for (int i = 0; i < p; i++)  // (4,7)
   {
      dof_map[1*dof3 + 0 + (i + p*p)*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (0,4)
   {
      dof_map[2*dof3 + 0 + (0 + i*(p + 1))*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (1,5)
   {
      dof_map[2*dof3 + p + (0 + i*(p + 1))*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (2,6)
   {
      dof_map[2*dof3 + p + (p + i*(p + 1))*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (3,7)
   {
      dof_map[2*dof3 + 0 + (p + i*(p + 1))*(p + 1)] = o++;
   }

   // faces
   // (3,2,1,0) -- bottom
   for (int j = 1; j < p; j++) // x - components
      for (int i = 0; i < p; i++)
      {
         dof_map[0*dof3 + i + ((p - j) + 0*(p + 1))*p] = o++;
      }
   for (int j = 0; j < p; j++) // y - components
      for (int i = 1; i < p; i++)
      {
         dof_map[1*dof3 + i + ((p - 1 - j) + 0*p)*(p + 1)] = -1 - (o++);
      }
   // (0,1,5,4) -- front
   for (int k = 1; k < p; k++) // x - components
      for (int i = 0; i < p; i++)
      {
         dof_map[0*dof3 + i + (0 + k*(p + 1))*p] = o++;
      }
   for (int k = 0; k < p; k++) // z - components
      for (int i = 1; i < p; i++ )
      {
         dof_map[2*dof3 + i + (0 + k*(p + 1))*(p + 1)] = o++;
      }
   // (1,2,6,5) -- right
   for (int k = 1; k < p; k++) // y - components
      for (int j = 0; j < p; j++)
      {
         dof_map[1*dof3 + p + (j + k*p)*(p + 1)] = o++;
      }
   for (int k = 0; k < p; k++) // z - components
      for (int j = 1; j < p; j++)
      {
         dof_map[2*dof3 + p + (j + k*(p + 1))*(p + 1)] = o++;
      }
   // (2,3,7,6) -- back
   for (int k = 1; k < p; k++) // x - components
      for (int i = 0; i < p; i++)
      {
         dof_map[0*dof3 + (p - 1 - i) + (p + k*(p + 1))*p] = -1 - (o++);
      }
   for (int k = 0; k < p; k++) // z - components
      for (int i = 1; i < p; i++)
      {
         dof_map[2*dof3 + (p - i) + (p + k*(p + 1))*(p + 1)] = o++;
      }
   // (3,0,4,7) -- left
   for (int k = 1; k < p; k++) // y - components
      for (int j = 0; j < p; j++)
      {
         dof_map[1*dof3 + 0 + ((p - 1 - j) + k*p)*(p + 1)] = -1 - (o++);
      }
   for (int k = 0; k < p; k++) // z - components
      for (int j = 1; j < p; j++)
      {
         dof_map[2*dof3 + 0 + ((p - j) + k*(p + 1))*(p + 1)] = o++;
      }
   // (4,5,6,7) -- top
   for (int j = 1; j < p; j++) // x - components
      for (int i = 0; i < p; i++)
      {
         dof_map[0*dof3 + i + (j + p*(p + 1))*p] = o++;
      }
   for (int j = 0; j < p; j++) // y - components
      for (int i = 1; i < p; i++)
      {
         dof_map[1*dof3 + i + (j + p*p)*(p + 1)] = o++;
      }

   // interior
   // x-components
   for (int k = 1; k < p; k++)
      for (int j = 1; j < p; j++)
         for (int i = 0; i < p; i++)
         {
            dof_map[0*dof3 + i + (j + k*(p + 1))*p] = o++;
         }
   // y-components
   for (int k = 1; k < p; k++)
      for (int j = 0; j < p; j++)
         for (int i = 1; i < p; i++)
         {
            dof_map[1*dof3 + i + (j + k*p)*(p + 1)] = o++;
         }
   // z-components
   for (int k = 0; k < p; k++)
      for (int j = 1; j < p; j++)
         for (int i = 1; i < p; i++)
         {
            dof_map[2*dof3 + i + (j + k*(p + 1))*(p + 1)] = o++;
         }

   // set dof2tk and Nodes
   o = 0;
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i < p; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               dof2tk[idx = -1 - idx] = 3;
            }
            else
            {
               dof2tk[idx] = 0;
            }
            Nodes.IntPoint(idx).Set3(op[i], cp[j], cp[k]);
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j < p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               dof2tk[idx = -1 - idx] = 4;
            }
            else
            {
               dof2tk[idx] = 1;
            }
            Nodes.IntPoint(idx).Set3(cp[i], op[j], cp[k]);
         }
   // z-components
   for (int k = 0; k < p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx;
            if ((idx = dof_map[o++]) < 0)
            {
               dof2tk[idx = -1 - idx] = 5;
            }
            else
            {
               dof2tk[idx] = 2;
            }
            Nodes.IntPoint(idx).Set3(cp[i], cp[j], op[k]);
         }
}

void ND_HexahedronElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(p + 1), shape_ox(p), shape_cy(p + 1), shape_oy(p);
   Vector shape_cz(p + 1), shape_oz(p);
#endif

   cbasis1d.Eval(ip.x, shape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy);
   obasis1d.Eval(ip.y, shape_oy);
   cbasis1d.Eval(ip.z, shape_cz);
   obasis1d.Eval(ip.z, shape_oz);

   int o = 0;
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i < p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = s*shape_ox(i)*shape_cy(j)*shape_cz(k);
            shape(idx,1) = 0.;
            shape(idx,2) = 0.;
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j < p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = 0.;
            shape(idx,1) = s*shape_cx(i)*shape_oy(j)*shape_cz(k);
            shape(idx,2) = 0.;
         }
   // z-components
   for (int k = 0; k < p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            shape(idx,0) = 0.;
            shape(idx,1) = 0.;
            shape(idx,2) = s*shape_cx(i)*shape_cy(j)*shape_oz(k);
         }
}

void ND_HexahedronElement::CalcCurlShape(const IntegrationPoint &ip,
                                         DenseMatrix &curl_shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(p + 1), shape_ox(p), shape_cy(p + 1), shape_oy(p);
   Vector shape_cz(p + 1), shape_oz(p);
   Vector dshape_cx(p + 1), dshape_cy(p + 1), dshape_cz(p + 1);
#endif

   cbasis1d.Eval(ip.x, shape_cx, dshape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy, dshape_cy);
   obasis1d.Eval(ip.y, shape_oy);
   cbasis1d.Eval(ip.z, shape_cz, dshape_cz);
   obasis1d.Eval(ip.z, shape_oz);

   int o = 0;
   // x-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i < p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            curl_shape(idx,0) = 0.;
            curl_shape(idx,1) =  s*shape_ox(i)* shape_cy(j)*dshape_cz(k);
            curl_shape(idx,2) = -s*shape_ox(i)*dshape_cy(j)* shape_cz(k);
         }
   // y-components
   for (int k = 0; k <= p; k++)
      for (int j = 0; j < p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            curl_shape(idx,0) = -s* shape_cx(i)*shape_oy(j)*dshape_cz(k);
            curl_shape(idx,1) = 0.;
            curl_shape(idx,2) =  s*dshape_cx(i)*shape_oy(j)* shape_cz(k);
         }
   // z-components
   for (int k = 0; k < p; k++)
      for (int j = 0; j <= p; j++)
         for (int i = 0; i <= p; i++)
         {
            int idx, s;
            if ((idx = dof_map[o++]) < 0)
            {
               idx = -1 - idx, s = -1;
            }
            else
            {
               s = +1;
            }
            curl_shape(idx,0) =   s* shape_cx(i)*dshape_cy(j)*shape_oz(k);
            curl_shape(idx,1) =  -s*dshape_cx(i)* shape_cy(j)*shape_oz(k);
            curl_shape(idx,2) = 0.;
         }
}

const DofToQuad &VectorTensorFiniteElement::GetDofToQuad(
   const IntegrationRule &ir,
   DofToQuad::Mode mode) const
{
   MFEM_VERIFY(mode != DofToQuad::FULL, "invalid mode requested");

   return GetTensorDofToQuad(ir, mode, true);
}

const DofToQuad &VectorTensorFiniteElement::GetDofToQuadOpen(
   const IntegrationRule &ir,
   DofToQuad::Mode mode) const
{
   MFEM_VERIFY(mode != DofToQuad::FULL, "invalid mode requested");

   return GetTensorDofToQuad(ir, mode, false);
}

const DofToQuad &VectorTensorFiniteElement::GetTensorDofToQuad(
   const IntegrationRule &ir,
   DofToQuad::Mode mode,
   const bool closed) const
{
   MFEM_VERIFY(mode == DofToQuad::TENSOR, "invalid mode requested");

   for (int i = 0;
        i < (closed ? dof2quad_array.Size() : dof2quad_array_open.Size());
        i++)
   {
      const DofToQuad &d2q = closed ? *dof2quad_array[i] : *dof2quad_array_open[i];
      if (d2q.IntRule == &ir && d2q.mode == mode) { return d2q; }
   }

   DofToQuad *d2q = new DofToQuad;
   const int ndof = closed ? order + 1 : order;
   const int nqpt = (int)floor(pow(ir.GetNPoints(), 1.0/dim) + 0.5);
   d2q->FE = this;
   d2q->IntRule = &ir;
   d2q->mode = mode;
   d2q->ndof = ndof;
   d2q->nqpt = nqpt;
   d2q->B.SetSize(nqpt*ndof);
   d2q->Bt.SetSize(ndof*nqpt);
   d2q->G.SetSize(nqpt*ndof);
   d2q->Gt.SetSize(ndof*nqpt);
   Vector val(ndof), grad(ndof);
   for (int i = 0; i < nqpt; i++)
   {
      // The first 'nqpt' points in 'ir' have the same x-coordinates as those
      // of the 1D rule.

      if (closed)
      {
         cbasis1d.Eval(ir.IntPoint(i).x, val, grad);
      }
      else
      {
         obasis1d.Eval(ir.IntPoint(i).x, val, grad);
      }

      for (int j = 0; j < ndof; j++)
      {
         d2q->B[i+nqpt*j] = d2q->Bt[j+ndof*i] = val(j);
         d2q->G[i+nqpt*j] = d2q->Gt[j+ndof*i] = grad(j);
      }
   }

   if (closed)
   {
      dof2quad_array.Append(d2q);
   }
   else
   {
      dof2quad_array_open.Append(d2q);
   }

   return *d2q;
}

VectorTensorFiniteElement::~VectorTensorFiniteElement()
{
   for (int i = 0; i < dof2quad_array_open.Size(); i++)
   {
      delete dof2quad_array_open[i];
   }
}

const double ND_QuadrilateralElement::tk[8] =
{ 1.,0.,  0.,1., -1.,0., 0.,-1. };

ND_QuadrilateralElement::ND_QuadrilateralElement(const int p,
                                                 const int cb_type,
                                                 const int ob_type)
   : VectorTensorFiniteElement(2, 2*p*(p + 1), p, cb_type, ob_type,
                               H_CURL, DofMapType::L2_DOF_MAP),
     dof2tk(dof)
{
   dof_map.SetSize(dof);

   const double *cp = poly1d.ClosedPoints(p, cb_type);
   const double *op = poly1d.OpenPoints(p - 1, ob_type);
   const int dof2 = dof/2;

#ifndef MFEM_THREAD_SAFE
   shape_cx.SetSize(p + 1);
   shape_ox.SetSize(p);
   shape_cy.SetSize(p + 1);
   shape_oy.SetSize(p);
   dshape_cx.SetSize(p + 1);
   dshape_cy.SetSize(p + 1);
#endif

   // edges
   int o = 0;
   for (int i = 0; i < p; i++)  // (0,1)
   {
      dof_map[0*dof2 + i + 0*p] = o++;
   }
   for (int j = 0; j < p; j++)  // (1,2)
   {
      dof_map[1*dof2 + p + j*(p + 1)] = o++;
   }
   for (int i = 0; i < p; i++)  // (2,3)
   {
      dof_map[0*dof2 + (p - 1 - i) + p*p] = -1 - (o++);
   }
   for (int j = 0; j < p; j++)  // (3,0)
   {
      dof_map[1*dof2 + 0 + (p - 1 - j)*(p + 1)] = -1 - (o++);
   }

   // interior
   // x-components
   for (int j = 1; j < p; j++)
      for (int i = 0; i < p; i++)
      {
         dof_map[0*dof2 + i + j*p] = o++;
      }
   // y-components
   for (int j = 0; j < p; j++)
      for (int i = 1; i < p; i++)
      {
         dof_map[1*dof2 + i + j*(p + 1)] = o++;
      }

   // set dof2tk and Nodes
   o = 0;
   // x-components
   for (int j = 0; j <= p; j++)
      for (int i = 0; i < p; i++)
      {
         int idx;
         if ((idx = dof_map[o++]) < 0)
         {
            dof2tk[idx = -1 - idx] = 2;
         }
         else
         {
            dof2tk[idx] = 0;
         }
         Nodes.IntPoint(idx).Set2(op[i], cp[j]);
      }
   // y-components
   for (int j = 0; j < p; j++)
      for (int i = 0; i <= p; i++)
      {
         int idx;
         if ((idx = dof_map[o++]) < 0)
         {
            dof2tk[idx = -1 - idx] = 3;
         }
         else
         {
            dof2tk[idx] = 1;
         }
         Nodes.IntPoint(idx).Set2(cp[i], op[j]);
      }
}

void ND_QuadrilateralElement::CalcVShape(const IntegrationPoint &ip,
                                         DenseMatrix &shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(p + 1), shape_ox(p), shape_cy(p + 1), shape_oy(p);
#endif

   cbasis1d.Eval(ip.x, shape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy);
   obasis1d.Eval(ip.y, shape_oy);

   int o = 0;
   // x-components
   for (int j = 0; j <= p; j++)
      for (int i = 0; i < p; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         shape(idx,0) = s*shape_ox(i)*shape_cy(j);
         shape(idx,1) = 0.;
      }
   // y-components
   for (int j = 0; j < p; j++)
      for (int i = 0; i <= p; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         shape(idx,0) = 0.;
         shape(idx,1) = s*shape_cx(i)*shape_oy(j);
      }
}

void ND_QuadrilateralElement::CalcCurlShape(const IntegrationPoint &ip,
                                            DenseMatrix &curl_shape) const
{
   const int p = order;

#ifdef MFEM_THREAD_SAFE
   Vector shape_cx(p + 1), shape_ox(p), shape_cy(p + 1), shape_oy(p);
   Vector dshape_cx(p + 1), dshape_cy(p + 1);
#endif

   cbasis1d.Eval(ip.x, shape_cx, dshape_cx);
   obasis1d.Eval(ip.x, shape_ox);
   cbasis1d.Eval(ip.y, shape_cy, dshape_cy);
   obasis1d.Eval(ip.y, shape_oy);

   int o = 0;
   // x-components
   for (int j = 0; j <= p; j++)
      for (int i = 0; i < p; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         curl_shape(idx,0) = -s*shape_ox(i)*dshape_cy(j);
      }
   // y-components
   for (int j = 0; j < p; j++)
      for (int i = 0; i <= p; i++)
      {
         int idx, s;
         if ((idx = dof_map[o++]) < 0)
         {
            idx = -1 - idx, s = -1;
         }
         else
         {
            s = +1;
         }
         curl_shape(idx,0) =  s*dshape_cx(i)*shape_oy(j);
      }
}


const double ND_TetrahedronElement::tk[18] =
{ 1.,0.,0.,  0.,1.,0.,  0.,0.,1.,  -1.,1.,0.,  -1.,0.,1.,  0.,-1.,1. };

const double ND_TetrahedronElement::c = 1./4.;

ND_TetrahedronElement::ND_TetrahedronElement(const int p)
   : VectorFiniteElement(3, Geometry::TETRAHEDRON, p*(p + 2)*(p + 3)/2, p,
                         H_CURL, FunctionSpace::Pk), dof2tk(dof)
{
   const double *eop = poly1d.OpenPoints(p - 1);
   const double *fop = (p > 1) ? poly1d.OpenPoints(p - 2) : NULL;
   const double *iop = (p > 2) ? poly1d.OpenPoints(p - 3) : NULL;

   const int pm1 = p - 1, pm2 = p - 2, pm3 = p - 3;

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p);
   shape_y.SetSize(p);
   shape_z.SetSize(p);
   shape_l.SetSize(p);
   dshape_x.SetSize(p);
   dshape_y.SetSize(p);
   dshape_z.SetSize(p);
   dshape_l.SetSize(p);
   u.SetSize(dof, dim);
#else
   Vector shape_x(p), shape_y(p), shape_z(p), shape_l(p);
#endif

   int o = 0;
   // edges
   for (int i = 0; i < p; i++) // (0,1)
   {
      Nodes.IntPoint(o).Set3(eop[i], 0., 0.);
      dof2tk[o++] = 0;
   }
   for (int i = 0; i < p; i++) // (0,2)
   {
      Nodes.IntPoint(o).Set3(0., eop[i], 0.);
      dof2tk[o++] = 1;
   }
   for (int i = 0; i < p; i++) // (0,3)
   {
      Nodes.IntPoint(o).Set3(0., 0., eop[i]);
      dof2tk[o++] = 2;
   }
   for (int i = 0; i < p; i++) // (1,2)
   {
      Nodes.IntPoint(o).Set3(eop[pm1-i], eop[i], 0.);
      dof2tk[o++] = 3;
   }
   for (int i = 0; i < p; i++) // (1,3)
   {
      Nodes.IntPoint(o).Set3(eop[pm1-i], 0., eop[i]);
      dof2tk[o++] = 4;
   }
   for (int i = 0; i < p; i++) // (2,3)
   {
      Nodes.IntPoint(o).Set3(0., eop[pm1-i], eop[i]);
      dof2tk[o++] = 5;
   }

   // faces
   for (int j = 0; j <= pm2; j++)  // (1,2,3)
      for (int i = 0; i + j <= pm2; i++)
      {
         double w = fop[i] + fop[j] + fop[pm2-i-j];
         Nodes.IntPoint(o).Set3(fop[pm2-i-j]/w, fop[i]/w, fop[j]/w);
         dof2tk[o++] = 3;
         Nodes.IntPoint(o).Set3(fop[pm2-i-j]/w, fop[i]/w, fop[j]/w);
         dof2tk[o++] = 4;
      }
   for (int j = 0; j <= pm2; j++)  // (0,3,2)
      for (int i = 0; i + j <= pm2; i++)
      {
         double w = fop[i] + fop[j] + fop[pm2-i-j];
         Nodes.IntPoint(o).Set3(0., fop[j]/w, fop[i]/w);
         dof2tk[o++] = 2;
         Nodes.IntPoint(o).Set3(0., fop[j]/w, fop[i]/w);
         dof2tk[o++] = 1;
      }
   for (int j = 0; j <= pm2; j++)  // (0,1,3)
      for (int i = 0; i + j <= pm2; i++)
      {
         double w = fop[i] + fop[j] + fop[pm2-i-j];
         Nodes.IntPoint(o).Set3(fop[i]/w, 0., fop[j]/w);
         dof2tk[o++] = 0;
         Nodes.IntPoint(o).Set3(fop[i]/w, 0., fop[j]/w);
         dof2tk[o++] = 2;
      }
   for (int j = 0; j <= pm2; j++)  // (0,2,1)
      for (int i = 0; i + j <= pm2; i++)
      {
         double w = fop[i] + fop[j] + fop[pm2-i-j];
         Nodes.IntPoint(o).Set3(fop[j]/w, fop[i]/w, 0.);
         dof2tk[o++] = 1;
         Nodes.IntPoint(o).Set3(fop[j]/w, fop[i]/w, 0.);
         dof2tk[o++] = 0;
      }

   // interior
   for (int k = 0; k <= pm3; k++)
      for (int j = 0; j + k <= pm3; j++)
         for (int i = 0; i + j + k <= pm3; i++)
         {
            double w = iop[i] + iop[j] + iop[k] + iop[pm3-i-j-k];
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2tk[o++] = 0;
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2tk[o++] = 1;
            Nodes.IntPoint(o).Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            dof2tk[o++] = 2;
         }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      const double *tm = tk + 3*dof2tk[m];
      o = 0;

      poly1d.CalcBasis(pm1, ip.x, shape_x);
      poly1d.CalcBasis(pm1, ip.y, shape_y);
      poly1d.CalcBasis(pm1, ip.z, shape_z);
      poly1d.CalcBasis(pm1, 1. - ip.x - ip.y - ip.z, shape_l);

      for (int k = 0; k <= pm1; k++)
         for (int j = 0; j + k <= pm1; j++)
            for (int i = 0; i + j + k <= pm1; i++)
            {
               double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(pm1-i-j-k);
               T(o++, m) = s * tm[0];
               T(o++, m) = s * tm[1];
               T(o++, m) = s * tm[2];
            }
      for (int k = 0; k <= pm1; k++)
         for (int j = 0; j + k <= pm1; j++)
         {
            double s = shape_x(pm1-j-k)*shape_y(j)*shape_z(k);
            T(o++, m) = s*((ip.y - c)*tm[0] - (ip.x - c)*tm[1]);
            T(o++, m) = s*((ip.z - c)*tm[0] - (ip.x - c)*tm[2]);
         }
      for (int k = 0; k <= pm1; k++)
      {
         T(o++, m) =
            shape_y(pm1-k)*shape_z(k)*((ip.z - c)*tm[1] - (ip.y - c)*tm[2]);
      }
   }

   Ti.Factor(T);
   // mfem::out << "ND_TetrahedronElement(" << p << ") : "; Ti.TestInversion();
}

void ND_TetrahedronElement::CalcVShape(const IntegrationPoint &ip,
                                       DenseMatrix &shape) const
{
   const int pm1 = order - 1;

#ifdef MFEM_THREAD_SAFE
   const int p = order;
   Vector shape_x(p), shape_y(p), shape_z(p), shape_l(p);
   DenseMatrix u(dof, dim);
#endif

   poly1d.CalcBasis(pm1, ip.x, shape_x);
   poly1d.CalcBasis(pm1, ip.y, shape_y);
   poly1d.CalcBasis(pm1, ip.z, shape_z);
   poly1d.CalcBasis(pm1, 1. - ip.x - ip.y - ip.z, shape_l);

   int n = 0;
   for (int k = 0; k <= pm1; k++)
      for (int j = 0; j + k <= pm1; j++)
         for (int i = 0; i + j + k <= pm1; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(pm1-i-j-k);
            u(n,0) =  s;  u(n,1) = 0.;  u(n,2) = 0.;  n++;
            u(n,0) = 0.;  u(n,1) =  s;  u(n,2) = 0.;  n++;
            u(n,0) = 0.;  u(n,1) = 0.;  u(n,2) =  s;  n++;
         }
   for (int k = 0; k <= pm1; k++)
      for (int j = 0; j + k <= pm1; j++)
      {
         double s = shape_x(pm1-j-k)*shape_y(j)*shape_z(k);
         u(n,0) = s*(ip.y - c);  u(n,1) = -s*(ip.x - c);  u(n,2) =  0.;  n++;
         u(n,0) = s*(ip.z - c);  u(n,1) =  0.;  u(n,2) = -s*(ip.x - c);  n++;
      }
   for (int k = 0; k <= pm1; k++)
   {
      double s = shape_y(pm1-k)*shape_z(k);
      u(n,0) = 0.;  u(n,1) = s*(ip.z - c);  u(n,2) = -s*(ip.y - c);  n++;
   }

   Ti.Mult(u, shape);
}

void ND_TetrahedronElement::CalcCurlShape(const IntegrationPoint &ip,
                                          DenseMatrix &curl_shape) const
{
   const int pm1 = order - 1;

#ifdef MFEM_THREAD_SAFE
   const int p = order;
   Vector shape_x(p), shape_y(p), shape_z(p), shape_l(p);
   Vector dshape_x(p), dshape_y(p), dshape_z(p), dshape_l(p);
   DenseMatrix u(dof, dim);
#endif

   poly1d.CalcBasis(pm1, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(pm1, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(pm1, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(pm1, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l);

   int n = 0;
   for (int k = 0; k <= pm1; k++)
      for (int j = 0; j + k <= pm1; j++)
         for (int i = 0; i + j + k <= pm1; i++)
         {
            int l = pm1-i-j-k;
            const double dx = (dshape_x(i)*shape_l(l) -
                               shape_x(i)*dshape_l(l))*shape_y(j)*shape_z(k);
            const double dy = (dshape_y(j)*shape_l(l) -
                               shape_y(j)*dshape_l(l))*shape_x(i)*shape_z(k);
            const double dz = (dshape_z(k)*shape_l(l) -
                               shape_z(k)*dshape_l(l))*shape_x(i)*shape_y(j);

            u(n,0) =  0.;  u(n,1) =  dz;  u(n,2) = -dy;  n++;
            u(n,0) = -dz;  u(n,1) =  0.;  u(n,2) =  dx;  n++;
            u(n,0) =  dy;  u(n,1) = -dx;  u(n,2) =  0.;  n++;
         }
   for (int k = 0; k <= pm1; k++)
      for (int j = 0; j + k <= pm1; j++)
      {
         int i = pm1 - j - k;
         // s = shape_x(i)*shape_y(j)*shape_z(k);
         // curl of s*(ip.y - c, -(ip.x - c), 0):
         u(n,0) =  shape_x(i)*(ip.x - c)*shape_y(j)*dshape_z(k);
         u(n,1) =  shape_x(i)*shape_y(j)*(ip.y - c)*dshape_z(k);
         u(n,2) =
            -((dshape_x(i)*(ip.x - c) + shape_x(i))*shape_y(j)*shape_z(k) +
              (dshape_y(j)*(ip.y - c) + shape_y(j))*shape_x(i)*shape_z(k));
         n++;
         // curl of s*(ip.z - c, 0, -(ip.x - c)):
         u(n,0) = -shape_x(i)*(ip.x - c)*dshape_y(j)*shape_z(k);
         u(n,1) = (shape_x(i)*shape_y(j)*(dshape_z(k)*(ip.z - c) + shape_z(k)) +
                   (dshape_x(i)*(ip.x - c) + shape_x(i))*shape_y(j)*shape_z(k));
         u(n,2) = -shape_x(i)*dshape_y(j)*shape_z(k)*(ip.z - c);
         n++;
      }
   for (int k = 0; k <= pm1; k++)
   {
      int j = pm1 - k;
      // curl of shape_y(j)*shape_z(k)*(0, ip.z - c, -(ip.y - c)):
      u(n,0) = -((dshape_y(j)*(ip.y - c) + shape_y(j))*shape_z(k) +
                 shape_y(j)*(dshape_z(k)*(ip.z - c) + shape_z(k)));
      u(n,1) = 0.;
      u(n,2) = 0.;  n++;
   }

   Ti.Mult(u, curl_shape);
}


const double ND_TriangleElement::tk[8] =
{ 1.,0.,  -1.,1.,  0.,-1.,  0.,1. };

const double ND_TriangleElement::c = 1./3.;

ND_TriangleElement::ND_TriangleElement(const int p)
   : VectorFiniteElement(2, Geometry::TRIANGLE, p*(p + 2), p,
                         H_CURL, FunctionSpace::Pk),
     dof2tk(dof)
{
   const double *eop = poly1d.OpenPoints(p - 1);
   const double *iop = (p > 1) ? poly1d.OpenPoints(p - 2) : NULL;

   const int pm1 = p - 1, pm2 = p - 2;

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p);
   shape_y.SetSize(p);
   shape_l.SetSize(p);
   dshape_x.SetSize(p);
   dshape_y.SetSize(p);
   dshape_l.SetSize(p);
   u.SetSize(dof, dim);
   curlu.SetSize(dof);
#else
   Vector shape_x(p), shape_y(p), shape_l(p);
#endif

   int n = 0;
   // edges
   for (int i = 0; i < p; i++) // (0,1)
   {
      Nodes.IntPoint(n).Set2(eop[i], 0.);
      dof2tk[n++] = 0;
   }
   for (int i = 0; i < p; i++) // (1,2)
   {
      Nodes.IntPoint(n).Set2(eop[pm1-i], eop[i]);
      dof2tk[n++] = 1;
   }
   for (int i = 0; i < p; i++) // (2,0)
   {
      Nodes.IntPoint(n).Set2(0., eop[pm1-i]);
      dof2tk[n++] = 2;
   }

   // interior
   for (int j = 0; j <= pm2; j++)
      for (int i = 0; i + j <= pm2; i++)
      {
         double w = iop[i] + iop[j] + iop[pm2-i-j];
         Nodes.IntPoint(n).Set2(iop[i]/w, iop[j]/w);
         dof2tk[n++] = 0;
         Nodes.IntPoint(n).Set2(iop[i]/w, iop[j]/w);
         dof2tk[n++] = 3;
      }

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      const double *tm = tk + 2*dof2tk[m];
      n = 0;

      poly1d.CalcBasis(pm1, ip.x, shape_x);
      poly1d.CalcBasis(pm1, ip.y, shape_y);
      poly1d.CalcBasis(pm1, 1. - ip.x - ip.y, shape_l);

      for (int j = 0; j <= pm1; j++)
         for (int i = 0; i + j <= pm1; i++)
         {
            double s = shape_x(i)*shape_y(j)*shape_l(pm1-i-j);
            T(n++, m) = s * tm[0];
            T(n++, m) = s * tm[1];
         }
      for (int j = 0; j <= pm1; j++)
      {
         T(n++, m) =
            shape_x(pm1-j)*shape_y(j)*((ip.y - c)*tm[0] - (ip.x - c)*tm[1]);
      }
   }

   Ti.Factor(T);
   // mfem::out << "ND_TriangleElement(" << p << ") : "; Ti.TestInversion();
}

void ND_TriangleElement::CalcVShape(const IntegrationPoint &ip,
                                    DenseMatrix &shape) const
{
   const int pm1 = order - 1;

#ifdef MFEM_THREAD_SAFE
   const int p = order;
   Vector shape_x(p), shape_y(p), shape_l(p);
   DenseMatrix u(dof, dim);
#endif

   poly1d.CalcBasis(pm1, ip.x, shape_x);
   poly1d.CalcBasis(pm1, ip.y, shape_y);
   poly1d.CalcBasis(pm1, 1. - ip.x - ip.y, shape_l);

   int n = 0;
   for (int j = 0; j <= pm1; j++)
      for (int i = 0; i + j <= pm1; i++)
      {
         double s = shape_x(i)*shape_y(j)*shape_l(pm1-i-j);
         u(n,0) = s;  u(n,1) = 0;  n++;
         u(n,0) = 0;  u(n,1) = s;  n++;
      }
   for (int j = 0; j <= pm1; j++)
   {
      double s = shape_x(pm1-j)*shape_y(j);
      u(n,0) =  s*(ip.y - c);
      u(n,1) = -s*(ip.x - c);
      n++;
   }

   Ti.Mult(u, shape);
}

void ND_TriangleElement::CalcCurlShape(const IntegrationPoint &ip,
                                       DenseMatrix &curl_shape) const
{
   const int pm1 = order - 1;

#ifdef MFEM_THREAD_SAFE
   const int p = order;
   Vector shape_x(p), shape_y(p), shape_l(p);
   Vector dshape_x(p), dshape_y(p), dshape_l(p);
   Vector curlu(dof);
#endif

   poly1d.CalcBasis(pm1, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(pm1, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(pm1, 1. - ip.x - ip.y, shape_l, dshape_l);

   int n = 0;
   for (int j = 0; j <= pm1; j++)
      for (int i = 0; i + j <= pm1; i++)
      {
         int l = pm1-i-j;
         const double dx = (dshape_x(i)*shape_l(l) -
                            shape_x(i)*dshape_l(l)) * shape_y(j);
         const double dy = (dshape_y(j)*shape_l(l) -
                            shape_y(j)*dshape_l(l)) * shape_x(i);

         curlu(n++) = -dy;
         curlu(n++) =  dx;
      }

   for (int j = 0; j <= pm1; j++)
   {
      int i = pm1 - j;
      // curl of shape_x(i)*shape_y(j) * (ip.y - c, -(ip.x - c), 0):
      curlu(n++) = -((dshape_x(i)*(ip.x - c) + shape_x(i)) * shape_y(j) +
                     (dshape_y(j)*(ip.y - c) + shape_y(j)) * shape_x(i));
   }

   Vector curl2d(curl_shape.Data(),dof);
   Ti.Mult(curlu, curl2d);
}


const double ND_SegmentElement::tk[1] = { 1. };

ND_SegmentElement::ND_SegmentElement(const int p, const int ob_type)
   : VectorFiniteElement(1, Geometry::SEGMENT, p, p - 1,
                         H_CURL, FunctionSpace::Pk),
     obasis1d(poly1d.GetBasis(p - 1, VerifyOpen(ob_type))),
     dof2tk(dof)
{
   const double *op = poly1d.OpenPoints(p - 1, ob_type);

   // set dof2tk and Nodes
   for (int i = 0; i < p; i++)
   {
      dof2tk[i] = 0;
      Nodes.IntPoint(i).x = op[i];
   }
}

void ND_SegmentElement::CalcVShape(const IntegrationPoint &ip,
                                   DenseMatrix &shape) const
{
   Vector vshape(shape.Data(), dof);

   obasis1d.Eval(ip.x, vshape);
}

void NURBS1DFiniteElement::SetOrder() const
{
   order = kv[0]->GetOrder();
   dof = order + 1;

   weights.SetSize(dof);
   shape_x.SetSize(dof);
}

void NURBS1DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   kv[0]->CalcShape(shape, ijk[0], ip.x);

   double sum = 0.0;
   for (int i = 0; i <= order; i++)
   {
      sum += (shape(i) *= weights(i));
   }

   shape /= sum;
}

void NURBS1DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   Vector grad(dshape.Data(), dof);

   kv[0]->CalcShape (shape_x, ijk[0], ip.x);
   kv[0]->CalcDShape(grad,    ijk[0], ip.x);

   double sum = 0.0, dsum = 0.0;
   for (int i = 0; i <= order; i++)
   {
      sum  += (shape_x(i) *= weights(i));
      dsum += (   grad(i) *= weights(i));
   }

   sum = 1.0/sum;
   add(sum, grad, -dsum*sum*sum, shape_x, grad);
}

void NURBS1DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                        DenseMatrix &hessian) const
{
   Vector grad(dof);
   Vector hess(hessian.Data(), dof);

   kv[0]->CalcShape (shape_x,  ijk[0], ip.x);
   kv[0]->CalcDShape(grad,     ijk[0], ip.x);
   kv[0]->CalcD2Shape(hess,    ijk[0], ip.x);

   double sum = 0.0, dsum = 0.0, d2sum = 0.0;
   for (int i = 0; i <= order; i++)
   {
      sum   += (shape_x(i) *= weights(i));
      dsum  += (   grad(i) *= weights(i));
      d2sum += (   hess(i) *= weights(i));
   }

   sum = 1.0/sum;
   add(sum, hess, -2*dsum*sum*sum, grad, hess);
   add(1.0, hess, (-d2sum + 2*dsum*dsum*sum)*sum*sum, shape_x, hess);
}


void NURBS2DFiniteElement::SetOrder() const
{
   orders[0] = kv[0]->GetOrder();
   orders[1] = kv[1]->GetOrder();
   shape_x.SetSize(orders[0]+1);
   shape_y.SetSize(orders[1]+1);
   dshape_x.SetSize(orders[0]+1);
   dshape_y.SetSize(orders[1]+1);
   d2shape_x.SetSize(orders[0]+1);
   d2shape_y.SetSize(orders[1]+1);

   order = max(orders[0], orders[1]);
   dof = (orders[0] + 1)*(orders[1] + 1);
   u.SetSize(dof);
   du.SetSize(dof);
   weights.SetSize(dof);
}

void NURBS2DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   kv[0]->CalcShape(shape_x, ijk[0], ip.x);
   kv[1]->CalcShape(shape_y, ijk[1], ip.y);

   double sum = 0.0;
   for (int o = 0, j = 0; j <= orders[1]; j++)
   {
      const double sy = shape_y(j);
      for (int i = 0; i <= orders[0]; i++, o++)
      {
         sum += ( shape(o) = shape_x(i)*sy*weights(o) );
      }
   }

   shape /= sum;
}

void NURBS2DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   double sum, dsum[2];

   kv[0]->CalcShape ( shape_x, ijk[0], ip.x);
   kv[1]->CalcShape ( shape_y, ijk[1], ip.y);

   kv[0]->CalcDShape(dshape_x, ijk[0], ip.x);
   kv[1]->CalcDShape(dshape_y, ijk[1], ip.y);

   sum = dsum[0] = dsum[1] = 0.0;
   for (int o = 0, j = 0; j <= orders[1]; j++)
   {
      const double sy = shape_y(j), dsy = dshape_y(j);
      for (int i = 0; i <= orders[0]; i++, o++)
      {
         sum += ( u(o) = shape_x(i)*sy*weights(o) );

         dsum[0] += ( dshape(o,0) = dshape_x(i)*sy *weights(o) );
         dsum[1] += ( dshape(o,1) =  shape_x(i)*dsy*weights(o) );
      }
   }

   sum = 1.0/sum;
   dsum[0] *= sum*sum;
   dsum[1] *= sum*sum;

   for (int o = 0; o < dof; o++)
   {
      dshape(o,0) = dshape(o,0)*sum - u(o)*dsum[0];
      dshape(o,1) = dshape(o,1)*sum - u(o)*dsum[1];
   }
}

void NURBS2DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                        DenseMatrix &hessian) const
{
   double sum, dsum[2], d2sum[3];

   kv[0]->CalcShape ( shape_x, ijk[0], ip.x);
   kv[1]->CalcShape ( shape_y, ijk[1], ip.y);

   kv[0]->CalcDShape(dshape_x, ijk[0], ip.x);
   kv[1]->CalcDShape(dshape_y, ijk[1], ip.y);

   kv[0]->CalcD2Shape(d2shape_x, ijk[0], ip.x);
   kv[1]->CalcD2Shape(d2shape_y, ijk[1], ip.y);

   sum = dsum[0] = dsum[1] = 0.0;
   d2sum[0] = d2sum[1] = d2sum[2] = 0.0;
   for (int o = 0, j = 0; j <= orders[1]; j++)
   {
      const double sy = shape_y(j), dsy = dshape_y(j), d2sy = d2shape_y(j);
      for (int i = 0; i <= orders[0]; i++, o++)
      {
         const double sx = shape_x(i), dsx = dshape_x(i), d2sx = d2shape_x(i);
         sum += ( u(o) = sx*sy*weights(o) );

         dsum[0] += ( du(o,0) = dsx*sy*weights(o) );
         dsum[1] += ( du(o,1) = sx*dsy*weights(o) );

         d2sum[0] += ( hessian(o,0) = d2sx*sy*weights(o) );
         d2sum[1] += ( hessian(o,1) = dsx*dsy*weights(o) );
         d2sum[2] += ( hessian(o,2) = sx*d2sy*weights(o) );
      }
   }

   sum = 1.0/sum;
   dsum[0] *= sum;
   dsum[1] *= sum;

   d2sum[0] *= sum;
   d2sum[1] *= sum;
   d2sum[2] *= sum;

   for (int o = 0; o < dof; o++)
   {
      hessian(o,0) = hessian(o,0)*sum
                     - 2*du(o,0)*sum*dsum[0]
                     + u[o]*sum*(2*dsum[0]*dsum[0] - d2sum[0]);

      hessian(o,1) = hessian(o,1)*sum
                     - du(o,0)*sum*dsum[1]
                     - du(o,1)*sum*dsum[0]
                     + u[o]*sum*(2*dsum[0]*dsum[1] - d2sum[1]);

      hessian(o,2) = hessian(o,2)*sum
                     - 2*du(o,1)*sum*dsum[1]
                     + u[o]*sum*(2*dsum[1]*dsum[1] - d2sum[2]);
   }
}


void NURBS3DFiniteElement::SetOrder() const
{
   orders[0] = kv[0]->GetOrder();
   orders[1] = kv[1]->GetOrder();
   orders[2] = kv[2]->GetOrder();
   shape_x.SetSize(orders[0]+1);
   shape_y.SetSize(orders[1]+1);
   shape_z.SetSize(orders[2]+1);

   dshape_x.SetSize(orders[0]+1);
   dshape_y.SetSize(orders[1]+1);
   dshape_z.SetSize(orders[2]+1);

   d2shape_x.SetSize(orders[0]+1);
   d2shape_y.SetSize(orders[1]+1);
   d2shape_z.SetSize(orders[2]+1);

   order = max(max(orders[0], orders[1]), orders[2]);
   dof = (orders[0] + 1)*(orders[1] + 1)*(orders[2] + 1);
   u.SetSize(dof);
   du.SetSize(dof);
   weights.SetSize(dof);
}

void NURBS3DFiniteElement::CalcShape(const IntegrationPoint &ip,
                                     Vector &shape) const
{
   kv[0]->CalcShape(shape_x, ijk[0], ip.x);
   kv[1]->CalcShape(shape_y, ijk[1], ip.y);
   kv[2]->CalcShape(shape_z, ijk[2], ip.z);

   double sum = 0.0;
   for (int o = 0, k = 0; k <= orders[2]; k++)
   {
      const double sz = shape_z(k);
      for (int j = 0; j <= orders[1]; j++)
      {
         const double sy_sz = shape_y(j)*sz;
         for (int i = 0; i <= orders[0]; i++, o++)
         {
            sum += ( shape(o) = shape_x(i)*sy_sz*weights(o) );
         }
      }
   }

   shape /= sum;
}

void NURBS3DFiniteElement::CalcDShape(const IntegrationPoint &ip,
                                      DenseMatrix &dshape) const
{
   double sum, dsum[3];

   kv[0]->CalcShape ( shape_x, ijk[0], ip.x);
   kv[1]->CalcShape ( shape_y, ijk[1], ip.y);
   kv[2]->CalcShape ( shape_z, ijk[2], ip.z);

   kv[0]->CalcDShape(dshape_x, ijk[0], ip.x);
   kv[1]->CalcDShape(dshape_y, ijk[1], ip.y);
   kv[2]->CalcDShape(dshape_z, ijk[2], ip.z);

   sum = dsum[0] = dsum[1] = dsum[2] = 0.0;
   for (int o = 0, k = 0; k <= orders[2]; k++)
   {
      const double sz = shape_z(k), dsz = dshape_z(k);
      for (int j = 0; j <= orders[1]; j++)
      {
         const double  sy_sz  =  shape_y(j)* sz;
         const double dsy_sz  = dshape_y(j)* sz;
         const double  sy_dsz =  shape_y(j)*dsz;
         for (int i = 0; i <= orders[0]; i++, o++)
         {
            sum += ( u(o) = shape_x(i)*sy_sz*weights(o) );

            dsum[0] += ( dshape(o,0) = dshape_x(i)* sy_sz *weights(o) );
            dsum[1] += ( dshape(o,1) =  shape_x(i)*dsy_sz *weights(o) );
            dsum[2] += ( dshape(o,2) =  shape_x(i)* sy_dsz*weights(o) );
         }
      }
   }

   sum = 1.0/sum;
   dsum[0] *= sum*sum;
   dsum[1] *= sum*sum;
   dsum[2] *= sum*sum;

   for (int o = 0; o < dof; o++)
   {
      dshape(o,0) = dshape(o,0)*sum - u(o)*dsum[0];
      dshape(o,1) = dshape(o,1)*sum - u(o)*dsum[1];
      dshape(o,2) = dshape(o,2)*sum - u(o)*dsum[2];
   }
}

void NURBS3DFiniteElement::CalcHessian (const IntegrationPoint &ip,
                                        DenseMatrix &hessian) const
{
   double sum, dsum[3], d2sum[6];

   kv[0]->CalcShape ( shape_x, ijk[0], ip.x);
   kv[1]->CalcShape ( shape_y, ijk[1], ip.y);
   kv[2]->CalcShape ( shape_z, ijk[2], ip.z);

   kv[0]->CalcDShape(dshape_x, ijk[0], ip.x);
   kv[1]->CalcDShape(dshape_y, ijk[1], ip.y);
   kv[2]->CalcDShape(dshape_z, ijk[2], ip.z);

   kv[0]->CalcD2Shape(d2shape_x, ijk[0], ip.x);
   kv[1]->CalcD2Shape(d2shape_y, ijk[1], ip.y);
   kv[2]->CalcD2Shape(d2shape_z, ijk[2], ip.z);

   sum = dsum[0] = dsum[1] = dsum[2] = 0.0;
   d2sum[0] = d2sum[1] = d2sum[2] = d2sum[3] = d2sum[4] = d2sum[5] = 0.0;

   for (int o = 0, k = 0; k <= orders[2]; k++)
   {
      const double sz = shape_z(k), dsz = dshape_z(k), d2sz = d2shape_z(k);
      for (int j = 0; j <= orders[1]; j++)
      {
         const double sy = shape_y(j), dsy = dshape_y(j), d2sy = d2shape_y(j);
         for (int i = 0; i <= orders[0]; i++, o++)
         {
            const double sx = shape_x(i), dsx = dshape_x(i), d2sx = d2shape_x(i);
            sum += ( u(o) = sx*sy*sz*weights(o) );

            dsum[0] += ( du(o,0) = dsx*sy*sz*weights(o) );
            dsum[1] += ( du(o,1) = sx*dsy*sz*weights(o) );
            dsum[2] += ( du(o,2) = sx*sy*dsz*weights(o) );

            d2sum[0] += ( hessian(o,0) = d2sx*sy*sz*weights(o) );
            d2sum[1] += ( hessian(o,1) = dsx*dsy*sz*weights(o) );
            d2sum[2] += ( hessian(o,2) = dsx*sy*dsz*weights(o) );

            d2sum[3] += ( hessian(o,3) = sx*dsy*dsz*weights(o) );

            d2sum[4] += ( hessian(o,4) = sx*sy*d2sz*weights(o) );
            d2sum[5] += ( hessian(o,5) = sx*d2sy*sz*weights(o) );
         }
      }
   }

   sum = 1.0/sum;
   dsum[0] *= sum;
   dsum[1] *= sum;
   dsum[2] *= sum;

   d2sum[0] *= sum;
   d2sum[1] *= sum;
   d2sum[2] *= sum;

   d2sum[3] *= sum;
   d2sum[4] *= sum;
   d2sum[5] *= sum;

   for (int o = 0; o < dof; o++)
   {
      hessian(o,0) = hessian(o,0)*sum
                     - 2*du(o,0)*sum*dsum[0]
                     + u[o]*sum*(2*dsum[0]*dsum[0] - d2sum[0]);

      hessian(o,1) = hessian(o,1)*sum
                     - du(o,0)*sum*dsum[1]
                     - du(o,1)*sum*dsum[0]
                     + u[o]*sum*(2*dsum[0]*dsum[1] - d2sum[1]);

      hessian(o,2) = hessian(o,2)*sum
                     - du(o,0)*sum*dsum[2]
                     - du(o,2)*sum*dsum[0]
                     + u[o]*sum*(2*dsum[0]*dsum[2] - d2sum[2]);

      hessian(o,3) = hessian(o,3)*sum
                     - du(o,1)*sum*dsum[2]
                     - du(o,2)*sum*dsum[1]
                     + u[o]*sum*(2*dsum[1]*dsum[2] - d2sum[3]);

      hessian(o,4) = hessian(o,4)*sum
                     - 2*du(o,2)*sum*dsum[2]
                     + u[o]*sum*(2*dsum[2]*dsum[2] - d2sum[4]);

      hessian(o,5) = hessian(o,5)*sum
                     - 2*du(o,1)*sum*dsum[1]
                     + u[o]*sum*(2*dsum[1]*dsum[1] - d2sum[5]);

   }
}

// Global object definitions

// Object declared in mesh/triangle.hpp.
// Defined here to ensure it is constructed before 'Geometries'.
Linear2DFiniteElement TriangleFE;

// Object declared in mesh/tetrahedron.hpp.
// Defined here to ensure it is constructed before 'Geometries'.
Linear3DFiniteElement TetrahedronFE;

// Object declared in mesh/wedge.hpp.
// Defined here to ensure it is constructed after 'poly1d' and before
// 'Geometries'.
// TODO: define as thread_local to prevent race conditions in GLVis, because
// there is no "LinearWedgeFiniteElement" and WedgeFE is in turn used from two
// different threads for different things in GLVis. We also don't want to turn
// MFEM_THREAD_SAFE on globally. (See PR #731)
H1_WedgeElement WedgeFE(1);

// Object declared in geom.hpp.
// Construct 'Geometries' after 'TriangleFE', 'TetrahedronFE', and 'WedgeFE'.
Geometry Geometries;

}
