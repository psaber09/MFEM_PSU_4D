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

// Implementation of class LinearForm

#include "fem.hpp"

namespace mfem
{

LinearForm::LinearForm(FiniteElementSpace *f, LinearForm *lf)
   : Vector(f->GetVSize())
{
   // Linear forms are stored on the device
   UseDevice(true);

   fes = f;
   extern_lfs = 1;

   // Copy the pointers to the integrators
   dlfi = lf->dlfi;

   dlfi_delta = lf->dlfi_delta;

   blfi = lf->blfi;

   flfi = lf->flfi;
   flfi_marker = lf->flfi_marker;
}

void LinearForm::AddDomainIntegrator(LinearFormIntegrator *lfi)
{
   DeltaLFIntegrator *maybe_delta =
      dynamic_cast<DeltaLFIntegrator *>(lfi);
   if (!maybe_delta || !maybe_delta->IsDelta())
   {
      dlfi.Append(lfi);
   }
   else
   {
      dlfi_delta.Append(maybe_delta);
   }
}

void LinearForm::AddBoundaryIntegrator (LinearFormIntegrator * lfi)
{
   blfi.Append (lfi);
   blfi_marker.Append(NULL); // NULL -> all attributes are active
}

void LinearForm::AddBoundaryIntegrator (LinearFormIntegrator * lfi,
                                        Array<int> &bdr_attr_marker)
{
   blfi.Append (lfi);
   blfi_marker.Append(&bdr_attr_marker);
}

void LinearForm::AddInteriorFaceIntegrator (LinearFormIntegrator * lfi)
{
   iflfi.Append (lfi);
}

void LinearForm::AddBdrFaceIntegrator (LinearFormIntegrator * lfi)
{
   flfi.Append(lfi);
   flfi_marker.Append(NULL); // NULL -> all attributes are active
}

void LinearForm::AddBdrFaceIntegrator(LinearFormIntegrator *lfi,
                                      Array<int> &bdr_attr_marker)
{
   flfi.Append(lfi);
   flfi_marker.Append(&bdr_attr_marker);
}

void LinearForm::Assemble()
{
   Array<int> vdofs;
   ElementTransformation *eltrans;
   Vector elemvect;

   int i;

   Vector::operator=(0.0);

   // The above operation is executed on device because of UseDevice().
   // The first use of AddElementVector() below will move it back to host
   // because both 'vdofs' and 'elemvect' are on host.

   if (dlfi.Size())
   {
      for (i = 0; i < fes -> GetNE(); i++)
      {
         fes -> GetElementVDofs (i, vdofs);
         eltrans = fes -> GetElementTransformation (i);
         for (int k=0; k < dlfi.Size(); k++)
         {
            dlfi[k]->AssembleRHSElementVect(*fes->GetFE(i), *eltrans, elemvect);
            AddElementVector (vdofs, elemvect);
         }
      }
   }
   AssembleDelta();

   if (blfi.Size())
   {
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      bdr_attr_marker = 0;
      for (int k = 0; k < blfi.Size(); k++)
      {
         if (blfi_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *blfi_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (i = 0; i < fes -> GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }
         fes -> GetBdrElementVDofs (i, vdofs);
         eltrans = fes -> GetBdrElementTransformation (i);
         for (int k=0; k < blfi.Size(); k++)
         {
            if (blfi_marker[k] &&
                (*blfi_marker[k])[bdr_attr-1] == 0) { continue; }

            blfi[k]->AssembleRHSElementVect(*fes->GetBE(i), *eltrans, elemvect);

            AddElementVector (vdofs, elemvect);
         }
      }
   }
   if (flfi.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      bdr_attr_marker = 0;
      for (int k = 0; k < flfi.Size(); k++)
      {
         if (flfi_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *flfi_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary face integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (i = 0; i < mesh->GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

         tr = mesh->GetBdrFaceTransformations(i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < flfi.Size(); k++)
            {
               if (flfi_marker[k] &&
                   (*flfi_marker[k])[bdr_attr-1] == 0) { continue; }

               flfi[k] -> AssembleRHSElementVect (*fes->GetFE(tr -> Elem1No),
                                                  *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }

   if (iflfi.Size())
   {
      FaceElementTransformations *tr;
      Array<int> vdofs2;
      Mesh *mesh = fes -> GetMesh();

      int nfaces = mesh->GetNumFaces();
      for (i = 0; i < nfaces; i++)
      {
         tr = mesh -> GetInteriorFaceTransformations (i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            fes -> GetElementVDofs (tr -> Elem2No, vdofs2);
            vdofs.Append (vdofs2);
            for (int k = 0; k < iflfi.Size(); k++)
            {
               iflfi[k] -> AssembleRHSElementVect (*fes -> GetFE (tr -> Elem1No),
                                                   *fes -> GetFE (tr -> Elem2No),
                                                   *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
}

void LinearForm::Update(FiniteElementSpace *f, Vector &v, int v_offset)
{
   fes = f;
   NewMemoryAndSize(Memory<double>(v.GetMemory(), v_offset, f->GetVSize()),
                    f->GetVSize(), false);
   ResetDeltaLocations();
}

void LinearForm::MakeRef(FiniteElementSpace *f, Vector &v, int v_offset)
{
   MFEM_ASSERT(v.Size() >= v_offset + f->GetVSize(), "");
   fes = f;
   v.UseDevice(true);
   this->Vector::MakeRef(v, v_offset, fes->GetVSize());
}

void LinearForm::AssembleDelta()
{
   if (dlfi_delta.Size() == 0) { return; }

   if (!HaveDeltaLocations())
   {
      int sdim = fes->GetMesh()->SpaceDimension();
      Vector center;
      DenseMatrix centers(sdim, dlfi_delta.Size());
      for (int i = 0; i < centers.Width(); i++)
      {
         centers.GetColumnReference(i, center);
         dlfi_delta[i]->GetDeltaCenter(center);
         MFEM_VERIFY(center.Size() == sdim,
                     "Point dim " << center.Size() <<
                     " does not match space dim " << sdim);
      }
      fes->GetMesh()->FindPoints(centers, dlfi_delta_elem_id, dlfi_delta_ip);
   }

   Array<int> vdofs;
   Vector elemvect;
   for (int i = 0; i < dlfi_delta.Size(); i++)
   {
      int elem_id = dlfi_delta_elem_id[i];
      // The delta center may be outside of this sub-domain, or
      // (Par)Mesh::FindPoints() failed to find this point:
      if (elem_id < 0) { continue; }

      const IntegrationPoint &ip = dlfi_delta_ip[i];
      ElementTransformation &Trans = *fes->GetElementTransformation(elem_id);
      Trans.SetIntPoint(&ip);

      fes->GetElementVDofs(elem_id, vdofs);
      dlfi_delta[i]->AssembleDeltaElementVect(*fes->GetFE(elem_id), Trans,
                                              elemvect);
      AddElementVector(vdofs, elemvect);
   }
}

LinearForm & LinearForm::operator=(double value)
{
   Vector::operator=(value);
   return *this;
}

LinearForm & LinearForm::operator=(const Vector &v)
{
   MFEM_ASSERT(fes && v.Size() == fes->GetVSize(), "");
   Vector::operator=(v);
   return *this;
}

LinearForm::~LinearForm()
{
   if (!extern_lfs)
   {
      int k;
      for (k=0; k < dlfi_delta.Size(); k++) { delete dlfi_delta[k]; }
      for (k=0; k < dlfi.Size(); k++) { delete dlfi[k]; }
      for (k=0; k < blfi.Size(); k++) { delete blfi[k]; }
      for (k=0; k < flfi.Size(); k++) { delete flfi[k]; }
      for (k=0; k < iflfi.Size(); k++) { delete iflfi[k]; }
   }
}

}
