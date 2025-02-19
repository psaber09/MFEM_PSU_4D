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

#include "fem.hpp"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#ifdef _WIN32
#define snprintf _snprintf_s
#endif

namespace mfem
{

using namespace std;

int FiniteElementCollection::HasFaceDofs(Geometry::Type geom, int p) const
{
   switch (geom)
   {
      case Geometry::TETRAHEDRON:
         return GetNumDof(Geometry::TRIANGLE, p);
      case Geometry::CUBE:
         return GetNumDof(Geometry::SQUARE, p);
      case Geometry::PRISM:
         return max(GetNumDof(Geometry::TRIANGLE, p),
                    GetNumDof(Geometry::SQUARE, p));
      case Geometry::PENTATOPE:   return DofForGeometry (Geometry::TETRAHEDRON);
      case Geometry::TESSERACT:   return DofForGeometry (Geometry::CUBE);
      default:
         MFEM_ABORT("unknown geometry type");
   }
   return 0;
}

int FiniteElementCollection::HasPlanarDofs(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::PENTATOPE:   return DofForGeometry (Geometry::TRIANGLE);
      case Geometry::TESSERACT:   return DofForGeometry (Geometry::SQUARE);
      default:
         mfem_error ("FiniteElementCollection::HasPlanarDofs:"
                     " unknown geometry type.");
   }
   return 0;
}

FiniteElementCollection *FiniteElementCollection::GetTraceCollection() const
{
   MFEM_ABORT("this method is not implemented in this derived class!");
   return NULL;
}

FiniteElementCollection *FiniteElementCollection::New(const char *name)
{
   FiniteElementCollection *fec = NULL;

   if (!strcmp(name, "Linear"))
   {
      fec = new LinearFECollection;
   }
   else if (!strcmp(name, "Quadratic"))
   {
      fec = new QuadraticFECollection;
   }
   else if (!strcmp(name, "QuadraticPos"))
   {
      fec = new QuadraticPosFECollection;
   }
   else if (!strcmp(name, "Cubic"))
   {
      fec = new CubicFECollection;
   }
   else if (!strcmp(name, "Const3D"))
   {
      fec = new Const3DFECollection;
   }
   else if (!strcmp(name, "Const2D"))
   {
      fec = new Const2DFECollection;
   }
   else if (!strcmp(name, "LinearDiscont2D"))
   {
      fec = new LinearDiscont2DFECollection;
   }
   else if (!strcmp(name, "GaussLinearDiscont2D"))
   {
      fec = new GaussLinearDiscont2DFECollection;
   }
   else if (!strcmp(name, "P1OnQuad"))
   {
      fec = new P1OnQuadFECollection;
   }
   else if (!strcmp(name, "QuadraticDiscont2D"))
   {
      fec = new QuadraticDiscont2DFECollection;
   }
   else if (!strcmp(name, "QuadraticPosDiscont2D"))
   {
      fec = new QuadraticPosDiscont2DFECollection;
   }
   else if (!strcmp(name, "GaussQuadraticDiscont2D"))
   {
      fec = new GaussQuadraticDiscont2DFECollection;
   }
   else if (!strcmp(name, "CubicDiscont2D"))
   {
      fec = new CubicDiscont2DFECollection;
   }
   else if (!strcmp(name, "LinearDiscont3D"))
   {
      fec = new LinearDiscont3DFECollection;
   }
   else if (!strcmp(name, "QuadraticDiscont3D"))
   {
      fec = new QuadraticDiscont3DFECollection;
   }
   else if (!strcmp(name, "LinearNonConf3D"))
   {
      fec = new LinearNonConf3DFECollection;
   }
   else if (!strcmp(name, "CrouzeixRaviart"))
   {
      fec = new CrouzeixRaviartFECollection;
   }
   else if (!strcmp(name, "ND1_3D"))
   {
      fec = new ND1_3DFECollection;
   }
   else if (!strcmp(name, "RT0_2D"))
   {
      fec = new RT0_2DFECollection;
   }
   else if (!strcmp(name, "RT1_2D"))
   {
      fec = new RT1_2DFECollection;
   }
   else if (!strcmp(name, "RT2_2D"))
   {
      fec = new RT2_2DFECollection;
   }
   else if (!strcmp(name, "RT0_3D"))
   {
      fec = new RT0_3DFECollection;
   }
   else if (!strcmp(name, "RT1_3D"))
   {
      fec = new RT1_3DFECollection;
   }
   else if (!strncmp(name, "H1_Trace_", 9))
   {
      fec = new H1_Trace_FECollection(atoi(name + 13), atoi(name + 9));
   }
   else if (!strncmp(name, "H1_Trace@", 9))
   {
      fec = new H1_Trace_FECollection(atoi(name + 15), atoi(name + 11),
                                      BasisType::GetType(name[9]));
   }
   else if (!strncmp(name, "H1_", 3))
   {
      fec = new H1_FECollection(atoi(name + 7), atoi(name + 3));
   }
   else if (!strncmp(name, "H1Pos_Trace_", 12))
   {
      fec = new H1_Trace_FECollection(atoi(name + 16), atoi(name + 12),
                                      BasisType::Positive);
   }
   else if (!strncmp(name, "H1Pos_", 6))
   {
      fec = new H1Pos_FECollection(atoi(name + 10), atoi(name + 6));
   }
   else if (!strncmp(name, "H1Ser_", 6))
   {
      fec = new H1Ser_FECollection(atoi(name + 10), atoi(name + 6));
   }
   else if (!strncmp(name, "H1@", 3))
   {
      fec = new H1_FECollection(atoi(name + 9), atoi(name + 5),
                                BasisType::GetType(name[3]));
   }
   else if (!strncmp(name, "L2_T", 4))
      fec = new L2_FECollection(atoi(name + 10), atoi(name + 6),
                                atoi(name + 4));
   else if (!strncmp(name, "L2_", 3))
   {
      fec = new L2_FECollection(atoi(name + 7), atoi(name + 3));
   }
   else if (!strncmp(name, "L2Int_T", 7))
   {
      fec = new L2_FECollection(atoi(name + 13), atoi(name + 9),
                                atoi(name + 7), FiniteElement::INTEGRAL);
   }
   else if (!strncmp(name, "L2Int_", 6))
   {
      fec = new L2_FECollection(atoi(name + 10), atoi(name + 6),
                                BasisType::GaussLegendre,
                                FiniteElement::INTEGRAL);
   }
   else if (!strncmp(name, "RT_Trace_", 9))
   {
      fec = new RT_Trace_FECollection(atoi(name + 13), atoi(name + 9));
   }
   else if (!strncmp(name, "RT_ValTrace_", 12))
   {
      fec = new RT_Trace_FECollection(atoi(name + 16), atoi(name + 12),
                                      FiniteElement::VALUE);
   }
   else if (!strncmp(name, "RT_Trace@", 9))
   {
      fec = new RT_Trace_FECollection(atoi(name + 15), atoi(name + 11),
                                      FiniteElement::INTEGRAL,
                                      BasisType::GetType(name[9]));
   }
   else if (!strncmp(name, "RT_ValTrace@", 12))
   {
      fec = new RT_Trace_FECollection(atoi(name + 18), atoi(name + 14),
                                      FiniteElement::VALUE,
                                      BasisType::GetType(name[12]));
   }
   else if (!strncmp(name, "DG_Iface_", 9))
   {
      fec = new DG_Interface_FECollection(atoi(name + 13), atoi(name + 9));
   }
   else if (!strncmp(name, "DG_Iface@", 9))
   {
      fec = new DG_Interface_FECollection(atoi(name + 15), atoi(name + 11),
                                          FiniteElement::VALUE,
                                          BasisType::GetType(name[9]));
   }
   else if (!strncmp(name, "DG_IntIface_", 12))
   {
      fec = new DG_Interface_FECollection(atoi(name + 16), atoi(name + 12),
                                          FiniteElement::INTEGRAL);
   }
   else if (!strncmp(name, "DG_IntIface@", 12))
   {
      fec = new DG_Interface_FECollection(atoi(name + 18), atoi(name + 14),
                                          FiniteElement::INTEGRAL,
                                          BasisType::GetType(name[12]));
   }
   else if (!strncmp(name, "RT_", 3))
   {
      fec = new RT_FECollection(atoi(name + 7), atoi(name + 3));
   }
   else if (!strncmp(name, "RT@", 3))
   {
      fec = new RT_FECollection(atoi(name + 10), atoi(name + 6),
                                BasisType::GetType(name[3]),
                                BasisType::GetType(name[4]));
   }
   else if (!strncmp(name, "ND_Trace_", 9))
   {
      fec = new ND_Trace_FECollection(atoi(name + 13), atoi(name + 9));
   }
   else if (!strncmp(name, "ND_Trace@", 9))
   {
      fec = new ND_Trace_FECollection(atoi(name + 16), atoi(name + 12),
                                      BasisType::GetType(name[9]),
                                      BasisType::GetType(name[10]));
   }
   else if (!strncmp(name, "ND_", 3))
   {
      fec = new ND_FECollection(atoi(name + 7), atoi(name + 3));
   }
   else if (!strncmp(name, "ND@", 3))
   {
      fec = new ND_FECollection(atoi(name + 10), atoi(name + 6),
                                BasisType::GetType(name[3]),
                                BasisType::GetType(name[4]));
   }
   else if (!strncmp(name, "Local_", 6))
   {
      fec = new Local_FECollection(name + 6);
   }
   else if (!strncmp(name, "NURBS", 5))
   {
      if (name[5] != '\0')
      {
         // "NURBS" + "number" --> fixed order nurbs collection
         fec = new NURBSFECollection(atoi(name + 5));
      }
      else
      {
         // "NURBS" --> variable order nurbs collection
         fec = new NURBSFECollection();
      }
   }
   else
   {
      MFEM_ABORT("unknown FiniteElementCollection: " << name);
   }
   MFEM_VERIFY(!strcmp(fec->Name(), name), "input name: \"" << name
               << "\" does not match the created collection name: \""
               << fec->Name() << '"');

   return fec;
}

FiniteElementCollection *FiniteElementCollection::Clone(int p) const
{
   // default implementation for collections that don't care about variable p
   MFEM_ABORT("Collection " << Name() << " does not support variable orders.");
   (void) p;
   return NULL;
}

void FiniteElementCollection::InitVarOrder(int p) const
{
   if (p >= var_orders.Size())
   {
      var_orders.SetSize(p+1, NULL);
   }
   var_orders[p] = Clone(p);
}

FiniteElementCollection::~FiniteElementCollection()
{
   for (int i = 0; i < var_orders.Size(); i++)
   {
      delete var_orders[i];
   }
}

template <Geometry::Type geom>
inline void FiniteElementCollection::GetNVE(int &nv, int &ne)
{
   typedef typename Geometry::Constants<geom> g_consts;

   nv = g_consts::NumVert;
   ne = g_consts::NumEdges;
}

template <Geometry::Type geom, typename v_t>
inline void FiniteElementCollection::
GetEdge(int &nv, v_t &v, int &ne, int &e, int &eo, const int edge_info)
{
   typedef typename Geometry::Constants<Geometry::SEGMENT> e_consts;
   typedef typename Geometry::Constants<geom> g_consts;

   nv = e_consts::NumVert;
   ne = 1;
   e = edge_info/64;
   eo = edge_info%64;
   MFEM_ASSERT(0 <= e && e < g_consts::NumEdges, "");
   MFEM_ASSERT(0 <= eo && eo < e_consts::NumOrient, "");
   v[0] = e_consts::Orient[eo][0];
   v[1] = e_consts::Orient[eo][1];
   v[0] = g_consts::Edges[e][v[0]];
   v[1] = g_consts::Edges[e][v[1]];
}

template <Geometry::Type geom, Geometry::Type f_geom,
          typename v_t, typename e_t, typename eo_t>
inline void FiniteElementCollection::
GetFace(int &nv, v_t &v, int &ne, e_t &e, eo_t &eo,
        int &nf, int &f, Geometry::Type &fg, int &fo, const int face_info)
{
   typedef typename Geometry::Constants<  geom> g_consts;
   typedef typename Geometry::Constants<f_geom> f_consts;

   nv = f_consts::NumVert;
   nf = 1;
   f = face_info/64;
   fg = f_geom;
   fo = face_info%64;
   MFEM_ASSERT(0 <= f && f < g_consts::NumFaces, "");
   MFEM_ASSERT(0 <= fo && fo < f_consts::NumOrient, "");
   for (int i = 0; i < f_consts::NumVert; i++)
   {
      v[i] = f_consts::Orient[fo][i];
      v[i] = g_consts::FaceVert[f][v[i]];
   }
   ne = f_consts::NumEdges;
   for (int i = 0; i < f_consts::NumEdges; i++)
   {
      int v0 = v[f_consts::Edges[i][0]];
      int v1 = v[f_consts::Edges[i][1]];
      int eor = 0;
      if (v0 > v1) { swap(v0, v1); eor = 1; }
      for (int j = g_consts::VertToVert::I[v0]; true; j++)
      {
         MFEM_ASSERT(j < g_consts::VertToVert::I[v0+1],
                     "internal error, edge not found");
         if (v1 == g_consts::VertToVert::J[j][0])
         {
            int en = g_consts::VertToVert::J[j][1];
            if (en < 0)
            {
               en = -1-en;
               eor = 1-eor;
            }
            e[i] = en;
            eo[i] = eor;
            break;
         }
      }
   }
}

//template <Geometry::Type geom, Geometry::Type>// f_geom,
        //typename v_t, typename e_t, typename eo_t, typename f_t, typename fo_t>
inline void FiniteElementCollection::
GetFacet()//int &nv, v_t &v, int &ne, e_t &e, eo_t &eo,
        //int &nf, f_t &f, Geometry::Type &fg, fo_t &fo, int &nfacet, int &facet, Geometry::Type &facetg, int &faceto,
         //const int facet_info)
{
//   typedef typename Geometry::Constants<  geom> g_consts;
//   typedef typename Geometry::Constants<f_geom> f_consts;
//
//   nv = f_consts::NumVert;
//   nf = 1;
//   f = face_info/64;
//   fg = f_geom;
//   fo = face_info%64;
//   MFEM_ASSERT(0 <= f && f < g_consts::NumFaces, "");
//   MFEM_ASSERT(0 <= fo && fo < f_consts::NumOrient, "");
//   for (int i = 0; i < f_consts::NumVert; i++)
//   {
//      v[i] = f_consts::Orient[fo][i];
//      v[i] = g_consts::FaceVert[f][v[i]];
//   }
//   ne = f_consts::NumEdges;
//   for (int i = 0; i < f_consts::NumEdges; i++)
//   {
//      int v0 = v[f_consts::Edges[i][0]];
//      int v1 = v[f_consts::Edges[i][1]];
//      int eor = 0;
//      if (v0 > v1) { swap(v0, v1); eor = 1; }
//      for (int j = g_consts::VertToVert::I[v0]; true; j++)
//      {
//         MFEM_ASSERT(j < g_consts::VertToVert::I[v0+1],
//                     "internal error, edge not found");
//         if (v1 == g_consts::VertToVert::J[j][0])
//         {
//            int en = g_consts::VertToVert::J[j][1];
//            if (en < 0)
//            {
//               en = -1-en;
//               eor = 1-eor;
//            }
//            e[i] = en;
//            eo[i] = eor;
//            break;
//         }
//      }
//   }
}

template <Geometry::Type geom>
inline void FiniteElementCollection::GetNF(int &nf)
{
   typedef typename Geometry::Constants<geom> g_consts;

   nf = g_consts::NumFaces;
}

void FiniteElementCollection::SubDofOrder(Geometry::Type Geom, int SDim,
                                          int Info,
                                          Array<int> &dofs) const
{
   // Info = 64 * SubIndex + SubOrientation
   MFEM_ASSERT(0 <= Geom && Geom < Geometry::NumGeom,
               "invalid Geom = " << Geom);
   MFEM_ASSERT(0 <= SDim && SDim <= Geometry::Dimension[Geom],
               "invalid SDim = " << SDim <<
               " for Geom = " << Geometry::Name[Geom]);

   const int nvd = DofForGeometry(Geometry::POINT);
   if (SDim == 0) // vertex
   {
      const int off = nvd*(Info/64);
      dofs.SetSize(nvd);
      for (int i = 0; i < nvd; i++)
      {
         dofs[i] = off + i;
      }
   }
   else
   {
      int v[4], e[4], eo[4], f[1], fo[1];
      // New param
      int facet[5], faceto[5];
       
      int av = 0, nv = 0, ae = 0, ne = 0, nf = 0;
      // New param
      int af = 0, nfacet = 0;
       
      Geometry::Type fg[1];
      // New
      //Geometry::Type facetg[5];


      switch (Geom)
      {
         case Geometry::SEGMENT:
         {
            GetNVE<Geometry::SEGMENT>(av, ae);
            GetEdge<Geometry::SEGMENT>(nv, v, ne, e[0], eo[0], Info);
            break;
         }

         case Geometry::TRIANGLE:
         {
            GetNVE<Geometry::TRIANGLE>(av, ae);
            switch (SDim)
            {
               case 1:
                  GetEdge<Geometry::TRIANGLE>(nv, v, ne, e[0], eo[0], Info);
                  break;
               case 2:
                  GetFace<Geometry::TRIANGLE,Geometry::TRIANGLE>(
                     nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                  break;
               default:
                  goto not_supp;
            }
            break;
         }

         case Geometry::SQUARE:
         {
            GetNVE<Geometry::SQUARE>(av, ae);
            switch (SDim)
            {
               case 1:
                  GetEdge<Geometry::SQUARE>(nv, v, ne, e[0], eo[0], Info);
                  break;
               case 2:
                  GetFace<Geometry::SQUARE,Geometry::SQUARE>(
                     nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                  break;
               default:
                  goto not_supp;
            }
            break;
         }

         case Geometry::TETRAHEDRON:
         {
            GetNVE<Geometry::TETRAHEDRON>(av, ae);
            switch (SDim)
            {
               case 1:
                  GetEdge<Geometry::TETRAHEDRON>(nv, v, ne, e[0], eo[0], Info);
                  break;
               case 2:
                  GetFace<Geometry::TETRAHEDRON,Geometry::TRIANGLE>(
                     nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                  break;
               default:
                  goto not_supp;
            }
            break;
         }

         case Geometry::CUBE:
         {
            GetNVE<Geometry::CUBE>(av, ae);
            switch (SDim)
            {
               case 1:
                  GetEdge<Geometry::CUBE>(nv, v, ne, e[0], eo[0], Info);
                  break;
               case 2:
                  GetFace<Geometry::CUBE,Geometry::SQUARE>(
                     nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                  break;
               default:
                  goto not_supp;
            }
            break;
         }
              
          case Geometry::PENTATOPE:
          {
             GetNVE<Geometry::PENTATOPE>(av, ae);
             GetNF<Geometry::PENTATOPE>(af);
             switch (SDim)
             {
                case 1:
                   //GetEdge<Geometry::PENTATOPE>(nv, v, ne, e[0], eo[0], Info);
                   break;
                case 2:
                   //GetFace<Geometry::PENTATOPE,Geometry::TRIANGLE>(
                      //nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                   break;
                 case 3:
                    //GetFacet();//nv, v, ne, e, eo, nf, f, fg, fo, nfacet, facet[0], facetg[0], faceto[0], Info);
                     //nv, v, ne, e, eo, nf, f[0], fg[0], fo[0], Info);
                    break;
                default:
                   goto not_supp;
             }
             break;
          }

         default:
            MFEM_ABORT("invalid Geom = " << Geom);
      }

      int ned = (ne > 0) ? DofForGeometry(Geometry::SEGMENT) : 0;

      // add vertex dofs
      dofs.SetSize(nv*nvd+ne*ned);
      for (int i = 0; i < nv; i++)
      {
         for (int j = 0; j < nvd; j++)
         {
            dofs[i*nvd+j] = v[i]*nvd+j;
         }
      }
      int l_off = nv*nvd, g_off = av*nvd;

      // add edge dofs
      if (ned > 0)
      {
         for (int i = 0; i < ne; i++)
         {
            const int *ed = DofOrderForOrientation(Geometry::SEGMENT,
                                                   eo[i] ? -1 : 1);
            for (int j = 0; j < ned; j++)
            {
               dofs[l_off+i*ned+j] =
                  ed[j] >= 0 ?
                  g_off+e[i]*ned+ed[j] :
                  -1-(g_off+e[i]*ned+(-1-ed[j]));
            }
         }
         l_off += ne*ned;
         g_off += ae*ned;
      }

      // add face dofs
      if (nf > 0)
      {
         const int nfd = DofForGeometry(fg[0]); // assume same face geometry
         dofs.SetSize(dofs.Size()+nf*nfd);
         for (int i = 0; i < nf; i++)
         {
            const int *fd = DofOrderForOrientation(fg[i], fo[i]);
            for (int j = 0; j < nfd; j++)
            {
               dofs[l_off+i*nfd+j] =
                  fd[j] >= 0 ?
                  g_off+f[i]*nfd+fd[j] :
                  -1-(g_off+f[i]*nfd+(-1-fd[j]));
            }
         }
          l_off += nf*nfd;
          g_off += af*nfd;
      }

      // add volume dofs ...
//       if (nfacet > 0)
//       {
//          const int nfacetd = DofForGeometry(facetg[0]); // assume same face geometry
//          dofs.SetSize(dofs.Size()+nfacet*nfacetd);
//          for (int i = 0; i < nfacet; i++)
//          {
//             const int *facetd = DofOrderForOrientation(facetg[i], faceto[i]);
//             for (int j = 0; j < nfacetd; j++)
//             {
//                dofs[l_off+i*nfacetd+j] =
//                   facetd[j] >= 0 ?
//                   g_off+facet[i]*nfacetd+facetd[j] :
//                   -1-(g_off+facet[i]*nfacetd+(-1-facetd[j]));
//             }
//          }
//       }
       
   }
   return;

not_supp:
   MFEM_ABORT("Geom = " << Geometry::Name[Geom] <<
              ", SDim = " << SDim << " is not supported");
}

const FiniteElement *
LinearFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return &PointFE;
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      case Geometry::PRISM:       return &WedgeFE;
      case Geometry::PENTATOPE:   return &PentatopeFE;
      case Geometry::TESSERACT:   return &TesseractFE;
      default:
         mfem_error ("LinearFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int LinearFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 1;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      case Geometry::PRISM:       return 0;
      case Geometry::PENTATOPE:   return 0;
      case Geometry::TESSERACT:   return 0;
      default:
         mfem_error ("LinearFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *LinearFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   return NULL;
}


const FiniteElement *
QuadraticFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return &PointFE;
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      case Geometry::PRISM:       return &WedgeFE;
      case Geometry::PENTATOPE:   return &PentatopeFE;
      default:
         mfem_error ("QuadraticFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int QuadraticFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 1;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 1;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 1;
      case Geometry::PRISM:       return 0;
      case Geometry::PENTATOPE:   return 0;
      default:
         mfem_error ("QuadraticFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *QuadraticFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   static int indexes[] = { 0 };

   return indexes;
}


const FiniteElement *
QuadraticPosFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("QuadraticPosFECollection: unknown geometry type.");
   }
   return NULL; // Make some compilers happy
}

int QuadraticPosFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 1;
      case Geometry::SEGMENT:     return 1;
      case Geometry::SQUARE:      return 1;
      default:
         mfem_error ("QuadraticPosFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *QuadraticPosFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   static int indexes[] = { 0 };

   return indexes;
}


const FiniteElement *
CubicFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return &PointFE;
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      case Geometry::PRISM:       return &WedgeFE;
      default:
         mfem_error ("CubicFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int CubicFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 1;
      case Geometry::SEGMENT:     return 2;
      case Geometry::TRIANGLE:    return 1;
      case Geometry::SQUARE:      return 4;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 8;
      case Geometry::PRISM:       return 2;
      default:
         mfem_error ("CubicFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *CubicFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                     int Or) const
{
   if (GeomType == Geometry::SEGMENT)
   {
      static int ind_pos[] = { 0, 1 };
      static int ind_neg[] = { 1, 0 };

      if (Or < 0)
      {
         return ind_neg;
      }
      return ind_pos;
   }
   else if (GeomType == Geometry::TRIANGLE)
   {
      static int indexes[] = { 0 };

      return indexes;
   }
   else if (GeomType == Geometry::SQUARE)
   {
      static int sq_ind[8][4] = {{0, 1, 2, 3}, {0, 2, 1, 3},
         {2, 0, 3, 1}, {1, 0, 3, 2},
         {3, 2, 1, 0}, {3, 1, 2, 0},
         {1, 3, 0, 2}, {2, 3, 0, 1}
      };
      return sq_ind[Or];
   }

   return NULL;
}


const FiniteElement *
CrouzeixRaviartFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("CrouzeixRaviartFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int CrouzeixRaviartFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      default:
         mfem_error ("CrouzeixRaviartFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *CrouzeixRaviartFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   static int indexes[] = { 0 };

   return indexes;
}


const FiniteElement *
RT0_2DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("RT0_2DFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int RT0_2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      default:
         mfem_error ("RT0_2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int * RT0_2DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or) const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}


const FiniteElement *
RT1_2DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("RT1_2DFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int RT1_2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 2;
      case Geometry::TRIANGLE:    return 2;
      case Geometry::SQUARE:      return 4;
      default:
         mfem_error ("RT1_2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *RT1_2DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   static int ind_pos[] = {  0,  1 };
   static int ind_neg[] = { -2, -1 };

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}

const FiniteElement *
RT2_2DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("RT2_2DFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int RT2_2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 3;
      case Geometry::TRIANGLE:    return 6;
      case Geometry::SQUARE:      return 12;
      default:
         mfem_error ("RT2_2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *RT2_2DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   static int ind_pos[] = { 0, 1, 2 };
   static int ind_neg[] = { -3, -2, -1 };

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}


const FiniteElement *
Const2DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("Const2DFECollection: unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int Const2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 1;
      case Geometry::SQUARE:      return 1;
      default:
         mfem_error ("Const2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *Const2DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or) const
{
   return NULL;
}


const FiniteElement *
LinearDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("LinearDiscont2DFECollection: unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int LinearDiscont2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 3;
      case Geometry::SQUARE:      return 4;
      default:
         mfem_error ("LinearDiscont2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int * LinearDiscont2DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
GaussLinearDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("GaussLinearDiscont2DFECollection:"
                     " unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int GaussLinearDiscont2DFECollection::DofForGeometry(Geometry::Type GeomType)
const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 3;
      case Geometry::SQUARE:      return 4;
      default:
         mfem_error ("GaussLinearDiscont2DFECollection:"
                     " unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *GaussLinearDiscont2DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
P1OnQuadFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   if (GeomType != Geometry::SQUARE)
   {
      mfem_error ("P1OnQuadFECollection: unknown geometry type.");
   }
   return &QuadrilateralFE;
}

int P1OnQuadFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::SQUARE:      return 3;
      default:
         mfem_error ("P1OnQuadFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *P1OnQuadFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
QuadraticDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("QuadraticDiscont2DFECollection: unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int QuadraticDiscont2DFECollection::DofForGeometry(Geometry::Type GeomType)
const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 6;
      case Geometry::SQUARE:      return 9;
      default:
         mfem_error ("QuadraticDiscont2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *QuadraticDiscont2DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
QuadraticPosDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SQUARE:  return &QuadrilateralFE;
      default:
         mfem_error ("QuadraticPosDiscont2DFECollection: unknown geometry type.");
   }
   return NULL; // Make some compilers happy
}

int QuadraticPosDiscont2DFECollection::DofForGeometry(Geometry::Type GeomType)
const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::SQUARE:      return 9;
      default:
         mfem_error ("QuadraticPosDiscont2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}


const FiniteElement *
GaussQuadraticDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType)
const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("GaussQuadraticDiscont2DFECollection:"
                     " unknown geometry type.");
   }
   return &QuadrilateralFE; // Make some compilers happy
}

int GaussQuadraticDiscont2DFECollection::DofForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 6;
      case Geometry::SQUARE:      return 9;
      default:
         mfem_error ("GaussQuadraticDiscont2DFECollection:"
                     " unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *GaussQuadraticDiscont2DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
CubicDiscont2DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      default:
         mfem_error ("CubicDiscont2DFECollection: unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int CubicDiscont2DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 10;
      case Geometry::SQUARE:      return 16;
      default:
         mfem_error ("CubicDiscont2DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *CubicDiscont2DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
LinearNonConf3DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      default:
         mfem_error ("LinearNonConf3DFECollection: unknown geometry type.");
   }
   return &TriangleFE; // Make some compilers happy
}

int LinearNonConf3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 1;
      case Geometry::SQUARE:      return 1;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      default:
         mfem_error ("LinearNonConf3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *LinearNonConf3DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   static int indexes[] = { 0 };

   return indexes;
}


const FiniteElement *
Const3DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      case Geometry::PRISM:       return &WedgeFE;
      default:
         mfem_error ("Const3DFECollection: unknown geometry type.");
   }
   return &TetrahedronFE; // Make some compilers happy
}

int Const3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 1;
      case Geometry::CUBE:        return 1;
      case Geometry::PRISM:       return 1;
      default:
         mfem_error ("Const3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *Const3DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or) const
{
   return NULL;
}


const FiniteElement *
LinearDiscont3DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      default:
         mfem_error ("LinearDiscont3DFECollection: unknown geometry type.");
   }
   return &TetrahedronFE; // Make some compilers happy
}

int LinearDiscont3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 4;
      case Geometry::CUBE:        return 8;
      default:
         mfem_error ("LinearDiscont3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *LinearDiscont3DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}


const FiniteElement *
QuadraticDiscont3DFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      default:
         mfem_error ("QuadraticDiscont3DFECollection: unknown geometry type.");
   }
   return &TetrahedronFE; // Make some compilers happy
}

int QuadraticDiscont3DFECollection::DofForGeometry(Geometry::Type GeomType)
const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 10;
      case Geometry::CUBE:        return 27;
      default:
         mfem_error ("QuadraticDiscont3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *QuadraticDiscont3DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   return NULL;
}

const FiniteElement *
RefinedLinearFECollection::FiniteElementForGeometry(
   Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return &PointFE;
      case Geometry::SEGMENT:     return &SegmentFE;
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::CUBE:        return &ParallelepipedFE;
      default:
         mfem_error ("RefinedLinearFECollection: unknown geometry type.");
   }
   return &SegmentFE; // Make some compilers happy
}

int RefinedLinearFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 1;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 1;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 1;
      default:
         mfem_error ("RefinedLinearFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *RefinedLinearFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or) const
{
   static int indexes[] = { 0 };

   return indexes;
}


const FiniteElement *
ND1_3DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::CUBE:        return &HexahedronFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      default:
         mfem_error ("ND1_3DFECollection: unknown geometry type.");
   }
   return &HexahedronFE; // Make some compilers happy
}

int ND1_3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      default:
         mfem_error ("ND1_3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *ND1_3DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}

const FiniteElement *
ND1_4DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::PENTATOPE:   return &NedPentatopFE;
      default:
         mfem_error ("ND1_4DFECollection: unknown geometry type.");
   }
   return &NedPentatopFE; // Make some compilers happy
}

int ND1_4DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 1;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      case Geometry::PENTATOPE:   return 0;
      default:
         mfem_error ("ND1_4DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int * ND1_4DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or)
const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}


const FiniteElement *
ND2_4DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::PENTATOPE:   return &NedPentatopFE;
      default:
         mfem_error ("ND2_4DFECollection: unknown geometry type.");
   }
   return &NedPentatopFE; // Make some compilers happy
}

int ND2_4DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 2;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      case Geometry::PENTATOPE:   return 0;
      default:
         mfem_error ("ND2_4DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int * ND2_4DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or)
const
{
   static int ind_pos[] = { 0, 1 };
   static int ind_neg[] = { -2, -1};

   if (Or > 0)
   {
      return ind_pos;
   }
   return ind_neg;
}

const FiniteElement *
DivSkew1_4DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::PENTATOPE:   return &DivSkew0PentatopFE;
      default:
         mfem_error ("DivSkew1_4DFECollection: unknown geometry type 1.");
   }
   return &DivSkew0PentatopFE; // Make some compilers happy
}

int DivSkew1_4DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 1;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      case Geometry::PENTATOPE:   return 0;
      default:
         mfem_error ("DivSkew1_4DFECollection: unknown geometry type 2.");
   }
   return 0; // Make some compilers happy
}

const int * DivSkew1_4DFECollection::DofOrderForOrientation(
   Geometry::Type GeomType, int Or)
const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if (Or %2 == 0)
   {
      return ind_pos;
   }
   return ind_neg;
}

const FiniteElement *
RT0_3DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::CUBE:        return &HexahedronFE;
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      default:
         mfem_error ("RT0_3DFECollection: unknown geometry type.");
   }
   return &HexahedronFE; // Make some compilers happy
}

int RT0_3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 1;
      case Geometry::SQUARE:      return 1;
      case Geometry::TETRAHEDRON: return 0;
      case Geometry::CUBE:        return 0;
      default:
         mfem_error ("RT0_3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *RT0_3DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if ((GeomType == Geometry::TRIANGLE) || (GeomType == Geometry::SQUARE))
   {
      if (Or % 2 == 0)
      {
         return ind_pos;
      }
      return ind_neg;
   }
   return NULL;
}

const FiniteElement *
RT1_3DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TRIANGLE:    return &TriangleFE;
      case Geometry::SQUARE:      return &QuadrilateralFE;
      case Geometry::CUBE:        return &HexahedronFE;
      default:
         mfem_error ("RT1_3DFECollection: unknown geometry type.");
   }
   return &HexahedronFE; // Make some compilers happy
}

int RT1_3DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 2;
      case Geometry::SQUARE:      return 4;
      case Geometry::CUBE:        return 12;
      default:
         mfem_error ("RT1_3DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int *RT1_3DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                      int Or) const
{
   if (GeomType == Geometry::SQUARE)
   {
      static int sq_ind[8][4] =
      {
         {0, 1, 2, 3}, {-1, -3, -2, -4},
         {2, 0, 3, 1}, {-2, -1, -4, -3},
         {3, 2, 1, 0}, {-4, -2, -3, -1},
         {1, 3, 0, 2}, {-3, -4, -1, -2}
      };

      return sq_ind[Or];
   }
   else
   {
      return NULL;
   }
}

const FiniteElement *
RT0_4DFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::TETRAHEDRON: return &TetrahedronFE;
      case Geometry::PENTATOPE: return &PentatopeFE;
      default:
         mfem_error ("RT0_4DFECollection: unknown geometry type.");
   }
   return &PentatopeFE; // Make some compilers happy
}

int RT0_4DFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return 0;
      case Geometry::SEGMENT:     return 0;
      case Geometry::TRIANGLE:    return 0;
      case Geometry::SQUARE:      return 0;
      case Geometry::TETRAHEDRON: return 1;
      case Geometry::CUBE:        return 0;
      case Geometry::PENTATOPE:   return 0;
      default:
         mfem_error ("RT0_4DFECollection: unknown geometry type.");
   }
   return 0; // Make some compilers happy
}

const int * RT0_4DFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                       int Or)
const
{
   static int ind_pos[] = { 0 };
   static int ind_neg[] = { -1 };

   if (GeomType == Geometry::TETRAHEDRON)
   {
      if (Or % 2 == 0) { return ind_pos; }
      return ind_neg;
   }
   return NULL;
}

H1_FECollection::H1_FECollection(const int p, const int dim, const int btype)
   : FiniteElementCollection(p)
   , dim(dim)
{
   MFEM_VERIFY(p >= 1, "H1_FECollection requires order >= 1.");
   MFEM_VERIFY(dim >= 0 && dim <= 4, "H1_FECollection requires 0 <= dim <= 4.");

   const int pm1 = p - 1, pm2 = pm1 - 1, pm3 = pm2 - 1, pm4 = pm3 - 1, pm5 = pm4 -1;

   int pt_type = BasisType::GetQuadrature1D(btype);
   b_type = BasisType::Check(btype);
   switch (btype)
   {
      case BasisType::GaussLobatto:
      {
         snprintf(h1_name, 32, "H1_%dD_P%d", dim, p);
         break;
      }
      case BasisType::Positive:
      {
         snprintf(h1_name, 32, "H1Pos_%dD_P%d", dim, p);
         break;
      }
      case BasisType::Serendipity:
      {
         snprintf(h1_name, 32, "H1Ser_%dD_P%d", dim, p);
         break;
      }
      default:
      {
         MFEM_VERIFY(Quadrature1D::CheckClosed(pt_type) !=
                     Quadrature1D::Invalid,
                     "unsupported BasisType: " << BasisType::Name(btype));

         snprintf(h1_name, 32, "H1@%c_%dD_P%d",
                  (int)BasisType::GetChar(btype), dim, p);
      }
   }

   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      H1_dof[g] = 0;
      H1_Elements[g] = NULL;
   }
   for (int i = 0; i < 2; i++)
   {
      SegDofOrd[i] = NULL;
   }
   for (int i = 0; i < 6; i++)
   {
      TriDofOrd[i] = NULL;
   }
   for (int i = 0; i < 8; i++)
   {
      QuadDofOrd[i] = NULL;
   }
   for (int i = 0; i < 24; i++)
   {
      TetDofOrd[i] = NULL;
   }

   H1_dof[Geometry::POINT] = 1;
   H1_Elements[Geometry::POINT] = new PointFiniteElement;

   if (dim >= 1)
   {
      H1_dof[Geometry::SEGMENT] = pm1;
      if (b_type == BasisType::Positive)
      {
         H1_Elements[Geometry::SEGMENT] = new H1Pos_SegmentElement(p);
      }
      else
      {
         H1_Elements[Geometry::SEGMENT] = new H1_SegmentElement(p, btype);
      }

      SegDofOrd[0] = new int[2*pm1];
      SegDofOrd[1] = SegDofOrd[0] + pm1;
      for (int i = 0; i < pm1; i++)
      {
         SegDofOrd[0][i] = i;
         SegDofOrd[1][i] = pm2 - i;
      }
   }

   if (dim >= 2)
   {
      H1_dof[Geometry::TRIANGLE] = (pm1*pm2)/2;
      H1_dof[Geometry::SQUARE] = pm1*pm1;
      if (b_type == BasisType::Positive)
      {
         H1_Elements[Geometry::TRIANGLE] = new H1Pos_TriangleElement(p);
         H1_Elements[Geometry::SQUARE] = new H1Pos_QuadrilateralElement(p);
      }
      else if (b_type == BasisType::Serendipity)
      {
         // Note: in fe_coll.hpp the DofForGeometry(Geometry::Type) method
         // returns H1_dof[GeomType], so we need to fix the value of H1_dof here
         // for the serendipity case.

         // formula for number of interior serendipity DoFs (when p>1)
         H1_dof[Geometry::SQUARE] = (pm3*pm2)/2;
         H1_Elements[Geometry::SQUARE] = new H1Ser_QuadrilateralElement(p);
         // allows for mixed tri/quad meshes
         H1_Elements[Geometry::TRIANGLE] = new H1Pos_TriangleElement(p);
      }
      else
      {
         H1_Elements[Geometry::TRIANGLE] = new H1_TriangleElement(p, btype);
         H1_Elements[Geometry::SQUARE] = new H1_QuadrilateralElement(p, btype);
      }

      const int &TriDof = H1_dof[Geometry::TRIANGLE];
      const int &QuadDof = H1_dof[Geometry::SQUARE];
      TriDofOrd[0] = new int[6*TriDof];
      for (int i = 1; i < 6; i++)
      {
         TriDofOrd[i] = TriDofOrd[i-1] + TriDof;
      }
      // see Mesh::GetTriOrientation in mesh/mesh.cpp
      for (int j = 0; j < pm2; j++)
      {
         for (int i = 0; i + j < pm2; i++)
         {
            int o = TriDof - ((pm1 - j)*(pm2 - j))/2 + i;
            //std::cout << "Orientation of Triangles" << std::endl;
            //std::cout << "i= " << i << " j= " << j << std::endl;
            //std::cout << "TriDof = " << TriDof << std::endl;
            //std::cout << "value o = " << o << std::endl;
            int k = pm3 - j - i;
            TriDofOrd[0][o] = o;  // (0,1,2)
            TriDofOrd[1][o] = TriDof - ((pm1-j)*(pm2-j))/2 + k;  // (1,0,2)
            TriDofOrd[2][o] = TriDof - ((pm1-i)*(pm2-i))/2 + k;  // (2,0,1)
            TriDofOrd[3][o] = TriDof - ((pm1-k)*(pm2-k))/2 + i;  // (2,1,0)
            TriDofOrd[4][o] = TriDof - ((pm1-k)*(pm2-k))/2 + j;  // (1,2,0)
            TriDofOrd[5][o] = TriDof - ((pm1-i)*(pm2-i))/2 + j;  // (0,2,1)
            //std::cout << "TriDofOrd = " << TriDofOrd[0][o] << ", " << TriDofOrd[1][o] << ", " << TriDofOrd[2][o] << ", " << TriDofOrd[3][o] << ", " << TriDofOrd[4][o] << ", " << TriDofOrd[5][o] << ", " << std::endl;
         }
      }

      QuadDofOrd[0] = new int[8*QuadDof];
      for (int i = 1; i < 8; i++)
      {
         QuadDofOrd[i] = QuadDofOrd[i-1] + QuadDof;
      }

      // For serendipity order >=4, the QuadDofOrd array must be re-defined. We
      // do this by computing the corresponding tensor product QuadDofOrd array
      // or two orders less, which contains enough DoFs for their serendipity
      // basis. This could be optimized.
      if (b_type == BasisType::Serendipity)
      {
         if (p < 4)
         {
            // no face dofs --> don't need to adjust QuadDofOrd
         }
         else  // p >= 4 --> have face dofs
         {
            // Exactly the same as tensor product case, but with all orders
            // reduced by 2 e.g. in case p=5 it builds a 2x2 array, even though
            // there are only 3 serendipity dofs.
            // In the tensor product case, the i and j index tensor directions,
            // and o index from 0 to (pm1)^2,

            for (int j = 0; j < pm3; j++)   // pm3 instead of pm1, etc
            {
               for (int i = 0; i < pm3; i++)
               {
                  int o = i + j*pm3;
                  QuadDofOrd[0][o] = i + j*pm3;  // (0,1,2,3)
                  QuadDofOrd[1][o] = j + i*pm3;  // (0,3,2,1)
                  QuadDofOrd[2][o] = j + (pm4 - i)*pm3;  // (1,2,3,0)
                  QuadDofOrd[3][o] = (pm4 - i) + j*pm3;  // (1,0,3,2)
                  QuadDofOrd[4][o] = (pm4 - i) + (pm4 - j)*pm3;  // (2,3,0,1)
                  QuadDofOrd[5][o] = (pm4 - j) + (pm4 - i)*pm3;  // (2,1,0,3)
                  QuadDofOrd[6][o] = (pm4 - j) + i*pm3;  // (3,0,1,2)
                  QuadDofOrd[7][o] = i + (pm4 - j)*pm3;  // (3,2,1,0)
               }
            }

         }
      }
      else // not serendipity
      {
         for (int j = 0; j < pm1; j++)
         {
            for (int i = 0; i < pm1; i++)
            {
               int o = i + j*pm1;
               QuadDofOrd[0][o] = i + j*pm1;  // (0,1,2,3)
               QuadDofOrd[1][o] = j + i*pm1;  // (0,3,2,1)
               QuadDofOrd[2][o] = j + (pm2 - i)*pm1;  // (1,2,3,0)
               QuadDofOrd[3][o] = (pm2 - i) + j*pm1;  // (1,0,3,2)
               QuadDofOrd[4][o] = (pm2 - i) + (pm2 - j)*pm1;  // (2,3,0,1)
               QuadDofOrd[5][o] = (pm2 - j) + (pm2 - i)*pm1;  // (2,1,0,3)
               QuadDofOrd[6][o] = (pm2 - j) + i*pm1;  // (3,0,1,2)
               QuadDofOrd[7][o] = i + (pm2 - j)*pm1;  // (3,2,1,0)
            }
         }
      }

      if (dim >= 3)
      {
         H1_dof[Geometry::TETRAHEDRON] = (TriDof*pm3)/3;
         H1_dof[Geometry::CUBE] = QuadDof*pm1;
         H1_dof[Geometry::PRISM] = TriDof*pm1;
         if (b_type == BasisType::Positive)
         {
            H1_Elements[Geometry::TETRAHEDRON] = new H1Pos_TetrahedronElement(p);
            H1_Elements[Geometry::CUBE] = new H1Pos_HexahedronElement(p);
            H1_Elements[Geometry::PRISM] = new H1Pos_WedgeElement(p);
         }
         else
         {
            H1_Elements[Geometry::TETRAHEDRON] = new H1_TetrahedronElement(p, btype);
            //H1_Elements[Geometry::TETRAHEDRON] = new H1_TetrahedronElement_Fuentes(p, btype);

            H1_Elements[Geometry::CUBE] = new H1_HexahedronElement(p, btype);
            H1_Elements[Geometry::PRISM] = new H1_WedgeElement(p, btype);
         }

         const int &TetDof = H1_dof[Geometry::TETRAHEDRON];
         TetDofOrd[0] = new int[24*TetDof];
         for (int i = 1; i < 24; i++)
         {
            TetDofOrd[i] = TetDofOrd[i-1] + TetDof;
         }
         // see Mesh::GetTetOrientation in mesh/mesh.cpp
         for (int k = 0; k < pm3; k++)
         {
            for (int j = 0; j + k < pm3; j++)
            {
               for (int i = 0; i + j + k < pm3; i++)
               {
//                  std::cout << "Orientation of Tets" << std::endl;
//                  std::cout << "i= " << i << " j= " << j << " k= " << k << std::endl;
//                  std::cout << "TetDof = " << TetDof << std::endl;

                  int l = pm4 - k - j - i;
                  int o   = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (j * (2 * p - 5 - j - 2 * k)) / 2 + i;
                  //std::cout << "o= " << o << std::endl;
                  int o1  = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (k * (2 * p - 5 - k - 2 * j)) / 2 + i;
                  int o2  = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (k * (2 * p - 5 - k - 2 * i)) / 2 + j;
                  int o3  = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (i * (2 * p - 5 - i - 2 * k)) / 2 + j;
                  int o4  = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (i * (2 * p - 5 - i - 2 * j)) / 2 + k;
                  int o5  = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (j * (2 * p - 5 - j - 2 * i)) / 2 + k;
                  int o6  = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (l * (2 * p - 5 - l - 2 * k)) / 2 + j;
                  int o7  = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (k * (2 * p - 5 - k - 2 * l)) / 2 + j;
                  int o8  = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (j * (2 * p - 5 - j - 2 * l)) / 2 + k;
                  int o9  = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (l * (2 * p - 5 - l - 2 * j)) / 2 + k;
                  int o10 = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (k * (2 * p - 5 - k - 2 * j)) / 2 + l;
                  int o11 = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (j * (2 * p - 5 - j - 2 * k)) / 2 + l;
                  int o12 = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (l * (2 * p - 5 - l - 2 * i)) / 2 + k;
                  int o13 = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (i * (2 * p - 5 - i - 2 * l)) / 2 + k;
                  int o14 = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (i * (2 * p - 5 - i - 2 * k)) / 2 + l;
                  int o15 = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (k * (2 * p - 5 - k - 2 * i)) / 2 + l;
                  int o16 = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (k * (2 * p - 5 - k - 2 * l)) / 2 + i;
                  int o17 = TetDof - ((pm1 - k) * (pm2 - k) * (pm3 - k)) / 6
                            + (l * (2 * p - 5 - l - 2 * k)) / 2 + i;
                  int o18 = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (j * (2 * p - 5 - j - 2 * i)) / 2 + l;
                  int o19 = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (i * (2 * p - 5 - i - 2 * j)) / 2 + l;
                  int o20 = TetDof - ((pm1 - j) * (pm2 - j) * (pm3 - j)) / 6
                            + (l * (2 * p - 5 - l - 2 * j)) / 2 + i;
                  int o21 = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (j * (2 * p - 5 - j - 2 * l)) / 2 + i;
                  int o22 = TetDof - ((pm1 - l) * (pm2 - l) * (pm3 - l)) / 6
                            + (i * (2 * p - 5 - i - 2 * l)) / 2 + j;
                  int o23 = TetDof - ((pm1 - i) * (pm2 - i) * (pm3 - i)) / 6
                            + (l * (2 * p - 5 - l - 2 * i)) / 2 + j;
                  TetDofOrd[ 0][o] = o;   // (0,1,2,3)
                  TetDofOrd[ 1][o] = o1;  // (0,1,3,2)
                  TetDofOrd[ 2][o] = o2;  // (0,2,3,1)
                  TetDofOrd[ 3][o] = o3;  // (0,2,1,3)
                  TetDofOrd[ 4][o] = o4;  // (0,3,1,2)
                  TetDofOrd[ 5][o] = o5;  // (0,3,2,1)
                  TetDofOrd[ 6][o] = o6;  // (1,2,0,3)
                  TetDofOrd[ 7][o] = o7;  // (1,2,3,0)
                  TetDofOrd[ 8][o] = o8;  // (1,3,2,0)
                  TetDofOrd[ 9][o] = o9;  // (1,3,0,2)
                  TetDofOrd[10][o] = o10; // (1,0,3,2)
                  TetDofOrd[11][o] = o11; // (1,0,2,3)
                  TetDofOrd[12][o] = o12; // (2,3,0,1)
                  TetDofOrd[13][o] = o13; // (2,3,1,0)
                  TetDofOrd[14][o] = o14; // (2,0,1,3)
                  TetDofOrd[15][o] = o15; // (2,0,3,1)
                  TetDofOrd[16][o] = o16; // (2,1,3,0)
                  TetDofOrd[17][o] = o17; // (2,1,0,3)
                  TetDofOrd[18][o] = o18; // (3,0,2,1)
                  TetDofOrd[19][o] = o19; // (3,0,1,2)
                  TetDofOrd[20][o] = o20; // (3,1,0,2)
                  TetDofOrd[21][o] = o21; // (3,1,2,0)
                  TetDofOrd[22][o] = o22; // (3,2,1,0)
                  TetDofOrd[23][o] = o23; // (3,2,0,1)
//                   for (int i = 0 ; i<24; i++) {
//                       std::cout << "TetDofOrd: " << TetDofOrd[i][o] << std::endl; //"," << TetDofOrd[i][o] << "," << TetDofOrd[i][o] << "," << TetDofOrd[i][o] << std::endl;
//                   }
               }
            }
         }


         if (dim >= 4)
         {
            H1_dof[Geometry::PENTATOPE] = (TriDof*pm3*pm4)/12;
            H1_dof[Geometry::TESSERACT] = QuadDof*pm1*pm1;
            if (b_type == BasisType::Positive)
            {
               mfem_error("H1_FECollection: BasisType::Positive not implemented");
            }
            else
            {
               //H1_Elements[Geometry::PENTATOPE] = new H1_PentatopeElement(p, pt_type);
               H1_Elements[Geometry::PENTATOPE] = new H1_PentatopeElement_Barycentric(p, pt_type);

            }
             // new code
             const int &PentDof = H1_dof[Geometry::PENTATOPE];
             PentDofOrd[0] = new int[120*PentDof];
             for (int i = 1; i < 120; i++)
             {
                PentDofOrd[i] = PentDofOrd[i-1] + PentDof;
             }
             // see Mesh::GetPentOrientation in mesh/mesh.cpp
             for (int n = 0; n < pm4; n++)
             {
                for (int k = 0; k + n < pm4; k++)
                {
                   for (int j = 0; j + k + n < pm4; j++)
                   {
                       for (int i = 0; i + j + k + n < pm4; i++)
                       {
                           // insert code here
                           // Alias p
                           int P = p - 5;
                           //std::cout << "value of p= " << P << std::endl;
                           int l = pm5 - n - k - j - i;
                           
                           int o = i + k*(2*P + P^2/2 + 11/6) - j*k - j*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + j*(P + 3/2) - j^2/2 + k^3/6 - n^4/24 - k*n*(P + 2) + 1;
                           int o1 = i + k*(2*P + P^2/2 + 11/6) - j*k - j*l + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + j*(P + 3/2) - j^2/2 + k^3/6 - l^4/24 - k*l*(P + 2) + 1;
                           int o2 = i + n*(2*P + P^2/2 + 11/6) - j*k - j*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + j*(P + 3/2) - j^2/2 - k^4/24 + n^3/6 - k*n*(P + 2) + 1;
                           int o3 = i + n*(2*P + P^2/2 + 11/6) - j*l - j*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + j*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - j^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o4 = i + l*(2*P + P^2/2 + 11/6) - j*k - j*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + j*(P + 3/2) - j^2/2 - k^4/24 + l^3/6 - k*l*(P + 2) + 1;
                           int o5 = i + l*(2*P + P^2/2 + 11/6) - j*l - j*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + j*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - j^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o6 = i + j*(2*P + P^2/2 + 11/6) - j*k - k*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + k*(P + 3/2) + j^3/6 - k^2/2 - n^4/24 - j*n*(P + 2) + 1;
                           int o7 = i + j*(2*P + P^2/2 + 11/6) - j*k - k*l - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + k*(P + 3/2) + j^3/6 - k^2/2 - l^4/24 - j*l*(P + 2) + 1;
                           int o8 = i + n*(2*P + P^2/2 + 11/6) - j*k - k*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + k*(P + 3/2) - j^4/24 - k^2/2 + n^3/6 - j*n*(P + 2) + 1;
                           int o9 = i + n*(2*P + P^2/2 + 11/6) - k*l - k*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + (l*n^2)/2 + (l^2*n)/2 + k*(P + 3/2) - k^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o10 = i + l*(2*P + P^2/2 + 11/6) - j*k - k*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + k*(P + 3/2) - j^4/24 - k^2/2 + l^3/6 - j*l*(P + 2) + 1;
                           int o11 = i + l*(2*P + P^2/2 + 11/6) - k*l - k*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (l*n^2)/2 + (l^2*n)/2 + k*(P + 3/2) - k^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o12 = i + j*(2*P + P^2/2 + 11/6) - j*n - k*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + n*(P + 3/2) + j^3/6 - k^4/24 - n^2/2 - j*k*(P + 2) + 1;
                           int o13 = i + j*(2*P + P^2/2 + 11/6) - j*n - l*n - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + n*(P + 3/2) + j^3/6 - l^4/24 - n^2/2 - j*l*(P + 2) + 1;
                           int o14 = i + k*(2*P + P^2/2 + 11/6) - j*n - k*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + n*(P + 3/2) - j^4/24 + k^3/6 - n^2/2 - j*k*(P + 2) + 1;
                           int o15 = i + k*(2*P + P^2/2 + 11/6) - k*n - l*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + n*(P + 3/2) + k^3/6 - l^4/24 - n^2/2 - k*l*(P + 2) + 1;
                           int o16 = i + l*(2*P + P^2/2 + 11/6) - j*n - l*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + n*(P + 3/2) - j^4/24 + l^3/6 - n^2/2 - j*l*(P + 2) + 1;
                           int o17 = i + l*(2*P + P^2/2 + 11/6) - k*n - l*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + n*(P + 3/2) - k^4/24 + l^3/6 - n^2/2 - k*l*(P + 2) + 1;
                           int o18 = i + j*(2*P + P^2/2 + 11/6) - j*l - k*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + l*(P + 3/2) + j^3/6 - k^4/24 - l^2/2 - j*k*(P + 2) + 1;
                           int o19 = i + j*(2*P + P^2/2 + 11/6) - j*l - l*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + l*(P + 3/2) + j^3/6 - l^2/2 - n^4/24 - j*n*(P + 2) + 1;
                           int o20 = i + k*(2*P + P^2/2 + 11/6) - j*l - k*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + l*(P + 3/2) - j^4/24 + k^3/6 - l^2/2 - j*k*(P + 2) + 1;
                           int o21 = i + k*(2*P + P^2/2 + 11/6) - k*l - l*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + l*(P + 3/2) + k^3/6 - l^2/2 - n^4/24 - k*n*(P + 2) + 1;
                           int o22 = i + n*(2*P + P^2/2 + 11/6) - j*l - l*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + l*(P + 3/2) - j^4/24 - l^2/2 + n^3/6 - j*n*(P + 2) + 1;
                           int o23 = i + n*(2*P + P^2/2 + 11/6) - k*l - l*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + l*(P + 3/2) - k^4/24 - l^2/2 + n^3/6 - k*n*(P + 2) + 1;
                           int o24 = j + k*(2*P + P^2/2 + 11/6) - i*k - i*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (k*n^2)/2 + (k^2*n)/2 - i^2/2 + k^3/6 - n^4/24 - k*n*(P + 2) + 1;
                           int o25 = j + k*(2*P + P^2/2 + 11/6) - i*k - i*l + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + i*(P + 3/2) - i^2/2 + k^3/6 - l^4/24 - k*l*(P + 2) + 1;
                           int o26 = j + n*(2*P + P^2/2 + 11/6) - i*k - i*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (k*n^2)/2 + (k^2*n)/2 - i^2/2 - k^4/24 + n^3/6 - k*n*(P + 2) + 1;
                           int o27 = j + n*(2*P + P^2/2 + 11/6) - i*l - i*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - i^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o28 = j + l*(2*P + P^2/2 + 11/6) - i*k - i*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + i*(P + 3/2) - i^2/2 - k^4/24 + l^3/6 - k*l*(P + 2) + 1;
                           int o29 = j + l*(2*P + P^2/2 + 11/6) - i*l - i*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - i^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o30 = j + i*(2*P + P^2/2 + 11/6) - i*k - k*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + k*(P + 3/2) + i^3/6 - k^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o31 = j + i*(2*P + P^2/2 + 11/6) - i*k - k*l - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + k*(P + 3/2) + i^3/6 - k^2/2 - l^4/24 - i*l*(P + 2) + 1;
                           int o32 = j + n*(2*P + P^2/2 + 11/6) - i*k - k*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + k*(P + 3/2) - i^4/24 - k^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o33 = j + n*(2*P + P^2/2 + 11/6) - k*l - k*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + (l*n^2)/2 + (l^2*n)/2 + k*(P + 3/2) - k^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o34 = j + l*(2*P + P^2/2 + 11/6) - i*k - k*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + k*(P + 3/2) - i^4/24 - k^2/2 + l^3/6 - i*l*(P + 2) + 1;
                           int o35 = j + l*(2*P + P^2/2 + 11/6) - k*l - k*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (l*n^2)/2 + (l^2*n)/2 + k*(P + 3/2) - k^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o36 = j + i*(2*P + P^2/2 + 11/6) - i*n - k*n - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + n*(P + 3/2) + i^3/6 - k^4/24 - n^2/2 - i*k*(P + 2) + 1;
                           int o37 = j + i*(2*P + P^2/2 + 11/6) - i*n - l*n - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + n*(P + 3/2) + i^3/6 - l^4/24 - n^2/2 - i*l*(P + 2) + 1;
                           int o38 = j + k*(2*P + P^2/2 + 11/6) - i*n - k*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + n*(P + 3/2) - i^4/24 + k^3/6 - n^2/2 - i*k*(P + 2) + 1;
                           int o39 = j + k*(2*P + P^2/2 + 11/6) - k*n - l*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + n*(P + 3/2) + k^3/6 - l^4/24 - n^2/2 - k*l*(P + 2) + 1;
                           int o40 = j + l*(2*P + P^2/2 + 11/6) - i*n - l*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + n*(P + 3/2) - i^4/24 + l^3/6 - n^2/2 - i*l*(P + 2) + 1;
                           int o41 = j + l*(2*P + P^2/2 + 11/6) - k*n - l*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + n*(P + 3/2) - k^4/24 + l^3/6 - n^2/2 - k*l*(P + 2) + 1;
                           int o42 = j + i*(2*P + P^2/2 + 11/6) - i*l - k*l - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + l*(P + 3/2) + i^3/6 - k^4/24 - l^2/2 - i*k*(P + 2) + 1;
                           int o43 = j + i*(2*P + P^2/2 + 11/6) - i*l - l*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + l*(P + 3/2) + i^3/6 - l^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o44 = j + k*(2*P + P^2/2 + 11/6) - i*l - k*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + l*(P + 3/2) - i^4/24 + k^3/6 - l^2/2 - i*k*(P + 2) + 1;
                           int o45 = j + k*(2*P + P^2/2 + 11/6) - k*l - l*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + l*(P + 3/2) + k^3/6 - l^2/2 - n^4/24 - k*n*(P + 2) + 1;
                           int o46 = j + n*(2*P + P^2/2 + 11/6) - i*l - l*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + l*(P + 3/2) - i^4/24 - l^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o47 = j + n*(2*P + P^2/2 + 11/6) - k*l - l*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + l*(P + 3/2) - k^4/24 - l^2/2 + n^3/6 - k*n*(P + 2) + 1;
                           int o48 = k + j*(2*P + P^2/2 + 11/6) - i*j - i*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - n^4/24 - j*n*(P + 2) + 1;
                           int o49 = k + j*(2*P + P^2/2 + 11/6) - i*j - i*l - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - l^4/24 - j*l*(P + 2) + 1;
                           int o50 = k + n*(2*P + P^2/2 + 11/6) - i*j - i*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + n^3/6 - j*n*(P + 2) + 1;
                           int o51 = k + n*(2*P + P^2/2 + 11/6) - i*l - i*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - i^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o52 = k + l*(2*P + P^2/2 + 11/6) - i*j - i*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + l^3/6 - j*l*(P + 2) + 1;
                           int o53 = k + l*(2*P + P^2/2 + 11/6) - i*l - i*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - i^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o54 = k + i*(2*P + P^2/2 + 11/6) - i*j - j*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o55 = k + i*(2*P + P^2/2 + 11/6) - i*j - j*l - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - l^4/24 - i*l*(P + 2) + 1;
                           int o56 = k + n*(2*P + P^2/2 + 11/6) - i*j - j*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o57 = k + n*(2*P + P^2/2 + 11/6) - j*l - j*n + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - n^2*(P/2 + 1) - l^2*((5*P)/4 + P^2/4 + 35/24) + j*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - j^2/2 - l^4/24 + n^3/6 - l*n*(P + 2) + 1;
                           int o58 = k + l*(2*P + P^2/2 + 11/6) - i*j - j*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + l^3/6 - i*l*(P + 2) + 1;
                           int o59 = k + l*(2*P + P^2/2 + 11/6) - j*l - j*n - l^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + j*(P + 3/2) + (l*n^2)/2 + (l^2*n)/2 - j^2/2 + l^3/6 - n^4/24 - l*n*(P + 2) + 1;
                           int o60 = k + i*(2*P + P^2/2 + 11/6) - i*n - j*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + n*(P + 3/2) + i^3/6 - j^4/24 - n^2/2 - i*j*(P + 2) + 1;
                           int o61 = k + i*(2*P + P^2/2 + 11/6) - i*n - l*n - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + n*(P + 3/2) + i^3/6 - l^4/24 - n^2/2 - i*l*(P + 2) + 1;
                           int o62 = k + j*(2*P + P^2/2 + 11/6) - i*n - j*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + n*(P + 3/2) - i^4/24 + j^3/6 - n^2/2 - i*j*(P + 2) + 1;
                           int o63 = k + j*(2*P + P^2/2 + 11/6) - j*n - l*n - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + n*(P + 3/2) + j^3/6 - l^4/24 - n^2/2 - j*l*(P + 2) + 1;
                           int o64 = k + l*(2*P + P^2/2 + 11/6) - i*n - l*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + n*(P + 3/2) - i^4/24 + l^3/6 - n^2/2 - i*l*(P + 2) + 1;
                           int o65 = k + l*(2*P + P^2/2 + 11/6) - j*n - l*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + n*(P + 3/2) - j^4/24 + l^3/6 - n^2/2 - j*l*(P + 2) + 1;
                           int o66 = k + i*(2*P + P^2/2 + 11/6) - i*l - j*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + l*(P + 3/2) + i^3/6 - j^4/24 - l^2/2 - i*j*(P + 2) + 1;
                           int o67 = k + i*(2*P + P^2/2 + 11/6) - i*l - l*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + l*(P + 3/2) + i^3/6 - l^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o68 = k + j*(2*P + P^2/2 + 11/6) - i*l - j*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + l*(P + 3/2) - i^4/24 + j^3/6 - l^2/2 - i*j*(P + 2) + 1;
                           int o69 = k + j*(2*P + P^2/2 + 11/6) - j*l - l*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + l*(P + 3/2) + j^3/6 - l^2/2 - n^4/24 - j*n*(P + 2) + 1;
                           int o70 = k + n*(2*P + P^2/2 + 11/6) - i*l - l*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + l*(P + 3/2) - i^4/24 - l^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o71 = k + n*(2*P + P^2/2 + 11/6) - j*l - l*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + l*(P + 3/2) - j^4/24 - l^2/2 + n^3/6 - j*n*(P + 2) + 1;
                           int o72 = n + j*(2*P + P^2/2 + 11/6) - i*j - i*k + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - k^4/24 - j*k*(P + 2) + 1;
                           int o73 = n + j*(2*P + P^2/2 + 11/6) - i*j - i*l - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - l^4/24 - j*l*(P + 2) + 1;
                           int o74 = n + k*(2*P + P^2/2 + 11/6) - i*j - i*k + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + k^3/6 - j*k*(P + 2) + 1;
                           int o75 = n + k*(2*P + P^2/2 + 11/6) - i*k - i*l + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + i*(P + 3/2) - i^2/2 + k^3/6 - l^4/24 - k*l*(P + 2) + 1;
                           int o76 = n + l*(2*P + P^2/2 + 11/6) - i*j - i*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + l^3/6 - j*l*(P + 2) + 1;
                           int o77 = n + l*(2*P + P^2/2 + 11/6) - i*k - i*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + i*(P + 3/2) - i^2/2 - k^4/24 + l^3/6 - k*l*(P + 2) + 1;
                           int o78 = n + i*(2*P + P^2/2 + 11/6) - i*j - j*k - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - k^4/24 - i*k*(P + 2) + 1;
                           int o79 = n + i*(2*P + P^2/2 + 11/6) - i*j - j*l - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - l^4/24 - i*l*(P + 2) + 1;
                           int o80 = n + k*(2*P + P^2/2 + 11/6) - i*j - j*k + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + k^3/6 - i*k*(P + 2) + 1;
                           int o81 = n + k*(2*P + P^2/2 + 11/6) - j*k - j*l + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - k^2*(P/2 + 1) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + j*(P + 3/2) - j^2/2 + k^3/6 - l^4/24 - k*l*(P + 2) + 1;
                           int o82 = n + l*(2*P + P^2/2 + 11/6) - i*j - j*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + l^3/6 - i*l*(P + 2) + 1;
                           int o83 = n + l*(2*P + P^2/2 + 11/6) - j*k - j*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - l^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*l^2)/2 + (k^2*l)/2 + j*(P + 3/2) - j^2/2 - k^4/24 + l^3/6 - k*l*(P + 2) + 1;
                           int o84 = n + i*(2*P + P^2/2 + 11/6) - i*k - j*k + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + k*(P + 3/2) + i^3/6 - j^4/24 - k^2/2 - i*j*(P + 2) + 1;
                           int o85 = n + i*(2*P + P^2/2 + 11/6) - i*k - k*l - i^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + k*(P + 3/2) + i^3/6 - k^2/2 - l^4/24 - i*l*(P + 2) + 1;
                           int o86 = n + j*(2*P + P^2/2 + 11/6) - i*k - j*k + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + k*(P + 3/2) - i^4/24 + j^3/6 - k^2/2 - i*j*(P + 2) + 1;
                           int o87 = n + j*(2*P + P^2/2 + 11/6) - j*k - k*l - j^2*(P/2 + 1) + l*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + l^3*(P/6 + 5/12) - l^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + k*(P + 3/2) + j^3/6 - k^2/2 - l^4/24 - j*l*(P + 2) + 1;
                           int o88 = n + l*(2*P + P^2/2 + 11/6) - i*k - k*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - l^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*l^2)/2 + (i^2*l)/2 + k*(P + 3/2) - i^4/24 - k^2/2 + l^3/6 - i*l*(P + 2) + 1;
                           int o89 = n + l*(2*P + P^2/2 + 11/6) - j*k - k*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - l^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*l^2)/2 + (j^2*l)/2 + k*(P + 3/2) - j^4/24 - k^2/2 + l^3/6 - j*l*(P + 2) + 1;
                           int o90 = n + i*(2*P + P^2/2 + 11/6) - i*l - j*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + l*(P + 3/2) + i^3/6 - j^4/24 - l^2/2 - i*j*(P + 2) + 1;
                           int o91 = n + i*(2*P + P^2/2 + 11/6) - i*l - k*l - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + l*(P + 3/2) + i^3/6 - k^4/24 - l^2/2 - i*k*(P + 2) + 1;
                           int o92 = n + j*(2*P + P^2/2 + 11/6) - i*l - j*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + l*(P + 3/2) - i^4/24 + j^3/6 - l^2/2 - i*j*(P + 2) + 1;
                           int o93 = n + j*(2*P + P^2/2 + 11/6) - j*l - k*l + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + l*(P + 3/2) + j^3/6 - k^4/24 - l^2/2 - j*k*(P + 2) + 1;
                           int o94 = n + k*(2*P + P^2/2 + 11/6) - i*l - k*l + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + l*(P + 3/2) - i^4/24 + k^3/6 - l^2/2 - i*k*(P + 2) + 1;
                           int o95 = n + k*(2*P + P^2/2 + 11/6) - j*l - k*l + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + l*(P + 3/2) - j^4/24 + k^3/6 - l^2/2 - j*k*(P + 2) + 1;
                           int o96 = l + j*(2*P + P^2/2 + 11/6) - i*j - i*k + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - k^4/24 - j*k*(P + 2) + 1;
                           int o97 = l + j*(2*P + P^2/2 + 11/6) - i*j - i*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + i*(P + 3/2) - i^2/2 + j^3/6 - n^4/24 - j*n*(P + 2) + 1;
                           int o98 = l + k*(2*P + P^2/2 + 11/6) - i*j - i*k + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + k^3/6 - j*k*(P + 2) + 1;
                           int o99 = l + k*(2*P + P^2/2 + 11/6) - i*k - i*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (k*n^2)/2 + (k^2*n)/2 - i^2/2 + k^3/6 - n^4/24 - k*n*(P + 2) + 1;
                           int o100 = l + n*(2*P + P^2/2 + 11/6) - i*j - i*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + i*(P + 3/2) - i^2/2 - j^4/24 + n^3/6 - j*n*(P + 2) + 1;
                           int o101 = l + n*(2*P + P^2/2 + 11/6) - i*k - i*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + i*(P + 3/2) + (k*n^2)/2 + (k^2*n)/2 - i^2/2 - k^4/24 + n^3/6 - k*n*(P + 2) + 1;
                           int o102 = l + i*(2*P + P^2/2 + 11/6) - i*j - j*k - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - k^4/24 - i*k*(P + 2) + 1;
                           int o103 = l + i*(2*P + P^2/2 + 11/6) - i*j - j*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + j*(P + 3/2) + i^3/6 - j^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o104 = l + k*(2*P + P^2/2 + 11/6) - i*j - j*k + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + k^3/6 - i*k*(P + 2) + 1;
                           int o105 = l + k*(2*P + P^2/2 + 11/6) - j*k - j*n - k^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + j*(P + 3/2) - j^2/2 + k^3/6 - n^4/24 - k*n*(P + 2) + 1;
                           int o106 = l + n*(2*P + P^2/2 + 11/6) - i*j - j*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + j*(P + 3/2) - i^4/24 - j^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o107 = l + n*(2*P + P^2/2 + 11/6) - j*k - j*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - n^2*(P/2 + 1) - k^2*((5*P)/4 + P^2/4 + 35/24) + (k*n^2)/2 + (k^2*n)/2 + j*(P + 3/2) - j^2/2 - k^4/24 + n^3/6 - k*n*(P + 2) + 1;
                           int o108 = l + i*(2*P + P^2/2 + 11/6) - i*k - j*k + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + k*(P + 3/2) + i^3/6 - j^4/24 - k^2/2 - i*j*(P + 2) + 1;
                           int o109 = l + i*(2*P + P^2/2 + 11/6) - i*k - k*n - i^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + k*(P + 3/2) + i^3/6 - k^2/2 - n^4/24 - i*n*(P + 2) + 1;
                           int o110 = l + j*(2*P + P^2/2 + 11/6) - i*k - j*k + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + k*(P + 3/2) - i^4/24 + j^3/6 - k^2/2 - i*j*(P + 2) + 1;
                           int o111 = l + j*(2*P + P^2/2 + 11/6) - j*k - k*n - j^2*(P/2 + 1) + n*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + n^3*(P/6 + 5/12) - n^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + k*(P + 3/2) + j^3/6 - k^2/2 - n^4/24 - j*n*(P + 2) + 1;
                           int o112 = l + n*(2*P + P^2/2 + 11/6) - i*k - k*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - n^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*n^2)/2 + (i^2*n)/2 + k*(P + 3/2) - i^4/24 - k^2/2 + n^3/6 - i*n*(P + 2) + 1;
                           int o113 = l + n*(2*P + P^2/2 + 11/6) - j*k - k*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - n^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*n^2)/2 + (j^2*n)/2 + k*(P + 3/2) - j^4/24 - k^2/2 + n^3/6 - j*n*(P + 2) + 1;
                           int o114 = l + i*(2*P + P^2/2 + 11/6) - i*n - j*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - i^2*(P/2 + 1) + j^3*(P/6 + 5/12) - j^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + n*(P + 3/2) + i^3/6 - j^4/24 - n^2/2 - i*j*(P + 2) + 1;
                           int o115 = l + i*(2*P + P^2/2 + 11/6) - i*n - k*n - i^2*(P/2 + 1) + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + n*(P + 3/2) + i^3/6 - k^4/24 - n^2/2 - i*k*(P + 2) + 1;
                           int o116 = l + j*(2*P + P^2/2 + 11/6) - i*n - j*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - j^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*j^2)/2 + (i^2*j)/2 + n*(P + 3/2) - i^4/24 + j^3/6 - n^2/2 - i*j*(P + 2) + 1;
                           int o117 = l + j*(2*P + P^2/2 + 11/6) - j*n - k*n + k*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) - j^2*(P/2 + 1) + k^3*(P/6 + 5/12) - k^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + n*(P + 3/2) + j^3/6 - k^4/24 - n^2/2 - j*k*(P + 2) + 1;
                           int o118 = l + k*(2*P + P^2/2 + 11/6) - i*n - k*n + i*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + i^3*(P/6 + 5/12) - k^2*(P/2 + 1) - i^2*((5*P)/4 + P^2/4 + 35/24) + (i*k^2)/2 + (i^2*k)/2 + n*(P + 3/2) - i^4/24 + k^3/6 - n^2/2 - i*k*(P + 2) + 1;
                           int o119 = l + k*(2*P + P^2/2 + 11/6) - j*n - k*n + j*((35*P)/12 + (5*P^2)/4 + P^3/6 + 25/12) + j^3*(P/6 + 5/12) - k^2*(P/2 + 1) - j^2*((5*P)/4 + P^2/4 + 35/24) + (j*k^2)/2 + (j^2*k)/2 + n*(P + 3/2) - j^4/24 + k^3/6 - n^2/2 - j*k*(P + 2) + 1;
                           

                           PentDofOrd[0][o] = o;
                           PentDofOrd[1][o] = o1;
                           PentDofOrd[2][o] = o2;
                           PentDofOrd[3][o] = o3;
                           PentDofOrd[4][o] = o4;
                           PentDofOrd[5][o] = o5;
                           PentDofOrd[6][o] = o6;
                           PentDofOrd[7][o] = o7;
                           PentDofOrd[8][o] = o8;
                           PentDofOrd[9][o] = o9;
                           PentDofOrd[10][o] = o10;
                           PentDofOrd[11][o] = o11;
                           PentDofOrd[12][o] = o12;
                           PentDofOrd[13][o] = o13;
                           PentDofOrd[14][o] = o14;
                           PentDofOrd[15][o] = o15;
                           PentDofOrd[16][o] = o16;
                           PentDofOrd[17][o] = o17;
                           PentDofOrd[18][o] = o18;
                           PentDofOrd[19][o] = o19;
                           PentDofOrd[20][o] = o20;
                           PentDofOrd[21][o] = o21;
                           PentDofOrd[22][o] = o22;
                           PentDofOrd[23][o] = o23;
                           PentDofOrd[24][o] = o24;
                           PentDofOrd[25][o] = o25;
                           PentDofOrd[26][o] = o26;
                           PentDofOrd[27][o] = o27;
                           PentDofOrd[28][o] = o28;
                           PentDofOrd[29][o] = o29;
                           PentDofOrd[30][o] = o30;
                           PentDofOrd[31][o] = o31;
                           PentDofOrd[32][o] = o32;
                           PentDofOrd[33][o] = o33;
                           PentDofOrd[34][o] = o34;
                           PentDofOrd[35][o] = o35;
                           PentDofOrd[36][o] = o36;
                           PentDofOrd[37][o] = o37;
                           PentDofOrd[38][o] = o38;
                           PentDofOrd[39][o] = o39;
                           PentDofOrd[40][o] = o40;
                           PentDofOrd[41][o] = o41;
                           PentDofOrd[42][o] = o42;
                           PentDofOrd[43][o] = o43;
                           PentDofOrd[44][o] = o44;
                           PentDofOrd[45][o] = o45;
                           PentDofOrd[46][o] = o46;
                           PentDofOrd[47][o] = o47;
                           PentDofOrd[48][o] = o48;
                           PentDofOrd[49][o] = o49;
                           PentDofOrd[50][o] = o50;
                           PentDofOrd[51][o] = o51;
                           PentDofOrd[52][o] = o52;
                           PentDofOrd[53][o] = o53;
                           PentDofOrd[54][o] = o54;
                           PentDofOrd[55][o] = o55;
                           PentDofOrd[56][o] = o56;
                           PentDofOrd[57][o] = o57;
                           PentDofOrd[58][o] = o58;
                           PentDofOrd[59][o] = o59;
                           PentDofOrd[60][o] = o60;
                           PentDofOrd[61][o] = o61;
                           PentDofOrd[62][o] = o62;
                           PentDofOrd[63][o] = o63;
                           PentDofOrd[64][o] = o64;
                           PentDofOrd[65][o] = o65;
                           PentDofOrd[66][o] = o66;
                           PentDofOrd[67][o] = o67;
                           PentDofOrd[68][o] = o68;
                           PentDofOrd[69][o] = o69;
                           PentDofOrd[70][o] = o70;
                           PentDofOrd[71][o] = o71;
                           PentDofOrd[72][o] = o72;
                           PentDofOrd[73][o] = o73;
                           PentDofOrd[74][o] = o74;
                           PentDofOrd[75][o] = o75;
                           PentDofOrd[76][o] = o76;
                           PentDofOrd[77][o] = o77;
                           PentDofOrd[78][o] = o78;
                           PentDofOrd[79][o] = o79;
                           PentDofOrd[80][o] = o80;
                           PentDofOrd[81][o] = o81;
                           PentDofOrd[82][o] = o82;
                           PentDofOrd[83][o] = o83;
                           PentDofOrd[84][o] = o84;
                           PentDofOrd[85][o] = o85;
                           PentDofOrd[86][o] = o86;
                           PentDofOrd[87][o] = o87;
                           PentDofOrd[88][o] = o88;
                           PentDofOrd[89][o] = o89;
                           PentDofOrd[90][o] = o90;
                           PentDofOrd[91][o] = o91;
                           PentDofOrd[92][o] = o92;
                           PentDofOrd[93][o] = o93;
                           PentDofOrd[94][o] = o94;
                           PentDofOrd[95][o] = o95;
                           PentDofOrd[96][o] = o96;
                           PentDofOrd[97][o] = o97;
                           PentDofOrd[98][o] = o98;
                           PentDofOrd[99][o] = o99;
                           PentDofOrd[100][o] = o100;
                           PentDofOrd[101][o] = o101;
                           PentDofOrd[102][o] = o102;
                           PentDofOrd[103][o] = o103;
                           PentDofOrd[104][o] = o104;
                           PentDofOrd[105][o] = o105;
                           PentDofOrd[106][o] = o106;
                           PentDofOrd[107][o] = o107;
                           PentDofOrd[108][o] = o108;
                           PentDofOrd[109][o] = o109;
                           PentDofOrd[110][o] = o110;
                           PentDofOrd[111][o] = o111;
                           PentDofOrd[112][o] = o112;
                           PentDofOrd[113][o] = o113;
                           PentDofOrd[114][o] = o114;
                           PentDofOrd[115][o] = o115;
                           PentDofOrd[116][o] = o116;
                           PentDofOrd[117][o] = o117;
                           PentDofOrd[118][o] = o118;
                           PentDofOrd[119][o] = o119;

                       }
                   }
                }
             }
             // end new code
             
         }
      }
   }
}

const int *H1_FECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                   int Or) const
{
   if (GeomType == Geometry::SEGMENT)
   {
      return (Or > 0) ? SegDofOrd[0] : SegDofOrd[1];
   }
   else if (GeomType == Geometry::TRIANGLE)
   {
      return TriDofOrd[Or%6];
   }
   else if (GeomType == Geometry::SQUARE)
   {
      return QuadDofOrd[Or%8];
   }
   else if (GeomType == Geometry::TETRAHEDRON)
   {
      return TetDofOrd[Or%24];
   }   
   else if (GeomType == Geometry::PENTATOPE)
   {
      return PentDofOrd[Or%120];
   }
   return NULL;
}

FiniteElementCollection *H1_FECollection::GetTraceCollection() const
{
   int p = H1_dof[Geometry::SEGMENT] + 1;
   int dim = -1;
   if (!strncmp(h1_name, "H1_", 3))
   {
      dim = atoi(h1_name + 3);
   }
   else if (!strncmp(h1_name, "H1Pos_", 6))
   {
      dim = atoi(h1_name + 6);
   }
   else if (!strncmp(h1_name, "H1@", 3))
   {
      dim = atoi(h1_name + 5);
   }
   return (dim < 0) ? NULL : new H1_Trace_FECollection(p, dim, b_type);
}

const int *H1_FECollection::GetDofMap(Geometry::Type GeomType) const
{
   const int *dof_map = NULL;
   const FiniteElement *fe = H1_Elements[GeomType];
   const NodalFiniteElement *nodal_fe =
      dynamic_cast<const NodalFiniteElement*>(fe);
   if (nodal_fe)
   {
      dof_map = nodal_fe->GetLexicographicOrdering().GetData();
   }
   else
   {
      MFEM_ABORT("Geometry type " << Geometry::Name[GeomType] << " is not "
                 "implemented");
   }
   return dof_map;
}

const int *H1_FECollection::GetDofMap(Geometry::Type GeomType, int p) const
{
   if (p == base_p) { return GetDofMap(GeomType); }
   if (p >= var_orders.Size() || !var_orders[p]) { InitVarOrder(p); }
   return ((H1_FECollection*) var_orders[p])->GetDofMap(GeomType);
}

H1_FECollection::~H1_FECollection()
{
   delete [] SegDofOrd[0];
   delete [] TriDofOrd[0];
   delete [] QuadDofOrd[0];
   delete [] TetDofOrd[0];
   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      delete H1_Elements[g];
   }
}


H1_Trace_FECollection::H1_Trace_FECollection(const int p, const int dim,
                                             const int btype)
   : H1_FECollection(p, dim-1, btype)
{
   if (btype == BasisType::GaussLobatto)
   {
      snprintf(h1_name, 32, "H1_Trace_%dD_P%d", dim, p);
   }
   else if (btype == BasisType::Positive)
   {
      snprintf(h1_name, 32, "H1Pos_Trace_%dD_P%d", dim, p);
   }
   else // base class checks that type is closed
   {
      snprintf(h1_name, 32, "H1_Trace@%c_%dD_P%d",
               (int)BasisType::GetChar(btype), dim, p);
   }
}


L2_FECollection::L2_FECollection(const int p, const int dim, const int btype,
                                 const int map_type)
   : FiniteElementCollection(p)
   , dim(dim)
   , m_type(map_type)
{
   MFEM_VERIFY(p >= 0, "L2_FECollection requires order >= 0.");

   b_type = BasisType::Check(btype);
   const char *prefix = NULL;
   switch (map_type)
   {
      case FiniteElement::VALUE:    prefix = "L2";    break;
      case FiniteElement::INTEGRAL: prefix = "L2Int"; break;
      default:
         MFEM_ABORT("invalid map_type: " << map_type);
   }
   switch (btype)
   {
      case BasisType::GaussLegendre:
         snprintf(d_name, 32, "%s_%dD_P%d", prefix, dim, p);
         break;
      default:
         snprintf(d_name, 32, "%s_T%d_%dD_P%d", prefix, btype, dim, p);
   }

   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      L2_Elements[g] = NULL;
      Tr_Elements[g] = NULL;
   }
   for (int i = 0; i < 2; i++)
   {
      SegDofOrd[i] = NULL;
   }
   for (int i = 0; i < 6; i++)
   {
      TriDofOrd[i] = NULL;
   }
   for (int i = 0; i < 24; i++)
   {
      TetDofOrd[i] = NULL;
   }
   OtherDofOrd = NULL;

   if (dim == 0)
   {
      L2_Elements[Geometry::POINT] = new PointFiniteElement;
   }
   else if (dim == 1)
   {
      if (b_type == BasisType::Positive)
      {
         L2_Elements[Geometry::SEGMENT] = new L2Pos_SegmentElement(p);
      }
      else
      {
         L2_Elements[Geometry::SEGMENT] = new L2_SegmentElement(p, btype);
      }
      L2_Elements[Geometry::SEGMENT]->SetMapType(map_type);

      Tr_Elements[Geometry::POINT] = new PointFiniteElement;
      // No need to set the map_type for Tr_Elements.

      const int pp1 = p + 1;
      SegDofOrd[0] = new int[2*pp1];
      SegDofOrd[1] = SegDofOrd[0] + pp1;
      for (int i = 0; i <= p; i++)
      {
         SegDofOrd[0][i] = i;
         SegDofOrd[1][i] = p - i;
      }
   }
   else if (dim == 2)
   {
      if (b_type == BasisType::Positive)
      {
         L2_Elements[Geometry::TRIANGLE] = new L2Pos_TriangleElement(p);
         L2_Elements[Geometry::SQUARE] = new L2Pos_QuadrilateralElement(p);
      }
      else
      {
         L2_Elements[Geometry::TRIANGLE] = new L2_TriangleElement(p, btype);
         L2_Elements[Geometry::SQUARE] = new L2_QuadrilateralElement(p, btype);
      }
      L2_Elements[Geometry::TRIANGLE]->SetMapType(map_type);
      L2_Elements[Geometry::SQUARE]->SetMapType(map_type);
      // Trace element use the default Gauss-Legendre nodal points for positive basis
      if (b_type == BasisType::Positive)
      {
         Tr_Elements[Geometry::SEGMENT] = new L2Pos_SegmentElement(p);
      }
      else
      {
         Tr_Elements[Geometry::SEGMENT] = new L2_SegmentElement(p, btype);
      }

      const int TriDof = L2_Elements[Geometry::TRIANGLE]->GetDof();
      TriDofOrd[0] = new int[6*TriDof];
      for (int i = 1; i < 6; i++)
      {
         TriDofOrd[i] = TriDofOrd[i-1] + TriDof;
      }
      const int pp1 = p + 1, pp2 = pp1 + 1;
      for (int j = 0; j <= p; j++)
      {
         for (int i = 0; i + j <= p; i++)
         {
            int o = TriDof - ((pp2 - j)*(pp1 - j))/2 + i;
            int k = p - j - i;
            TriDofOrd[0][o] = o;  // (0,1,2)
            TriDofOrd[1][o] = TriDof - ((pp2-j)*(pp1-j))/2 + k;  // (1,0,2)
            TriDofOrd[2][o] = TriDof - ((pp2-i)*(pp1-i))/2 + k;  // (2,0,1)
            TriDofOrd[3][o] = TriDof - ((pp2-k)*(pp1-k))/2 + i;  // (2,1,0)
            TriDofOrd[4][o] = TriDof - ((pp2-k)*(pp1-k))/2 + j;  // (1,2,0)
            TriDofOrd[5][o] = TriDof - ((pp2-i)*(pp1-i))/2 + j;  // (0,2,1)
         }
      }
      const int QuadDof = L2_Elements[Geometry::SQUARE]->GetDof();
      OtherDofOrd = new int[QuadDof];
      for (int j = 0; j < QuadDof; j++)
      {
         OtherDofOrd[j] = j; // for Or == 0
      }
   }
   else if (dim == 3)
   {
      if (b_type == BasisType::Positive)
      {
         L2_Elements[Geometry::TETRAHEDRON] = new L2Pos_TetrahedronElement(p);
         L2_Elements[Geometry::CUBE] = new L2Pos_HexahedronElement(p);
         L2_Elements[Geometry::PRISM] = new L2Pos_WedgeElement(p);
      }
      else
      {
         L2_Elements[Geometry::TETRAHEDRON] =
            new L2_TetrahedronElement(p, btype);
         L2_Elements[Geometry::CUBE] = new L2_HexahedronElement(p, btype);
         L2_Elements[Geometry::PRISM] = new L2_WedgeElement(p, btype);
      }
      L2_Elements[Geometry::TETRAHEDRON]->SetMapType(map_type);
      L2_Elements[Geometry::CUBE]->SetMapType(map_type);
      L2_Elements[Geometry::PRISM]->SetMapType(map_type);
      // Trace element use the default Gauss-Legendre nodal points for positive basis
      if (b_type == BasisType::Positive)
      {
         Tr_Elements[Geometry::TRIANGLE] = new L2Pos_TriangleElement(p);
         Tr_Elements[Geometry::SQUARE] = new L2Pos_QuadrilateralElement(p);
      }
      else
      {
         Tr_Elements[Geometry::TRIANGLE] = new L2_TriangleElement(p, btype);
         Tr_Elements[Geometry::SQUARE] = new L2_QuadrilateralElement(p, btype);
      }

      const int TetDof = L2_Elements[Geometry::TETRAHEDRON]->GetDof();
      const int HexDof = L2_Elements[Geometry::CUBE]->GetDof();
      const int PriDof = L2_Elements[Geometry::PRISM]->GetDof();
      const int MaxDof = std::max(TetDof, std::max(PriDof, HexDof));

      TetDofOrd[0] = new int[24*TetDof];
      for (int i = 1; i < 24; i++)
      {
         TetDofOrd[i] = TetDofOrd[i-1] + TetDof;
      }
      // see Mesh::GetTetOrientation in mesh/mesh.cpp
      const int pp1 = p + 1, pp2 = pp1 + 1, pp3 = pp2 + 1;
      for (int k = 0; k <= p; k++)
      {
         for (int j = 0; j + k <= p; j++)
         {
            for (int i = 0; i + j + k <= p; i++)
            {
               int l = p - k - j - i;
               int o   = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (j * (2 * p + 3 - j - 2 * k)) / 2 + i;
               int o1  = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (k * (2 * p + 3 - k - 2 * j)) / 2 + i;
               int o2  = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (k * (2 * p + 3 - k - 2 * i)) / 2 + j;
               int o3  = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (i * (2 * p + 3 - i - 2 * k)) / 2 + j;
               int o4  = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (i * (2 * p + 3 - i - 2 * j)) / 2 + k;
               int o5  = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (j * (2 * p + 3 - j - 2 * i)) / 2 + k;
               int o6  = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (l * (2 * p + 3 - l - 2 * k)) / 2 + j;
               int o7  = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (k * (2 * p + 3 - k - 2 * l)) / 2 + j;
               int o8  = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (j * (2 * p + 3 - j - 2 * l)) / 2 + k;
               int o9  = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (l * (2 * p + 3 - l - 2 * j)) / 2 + k;
               int o10 = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (k * (2 * p + 3 - k - 2 * j)) / 2 + l;
               int o11 = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (j * (2 * p + 3 - j - 2 * k)) / 2 + l;
               int o12 = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (l * (2 * p + 3 - l - 2 * i)) / 2 + k;
               int o13 = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (i * (2 * p + 3 - i - 2 * l)) / 2 + k;
               int o14 = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (i * (2 * p + 3 - i - 2 * k)) / 2 + l;
               int o15 = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (k * (2 * p + 3 - k - 2 * i)) / 2 + l;
               int o16 = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (k * (2 * p + 3 - k - 2 * l)) / 2 + i;
               int o17 = TetDof - ((pp1 - k) * (pp2 - k) * (pp3 - k)) / 6
                         + (l * (2 * p + 3 - l - 2 * k)) / 2 + i;
               int o18 = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (j * (2 * p + 3 - j - 2 * i)) / 2 + l;
               int o19 = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (i * (2 * p + 3 - i - 2 * j)) / 2 + l;
               int o20 = TetDof - ((pp1 - j) * (pp2 - j) * (pp3 - j)) / 6
                         + (l * (2 * p + 3 - l - 2 * j)) / 2 + i;
               int o21 = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (j * (2 * p + 3 - j - 2 * l)) / 2 + i;
               int o22 = TetDof - ((pp1 - l) * (pp2 - l) * (pp3 - l)) / 6
                         + (i * (2 * p + 3 - i - 2 * l)) / 2 + j;
               int o23 = TetDof - ((pp1 - i) * (pp2 - i) * (pp3 - i)) / 6
                         + (l * (2 * p + 3 - l - 2 * i)) / 2 + j;
               TetDofOrd[ 0][o] = o;   // (0,1,2,3)
               TetDofOrd[ 1][o] = o1;  // (0,1,3,2)
               TetDofOrd[ 2][o] = o2;  // (0,2,3,1)
               TetDofOrd[ 3][o] = o3;  // (0,2,1,3)
               TetDofOrd[ 4][o] = o4;  // (0,3,1,2)
               TetDofOrd[ 5][o] = o5;  // (0,3,2,1)
               TetDofOrd[ 6][o] = o6;  // (1,2,0,3)
               TetDofOrd[ 7][o] = o7;  // (1,2,3,0)
               TetDofOrd[ 8][o] = o8;  // (1,3,2,0)
               TetDofOrd[ 9][o] = o9;  // (1,3,0,2)
               TetDofOrd[10][o] = o10; // (1,0,3,2)
               TetDofOrd[11][o] = o11; // (1,0,2,3)
               TetDofOrd[12][o] = o12; // (2,3,0,1)
               TetDofOrd[13][o] = o13; // (2,3,1,0)
               TetDofOrd[14][o] = o14; // (2,0,1,3)
               TetDofOrd[15][o] = o15; // (2,0,3,1)
               TetDofOrd[16][o] = o16; // (2,1,3,0)
               TetDofOrd[17][o] = o17; // (2,1,0,3)
               TetDofOrd[18][o] = o18; // (3,0,2,1)
               TetDofOrd[19][o] = o19; // (3,0,1,2)
               TetDofOrd[20][o] = o20; // (3,1,0,2)
               TetDofOrd[21][o] = o21; // (3,1,2,0)
               TetDofOrd[22][o] = o22; // (3,2,1,0)
               TetDofOrd[23][o] = o23; // (3,2,0,1)
            }
         }
      }
      OtherDofOrd = new int[MaxDof];
      for (int j = 0; j < MaxDof; j++)
      {
         OtherDofOrd[j] = j; // for Or == 0
      }
   }
   else if (dim == 4)
   {
      if (b_type == BasisType::Positive)
      {
         mfem::err <<
                   "L2_FECollection::L2_FECollection : BasisType::Positive not implemented" <<
                   endl;
         mfem_error();
      }
      else
      {
         L2_Elements[Geometry::PENTATOPE] =
            new L2_PentatopeElement(p, btype);
         L2_Elements[Geometry::TESSERACT] = new L2_HexahedronElement(p, btype);
      }
      L2_Elements[Geometry::PENTATOPE]->SetMapType(map_type);
      L2_Elements[Geometry::TESSERACT]->SetMapType(map_type);
      // All trace element use the default Gauss-Legendre nodal points
      Tr_Elements[Geometry::TETRAHEDRON] = new L2_TetrahedronElement(p);
      Tr_Elements[Geometry::CUBE] = new L2_HexahedronElement(p);

      const int PentDof = L2_Elements[Geometry::PENTATOPE]->GetDof();
      const int TessDof = L2_Elements[Geometry::TESSERACT]->GetDof();
      const int MaxDof = std::max(PentDof, TessDof);
      OtherDofOrd = new int[MaxDof];
      for (int j = 0; j < MaxDof; j++)
      {
         OtherDofOrd[j] = j; // for Or == 0
      }
   }
   else
   {
      mfem::err << "L2_FECollection::L2_FECollection : dim = "
                << dim << endl;
      mfem_error();
   }
}

const int *L2_FECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                   int Or) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:
         return (Or > 0) ? SegDofOrd[0] : SegDofOrd[1];

      case Geometry::TRIANGLE:
         return TriDofOrd[Or%6];

      case Geometry::TETRAHEDRON:
         return TetDofOrd[Or%24];

      default:
         return (Or == 0) ? OtherDofOrd : NULL;
   }
}

L2_FECollection::~L2_FECollection()
{
   delete [] OtherDofOrd;
   delete [] SegDofOrd[0];
   delete [] TriDofOrd[0];
   delete [] TetDofOrd[0];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      delete L2_Elements[i];
      delete Tr_Elements[i];
   }
}


RT_FECollection::RT_FECollection(const int order, const int dim,
                                 const int cb_type, const int ob_type)
   : FiniteElementCollection(order + 1)
   , ob_type(ob_type)
{
   int p = order;
   MFEM_VERIFY(p >= 0, "RT_FECollection requires order >= 0.");

   int cp_type = BasisType::GetQuadrature1D(cb_type);
   int op_type = BasisType::GetQuadrature1D(ob_type);

   if (Quadrature1D::CheckClosed(cp_type) == Quadrature1D::Invalid)
   {
      const char *cb_name = BasisType::Name(cb_type); // this may abort
      MFEM_ABORT("unknown closed BasisType: " << cb_name);
   }
   if (Quadrature1D::CheckOpen(op_type) == Quadrature1D::Invalid)
   {
      const char *ob_name = BasisType::Name(ob_type); // this may abort
      MFEM_ABORT("unknown open BasisType: " << ob_name);
   }

   InitFaces(p, dim, FiniteElement::INTEGRAL, true);

   if (cb_type == BasisType::GaussLobatto &&
       ob_type == BasisType::GaussLegendre)
   {
      snprintf(rt_name, 32, "RT_%dD_P%d", dim, p);
   }
   else
   {
      snprintf(rt_name, 32, "RT@%c%c_%dD_P%d", (int)BasisType::GetChar(cb_type),
               (int)BasisType::GetChar(ob_type), dim, p);
   }

   const int pp1 = p + 1;
   if (dim == 2)
   {
      // TODO: cb_type, ob_type for triangles
      RT_Elements[Geometry::TRIANGLE] = new RT_TriangleElement(p);
      RT_dof[Geometry::TRIANGLE] = p*pp1;

      RT_Elements[Geometry::SQUARE] = new RT_QuadrilateralElement(p, cb_type,
                                                                  ob_type);
      // two vector components * n_unk_face *
      RT_dof[Geometry::SQUARE] = 2*p*pp1;
   }
   else if (dim == 3)
   {
      // TODO: cb_type, ob_type for tets
      RT_Elements[Geometry::TETRAHEDRON] = new RT_TetrahedronElement(p);
      RT_dof[Geometry::TETRAHEDRON] = p*pp1*(p + 2)/2;

      RT_Elements[Geometry::CUBE] = new RT_HexahedronElement(p, cb_type, ob_type);
      RT_dof[Geometry::CUBE] = 3*p*pp1*pp1;
   }
   else if (dim == 4)
   {
      //RT_Elements[Geometry::PENTATOPE] = new RT_PentatopeElement(p);
      //RT_Elements[Geometry::PENTATOPE] = new RT0PentFiniteElement();
      RT_Elements[Geometry::PENTATOPE] = new Hdiv_PentatopeElement(p);
       
      RT_dof[Geometry::PENTATOPE] = p*pp1*(p + 2)*(p + 3)/6;

      //TODO: tesseracts
   }
   else
   {
      MFEM_ABORT("invalid dim = " << dim);
   }
}

// This is a special protected constructor only used by RT_Trace_FECollection
// and DG_Interface_FECollection
RT_FECollection::RT_FECollection(const int order, const int dim,
                                 const int map_type, const bool signs,
                                 const int ob_type)
   : FiniteElementCollection(order + 1)
   , ob_type(ob_type)
{
   if (Quadrature1D::CheckOpen(BasisType::GetQuadrature1D(ob_type)) ==
       Quadrature1D::Invalid)
   {
      const char *ob_name = BasisType::Name(ob_type); // this may abort
      MFEM_ABORT("Invalid open basis type: " << ob_name);
   }
   InitFaces(order, dim, map_type, signs);
}

void RT_FECollection::InitFaces(const int order, const int dim,
                                const int map_type,
                                const bool signs)
{
   int p = order;
   int op_type = BasisType::GetQuadrature1D(ob_type);

   MFEM_VERIFY(Quadrature1D::CheckOpen(op_type) != Quadrature1D::Invalid,
               "invalid open point type");

   const int pp1 = p + 1, pp2 = p + 2, pp3 = p + 3;

   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      RT_Elements[g] = NULL;
      RT_dof[g] = 0;
   }
   // Degree of Freedom orderings
   for (int i = 0; i < 2; i++)
   {
      SegDofOrd[i] = NULL;
   }
   for (int i = 0; i < 6; i++)
   {
      TriDofOrd[i] = NULL;
   }
   for (int i = 0; i < 8; i++)
   {
      QuadDofOrd[i] = NULL;
   }
   for (int i = 0; i < 24; i++)
   {
      TetDofOrd[i] = NULL;
   }

   if (dim == 2)
   {
      L2_SegmentElement *l2_seg = new L2_SegmentElement(p, ob_type);
      l2_seg->SetMapType(map_type);
      RT_Elements[Geometry::SEGMENT] = l2_seg;
      RT_dof[Geometry::SEGMENT] = pp1;

      SegDofOrd[0] = new int[2*pp1];
      SegDofOrd[1] = SegDofOrd[0] + pp1;
      for (int i = 0; i <= p; i++)
      {
         SegDofOrd[0][i] = i;
         SegDofOrd[1][i] = signs ? (-1 - (p - i)) : (p - i);
      }
   }
   else if (dim == 3)
   {
      L2_TriangleElement *l2_tri = new L2_TriangleElement(p, ob_type);
      l2_tri->SetMapType(map_type);
      RT_Elements[Geometry::TRIANGLE] = l2_tri;
      RT_dof[Geometry::TRIANGLE] = pp1*pp2/2;

      L2_QuadrilateralElement *l2_quad = new L2_QuadrilateralElement(p, ob_type);
      l2_quad->SetMapType(map_type);
      RT_Elements[Geometry::SQUARE] = l2_quad;
      RT_dof[Geometry::SQUARE] = pp1*pp1;

      int TriDof = RT_dof[Geometry::TRIANGLE];
      TriDofOrd[0] = new int[6*TriDof];
      for (int i = 1; i < 6; i++)
      {
         TriDofOrd[i] = TriDofOrd[i-1] + TriDof;
      }
      // see Mesh::GetTriOrientation in mesh/mesh.cpp,
      // the constructor of H1_FECollection
      for (int j = 0; j <= p; j++)
      {
         for (int i = 0; i + j <= p; i++)
         {
            int o = TriDof - ((pp2 - j)*(pp1 - j))/2 + i;
            int k = p - j - i;
            TriDofOrd[0][o] = o; // (0,1,2)
             std::cout << "o = " << o << std::endl;
            TriDofOrd[1][o] = -1-(TriDof-((pp2-j)*(pp1-j))/2+k);  // (1,0,2)
            TriDofOrd[2][o] =     TriDof-((pp2-i)*(pp1-i))/2+k;   // (2,0,1)
            TriDofOrd[3][o] = -1-(TriDof-((pp2-k)*(pp1-k))/2+i);  // (2,1,0)
            TriDofOrd[4][o] =     TriDof-((pp2-k)*(pp1-k))/2+j;   // (1,2,0)
            TriDofOrd[5][o] = -1-(TriDof-((pp2-i)*(pp1-i))/2+j);  // (0,2,1)
//             for (int ii=0; ii<6; ii++) {
//                 std::cout << "TriDofOrd = " << TriDofOrd[ii][o] << std::endl;
//             }
            if (!signs)
            {
               for (int k = 1; k < 6; k += 2)
               {
                  TriDofOrd[k][o] = -1 - TriDofOrd[k][o];
                   std::cout << "Inside IF statement" << std::endl;
               }
            }
         }
      }

      int QuadDof = RT_dof[Geometry::SQUARE];
      QuadDofOrd[0] = new int[8*QuadDof];
      for (int i = 1; i < 8; i++)
      {
         QuadDofOrd[i] = QuadDofOrd[i-1] + QuadDof;
      }
      // see Mesh::GetQuadOrientation in mesh/mesh.cpp
      for (int j = 0; j <= p; j++)
      {
         for (int i = 0; i <= p; i++)
         {
            int o = i + j*pp1;
            QuadDofOrd[0][o] = i + j*pp1;                    // (0,1,2,3)
            QuadDofOrd[1][o] = -1 - (j + i*pp1);             // (0,3,2,1)
            QuadDofOrd[2][o] = j + (p - i)*pp1;              // (1,2,3,0)
            QuadDofOrd[3][o] = -1 - ((p - i) + j*pp1);       // (1,0,3,2)
            QuadDofOrd[4][o] = (p - i) + (p - j)*pp1;        // (2,3,0,1)
            QuadDofOrd[5][o] = -1 - ((p - j) + (p - i)*pp1); // (2,1,0,3)
            QuadDofOrd[6][o] = (p - j) + i*pp1;              // (3,0,1,2)
            QuadDofOrd[7][o] = -1 - (i + (p - j)*pp1);       // (3,2,1,0)
            if (!signs)
            {
               for (int k = 1; k < 8; k += 2)
               {
                  QuadDofOrd[k][o] = -1 - QuadDofOrd[k][o];
               }
            }
         }
      }
   }
   else if (dim == 4)
   {
      L2_TetrahedronElement *l2_tet = new L2_TetrahedronElement(p, ob_type);
      l2_tet->SetMapType(map_type);
      RT_Elements[Geometry::TETRAHEDRON] = l2_tet;
      RT_dof[Geometry::TETRAHEDRON] = pp1*pp2*pp3/6;

      int TetDof = RT_dof[Geometry::TETRAHEDRON];
      //std::cout << "TetDof = " << TetDof << " p = " << p << std::endl;
      int TriDof2 = pp2*pp1/2;
      TetDofOrd[0] = new int[24*TetDof];
      for (int i = 1; i < 24; i++)
      {
         TetDofOrd[i] = TetDofOrd[i-1] + TetDof;
      }
      // see Mesh::GetTriOrientation in mesh/mesh.cpp,
      // the constructor of H1_FECollection
      for (int k=0; k<=p; k++)
      {
         for (int j=0; j+k<=p; j++)
         {
            for (int i=0; i+j+k<=p; i++)
            {
               int o = TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 - (pp2-j)*
                       (pp1-j)/2 - k*j + i;
               //std::cout << "o = " << o << std::endl;
               //std::cout << "i = " << i << "j = " << j << "k = " << k << std::endl;

               int l = p-k-j-i;
               TetDofOrd[0][o] = o;
                
               TetDofOrd[1][o] = -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-k)*(pp1-k)/2 - j*k + i);
                
               TetDofOrd[2][o] =       TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-k)*(pp1-k)/2 - i*k + j;
                
               TetDofOrd[3][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                         (pp2-i)*(pp1-i)/2 - k*i + j);
                
               TetDofOrd[4][o] =        TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-i)*(pp1-i)/2 - j*i + k;
                
               TetDofOrd[5][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-j)*(pp1-j)/2 - i*j + k);
                
               TetDofOrd[6][o] =        TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                         (pp2-l)*(pp1-l)/2 - k*l + j;
                
               TetDofOrd[7][o] = -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-k)*(pp1-k)/2 - l*k + j);
                
               TetDofOrd[8][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-j)*(pp1-j)/2 - l*j + k;
                
               TetDofOrd[9][o] =  -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-l)*(pp1-l)/2 - j*l + k);
                
               TetDofOrd[10][o] =       TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-k)*(pp1-k)/2 - j*k + l;
            
               TetDofOrd[11][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                        (pp2-j)*(pp1-j)/2 - k*j + l);
                
               TetDofOrd[12][o] =        TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-l)*(pp1-l)/2 - i*l + k;
                
               TetDofOrd[13][o] =  -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-i)*(pp1-i)/2 - l*i + k);
                
               TetDofOrd[14][o] =        TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                        (pp2-i)*(pp1-i)/2 - k*i + l;
                
               TetDofOrd[15][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-k)*(pp1-k)/2 - i*k + l);
                
               TetDofOrd[16][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-k)*(pp1-k)/2 - l*k + i;
                
               TetDofOrd[17][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                        (pp2-l)*(pp1-l)/2 - k*l + i);
                
               TetDofOrd[18][o] =       TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-j)*(pp1-j)/2 - i*j + l;
                
               TetDofOrd[19][o] = -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-i)*(pp1-i)/2 - j*i + l);
                
               TetDofOrd[20][o] =       TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                         (pp2-l)*(pp1-l)/2 - j*l + i;
                
               TetDofOrd[21][o] = -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-j)*(pp1-j)/2 - l*j + i);
                
               TetDofOrd[22][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                         (pp2-i)*(pp1-i)/2 - l*i + j;
                
               TetDofOrd[23][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                         (pp2-l)*(pp1-l)/2 - i*l + j);
                
               // OLD ORDER
                /*TetDofOrd[0][o] = o;
                 TetDofOrd[1][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                          (pp2-j)*(pp1-j)/2 - k*j + l);
                 TetDofOrd[2][o] =        TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                          (pp2-i)*(pp1-i)/2 - k*i + l;
                 TetDofOrd[3][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                          (pp2-l)*(pp1-l)/2 - k*l + i);
                 TetDofOrd[4][o] =        TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                          (pp2-l)*(pp1-l)/2 - k*l + j;
                 TetDofOrd[5][o] =  -1 - (TetDof + TriDof2 - ((pp3-k)*(pp2-k)*(pp1-k))/6 -
                                          (pp2-i)*(pp1-i)/2 - k*i + j);
                 TetDofOrd[6][o] =        TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-i)*(pp1-i)/2 - j*i + k;
                 TetDofOrd[7][o] =  -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-l)*(pp1-l)/2 - j*l + k);
                 TetDofOrd[8][o] =        TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-l)*(pp1-l)/2 - i*l + k;
                 TetDofOrd[9][o] =  -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-i)*(pp1-i)/2 - l*i + k);
                 TetDofOrd[10][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-j)*(pp1-j)/2 - l*j + k;
                 TetDofOrd[11][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-j)*(pp1-j)/2 - i*j + k);
                 TetDofOrd[12][o] =       TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-k)*(pp1-k)/2 - i*k + j;
                 TetDofOrd[13][o] = -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-k)*(pp1-k)/2 - l*k + j);
                 TetDofOrd[14][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-k)*(pp1-k)/2 - l*k + i;
                 TetDofOrd[15][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-k)*(pp1-k)/2 - i*k + l);
                 TetDofOrd[16][o] =       TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-k)*(pp1-k)/2 - j*k + l;
                 TetDofOrd[17][o] = -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-k)*(pp1-k)/2 - j*k + i);
                 TetDofOrd[18][o] =       TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-l)*(pp1-l)/2 - j*l + i;
                 TetDofOrd[19][o] = -1 - (TetDof + TriDof2 - ((pp3-j)*(pp2-j)*(pp1-j))/6 -
                                          (pp2-i)*(pp1-i)/2 - j*i + l);
                 TetDofOrd[20][o] =       TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-j)*(pp1-j)/2 - i*j + l;
                 TetDofOrd[21][o] = -1 - (TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-j)*(pp1-j)/2 - l*j + i);
                 TetDofOrd[22][o] =       TetDof + TriDof2 - ((pp3-l)*(pp2-l)*(pp1-l))/6 -
                                          (pp2-i)*(pp1-i)/2 - l*i + j;
                 TetDofOrd[23][o] = -1 - (TetDof + TriDof2 - ((pp3-i)*(pp2-i)*(pp1-i))/6 -
                                          (pp2-l)*(pp1-l)/2 - i*l + j);*/
//                for (int ii=0; ii<24; ii++) {
//                    std::cout << "TetDofOrd = " << TetDofOrd[ii][o] << std::endl;
//                }
               if (!signs)
               {
                  for (int m = 0; m < 24; m+=2)
                  {
                     TetDofOrd[m][o] = -1 - TetDofOrd[m][o];
                     std::cout << "Inside IF statement" << std::endl;

                  }
               }
            }
         }
      }
   }
}

const int *RT_FECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                   int Or) const
{
   if (GeomType == Geometry::SEGMENT)
   {
      return (Or > 0) ? SegDofOrd[0] : SegDofOrd[1];
   }
   else if (GeomType == Geometry::TRIANGLE)
   {
      return TriDofOrd[Or%6];
   }
   else if (GeomType == Geometry::SQUARE)
   {
      return QuadDofOrd[Or%8];
   }
   else if (GeomType == Geometry::TETRAHEDRON)
   {
      return TetDofOrd[Or%24];
   }
   return NULL;
}

FiniteElementCollection *RT_FECollection::GetTraceCollection() const
{
   int dim, p;
   if (!strncmp(rt_name, "RT_", 3))
   {
      dim = atoi(rt_name + 3);
      p = atoi(rt_name + 7);
   }
   else // rt_name = RT@.._.D_P*
   {
      dim = atoi(rt_name + 6);
      p = atoi(rt_name + 10);
   }
   return new RT_Trace_FECollection(p, dim, FiniteElement::INTEGRAL, ob_type);
}

RT_FECollection::~RT_FECollection()
{
   delete [] SegDofOrd[0];
   delete [] TriDofOrd[0];
   delete [] QuadDofOrd[0];
   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      delete RT_Elements[g];
   }
}


RT_Trace_FECollection::RT_Trace_FECollection(const int p, const int dim,
                                             const int map_type,
                                             const int ob_type)
   : RT_FECollection(p, dim, map_type, true, ob_type)
{
   const char *prefix =
      (map_type == FiniteElement::INTEGRAL) ? "RT_Trace" : "RT_ValTrace";
   char ob_str[3] = { '\0', '\0', '\0' };

   if (ob_type != BasisType::GaussLegendre)
   {
      ob_str[0] = '@';
      ob_str[1] = BasisType::GetChar(ob_type);
   }
   snprintf(rt_name, 32, "%s%s_%dD_P%d", prefix, ob_str, dim, p);

   MFEM_VERIFY(dim == 2 || dim == 3, "Wrong dimension, dim = " << dim);
}


DG_Interface_FECollection::DG_Interface_FECollection(const int p, const int dim,
                                                     const int map_type,
                                                     const int ob_type)
   : RT_FECollection(p, dim, map_type, false, ob_type)
{
   MFEM_VERIFY(dim == 2 || dim == 3, "Wrong dimension, dim = " << dim);

   const char *prefix =
      (map_type == FiniteElement::VALUE) ? "DG_Iface" : "DG_IntIface";
   if (ob_type == BasisType::GaussLegendre)
   {
      snprintf(rt_name, 32, "%s_%dD_P%d", prefix, dim, p);
   }
   else
   {
      snprintf(rt_name, 32, "%s@%c_%dD_P%d", prefix,
               (int)BasisType::GetChar(ob_type), dim, p);
   }
}

ND_FECollection::ND_FECollection(const int p, const int dim,
                                 const int cb_type, const int ob_type)
   : FiniteElementCollection(p)
{
   MFEM_VERIFY(p >= 1, "ND_FECollection requires order >= 1.");
   MFEM_VERIFY(dim >= 1 && dim <= 3, "ND_FECollection requires 1 <= dim <= 3.");

   const int pm1 = p - 1, pm2 = p - 2;

   if (cb_type == BasisType::GaussLobatto &&
       ob_type == BasisType::GaussLegendre)
   {
      snprintf(nd_name, 32, "ND_%dD_P%d", dim, p);
   }
   else
   {
      snprintf(nd_name, 32, "ND@%c%c_%dD_P%d", (int)BasisType::GetChar(cb_type),
               (int)BasisType::GetChar(ob_type), dim, p);
   }

   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      ND_Elements[g] = NULL;
      ND_dof[g] = 0;
   }
   for (int i = 0; i < 2; i++)
   {
      SegDofOrd[i] = NULL;
   }
   for (int i = 0; i < 6; i++)
   {
      TriDofOrd[i] = NULL;
   }
   for (int i = 0; i < 8; i++)
   {
      QuadDofOrd[i] = NULL;
   }

   int op_type = BasisType::GetQuadrature1D(ob_type);
   int cp_type = BasisType::GetQuadrature1D(cb_type);

   // Error checking
   if (Quadrature1D::CheckOpen(op_type) == Quadrature1D::Invalid)
   {
      const char *ob_name = BasisType::Name(ob_type);
      MFEM_ABORT("Invalid open basis point type: " << ob_name);
   }
   if (Quadrature1D::CheckClosed(cp_type) == Quadrature1D::Invalid)
   {
      const char *cb_name = BasisType::Name(cb_type);
      MFEM_ABORT("Invalid closed basis point type: " << cb_name);
   }

   if (dim >= 1)
   {
      ND_Elements[Geometry::SEGMENT] = new ND_SegmentElement(p, ob_type);
      ND_dof[Geometry::SEGMENT] = p;

      SegDofOrd[0] = new int[2*p];
      SegDofOrd[1] = SegDofOrd[0] + p;
      for (int i = 0; i < p; i++)
      {
         SegDofOrd[0][i] = i;
         SegDofOrd[1][i] = -1 - (pm1 - i);
      }
   }

   if (dim >= 2)
   {
      ND_Elements[Geometry::SQUARE] = new ND_QuadrilateralElement(p, cb_type,
                                                                  ob_type);
      ND_dof[Geometry::SQUARE] = 2*p*pm1;

      // TODO: cb_type and ob_type for triangles
      ND_Elements[Geometry::TRIANGLE] = new ND_TriangleElement(p);
      ND_dof[Geometry::TRIANGLE] = p*pm1;

      int QuadDof = ND_dof[Geometry::SQUARE];
      QuadDofOrd[0] = new int[8*QuadDof];
      for (int i = 1; i < 8; i++)
      {
         QuadDofOrd[i] = QuadDofOrd[i-1] + QuadDof;
      }
      // see Mesh::GetQuadOrientation in mesh/mesh.cpp
      for (int j = 0; j < pm1; j++)
      {
         for (int i = 0; i < p; i++)
         {
            int d1 = i + j*p;            // x-component
            int d2 = p*pm1 + j + i*pm1;  // y-component
            // (0,1,2,3)
            QuadDofOrd[0][d1] = d1;
            QuadDofOrd[0][d2] = d2;
            // (0,3,2,1)
            QuadDofOrd[1][d1] = d2;
            QuadDofOrd[1][d2] = d1;
            // (1,2,3,0)
            // QuadDofOrd[2][d1] = p*pm1 + (pm2 - j) + i*pm1;
            // QuadDofOrd[2][d2] = -1 - ((pm1 - i) + j*p);
            QuadDofOrd[2][d1] = -1 - (p*pm1 + j + (pm1 - i)*pm1);
            QuadDofOrd[2][d2] = i + (pm2 - j)*p;
            // (1,0,3,2)
            QuadDofOrd[3][d1] = -1 - ((pm1 - i) + j*p);
            QuadDofOrd[3][d2] = p*pm1 + (pm2 - j) + i*pm1;
            // (2,3,0,1)
            QuadDofOrd[4][d1] = -1 - ((pm1 - i) + (pm2 - j)*p);
            QuadDofOrd[4][d2] = -1 - (p*pm1 + (pm2 - j) + (pm1 - i)*pm1);
            // (2,1,0,3)
            QuadDofOrd[5][d1] = -1 - (p*pm1 + (pm2 - j) + (pm1 - i)*pm1);
            QuadDofOrd[5][d2] = -1 - ((pm1 - i) + (pm2 - j)*p);
            // (3,0,1,2)
            // QuadDofOrd[6][d1] = -1 - (p*pm1 + j + (pm1 - i)*pm1);
            // QuadDofOrd[6][d2] = i + (pm2 - j)*p;
            QuadDofOrd[6][d1] = p*pm1 + (pm2 - j) + i*pm1;
            QuadDofOrd[6][d2] = -1 - ((pm1 - i) + j*p);
            // (3,2,1,0)
            QuadDofOrd[7][d1] = i + (pm2 - j)*p;
            QuadDofOrd[7][d2] = -1 - (p*pm1 + j + (pm1 - i)*pm1);
         }
      }

      int TriDof = ND_dof[Geometry::TRIANGLE];
      TriDofOrd[0] = new int[6*TriDof];
      for (int i = 1; i < 6; i++)
      {
         TriDofOrd[i] = TriDofOrd[i-1] + TriDof;
      }
      // see Mesh::GetTriOrientation in mesh/mesh.cpp,
      // the constructor of H1_FECollection
      for (int j = 0; j <= pm2; j++)
      {
         for (int i = 0; i + j <= pm2; i++)
         {
            int k1 = p*pm1 - (p - j)*(pm1 - j) + 2*i;
            int k2 = p*pm1 - (p - i)*(pm1 - i) + 2*j;
            // (0,1,2)
            TriDofOrd[0][k1  ] = k1;
            TriDofOrd[0][k1+1] = k1 + 1;
            // (0,2,1)
            TriDofOrd[5][k1  ] = k2 + 1;
            TriDofOrd[5][k1+1] = k2;

            // The other orientations can not be supported with the current
            // interface. The method Mesh::ReorientTetMesh will ensure that
            // only orientations 0 and 5 are generated.
         }
      }
   }

   if (dim >= 3)
   {
      ND_Elements[Geometry::CUBE] = new ND_HexahedronElement(p, cb_type, ob_type);
      ND_dof[Geometry::CUBE] = 3*p*pm1*pm1;

      // TODO: cb_type and ob_type for tets
      ND_Elements[Geometry::TETRAHEDRON] = new ND_TetrahedronElement(p);
      ND_dof[Geometry::TETRAHEDRON] = p*pm1*pm2/2;
   }
}

const int *ND_FECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                   int Or) const
{
   if (GeomType == Geometry::SEGMENT)
   {
      return (Or > 0) ? SegDofOrd[0] : SegDofOrd[1];
   }
   else if (GeomType == Geometry::TRIANGLE)
   {
      if (Or != 0 && Or != 5)
      {
         MFEM_ABORT("triangle face orientation " << Or << " is not supported! "
                    "Use Mesh::ReorientTetMesh to fix it.");
      }
      return TriDofOrd[Or%6];
   }
   else if (GeomType == Geometry::SQUARE)
   {
      return QuadDofOrd[Or%8];
   }
   return NULL;
}

FiniteElementCollection *ND_FECollection::GetTraceCollection() const
{
   int p, dim, cb_type, ob_type;

   p = ND_dof[Geometry::SEGMENT];
   if (nd_name[2] == '_') // ND_
   {
      dim = atoi(nd_name + 3);
      cb_type = BasisType::GaussLobatto;
      ob_type = BasisType::GaussLegendre;
   }
   else // ND@
   {
      dim = atoi(nd_name + 6);
      cb_type = BasisType::GetType(nd_name[3]);
      ob_type = BasisType::GetType(nd_name[4]);
   }
   return new ND_Trace_FECollection(p, dim, cb_type, ob_type);
}

ND_FECollection::~ND_FECollection()
{
   delete [] SegDofOrd[0];
   delete [] TriDofOrd[0];
   delete [] QuadDofOrd[0];
   for (int g = 0; g < Geometry::NumGeom; g++)
   {
      delete ND_Elements[g];
   }
}


ND_Trace_FECollection::ND_Trace_FECollection(const int p, const int dim,
                                             const int cb_type,
                                             const int ob_type)
   : ND_FECollection(p, dim-1, cb_type, ob_type)
{
   if (cb_type == BasisType::GaussLobatto &&
       ob_type == BasisType::GaussLegendre)
   {
      snprintf(nd_name, 32, "ND_Trace_%dD_P%d", dim, p);
   }
   else
   {
      snprintf(nd_name, 32, "ND_Trace@%c%c_%dD_P%d",
               (int)BasisType::GetChar(cb_type),
               (int)BasisType::GetChar(ob_type), dim, p);
   }
}


Local_FECollection::Local_FECollection(const char *fe_name)
{
   snprintf(d_name, 32, "Local_%s", fe_name);

   Local_Element = NULL;

   if (!strcmp(fe_name, "BiCubic2DFiniteElement") ||
       !strcmp(fe_name, "Quad_Q3"))
   {
      GeomType = Geometry::SQUARE;
      Local_Element = new BiCubic2DFiniteElement;
   }
   else if (!strcmp(fe_name, "Nedelec1HexFiniteElement") ||
            !strcmp(fe_name, "Hex_ND1"))
   {
      GeomType = Geometry::CUBE;
      Local_Element = new Nedelec1HexFiniteElement;
   }
   else if (!strncmp(fe_name, "H1_", 3))
   {
      GeomType = Geometry::SQUARE;
      Local_Element = new H1_QuadrilateralElement(atoi(fe_name + 7));
   }
   else if (!strncmp(fe_name, "H1Pos_", 6))
   {
      GeomType = Geometry::SQUARE;
      Local_Element = new H1Pos_QuadrilateralElement(atoi(fe_name + 10));
   }
   else if (!strncmp(fe_name, "L2_", 3))
   {
      GeomType = Geometry::SQUARE;
      Local_Element = new L2_QuadrilateralElement(atoi(fe_name + 7));
   }
   else
   {
      mfem::err << "Local_FECollection::Local_FECollection : fe_name = "
                << fe_name << endl;
      mfem_error();
   }
}


NURBSFECollection::NURBSFECollection(int Order)
   : FiniteElementCollection((Order == VariableOrder) ? 1 : Order)
{
   const int order = (Order == VariableOrder) ? 1 : Order;
   SegmentFE        = new NURBS1DFiniteElement(order);
   QuadrilateralFE  = new NURBS2DFiniteElement(order);
   ParallelepipedFE = new NURBS3DFiniteElement(order);

   SetOrder(Order);
}

void NURBSFECollection::SetOrder(int Order) const
{
   mOrder = Order;
   if (Order != VariableOrder)
   {
      snprintf(name, 16, "NURBS%i", Order);
   }
   else
   {
      snprintf(name, 16, "NURBS");
   }
}

NURBSFECollection::~NURBSFECollection()
{
   delete ParallelepipedFE;
   delete QuadrilateralFE;
   delete SegmentFE;
}

const FiniteElement *
NURBSFECollection::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::SEGMENT:     return SegmentFE;
      case Geometry::SQUARE:      return QuadrilateralFE;
      case Geometry::CUBE:        return ParallelepipedFE;
      default:
         mfem_error ("NURBSFECollection: unknown geometry type.");
   }
   return SegmentFE; // Make some compilers happy
}

int NURBSFECollection::DofForGeometry(Geometry::Type GeomType) const
{
   mfem_error("NURBSFECollection::DofForGeometry");
   return 0; // Make some compilers happy
}

const int *NURBSFECollection::DofOrderForOrientation(Geometry::Type GeomType,
                                                     int Or) const
{
   mfem_error("NURBSFECollection::DofOrderForOrientation");
   return NULL;
}

FiniteElementCollection *NURBSFECollection::GetTraceCollection() const
{
   MFEM_ABORT("NURBS finite elements can not be statically condensed!");
   return NULL;
}

}
