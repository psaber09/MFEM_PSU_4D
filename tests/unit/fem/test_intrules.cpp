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

#include "mfem.hpp"
using namespace mfem;

#include "unit_tests.hpp"

//You typically want to start by testing things one object at a time.
TEST_CASE("Integration rule container with no refinement", "[IntegrationRules]")
{
   //This code is automatically re-executed for all of the sections.
   IntegrationRules my_intrules(0, Quadrature1D::GaussLegendre);
   IntegrationRule single_point_rule(1);
   const IntegrationRule *ir;

   //The tests will be reported in these sections.
   //Each REQUIRE counts as an assertion.
   //true = pass, false = fail
   SECTION("can set int rules in empty container")
   {
      single_point_rule[0].weight = 1.0;
      my_intrules.SetOwnRules(0);

      my_intrules.Set(Geometry::SEGMENT, 0, single_point_rule);
      ir = &my_intrules.Get(Geometry::SEGMENT, 0);

      REQUIRE(ir->Size() == 1);
      REQUIRE(ir->IntPoint(0).weight == 1.0);
   }

   SECTION("user set int rules really are owned by the user")
   {
      //Set up a custom int rule and put it in the container
      single_point_rule[0].weight = 1.0;
      my_intrules.SetOwnRules(0);
      my_intrules.Set(Geometry::SEGMENT, 0, single_point_rule);

      //Alter the int rule
      single_point_rule[0].weight = 2.0;

      //Test ownership by making sure that the int rule has changed
      ir = &my_intrules.Get(Geometry::SEGMENT, 0);
      REQUIRE(ir->IntPoint(0).weight == 2.0);
   }

   SECTION("point int rules 0, 1 accessible")
   {
      ir = &my_intrules.Get(Geometry::POINT, 0);
      REQUIRE(ir->Size() == 1);

      ir = &my_intrules.Get(Geometry::POINT, 1);
      REQUIRE(ir->Size() == 1);
   }

   //Can't really unit test these because it will crash due to
   //a null pointer dereference if it doesn't work.  This will
   //force the crash in the unit tests though.
   SECTION("resize of the SEGMENT int rule array")
   {
      ir = &my_intrules.Get(Geometry::SEGMENT, 100);
      REQUIRE(true);
   }

   SECTION("setting the integration point index works")
   {
      ir = &my_intrules.Get(Geometry::CUBE, 5);
      for (int i = 0; i < ir->Size(); i++)
      {
         REQUIRE(ir->IntPoint(i).index == i);
      }
   }

   SECTION("intrules up to order 16 accessible")
   {
      for (int order = 0; order <= 16; order ++)
      {
         //Do this in reverse the usual order to make sure that
         //the higher dimension cases are causing their constituent
         //lower dimension cases to lazily create properly.
         my_intrules.Get(Geometry::CUBE,        order);
         my_intrules.Get(Geometry::TETRAHEDRON, order);
         my_intrules.Get(Geometry::SQUARE,      order);
         my_intrules.Get(Geometry::TRIANGLE,    order);
         my_intrules.Get(Geometry::SEGMENT,     order);
      }
      REQUIRE(true);
   }
}


double poly2d(const IntegrationPoint &ip, int m, int n)
{
   return pow(ip.x, m)*pow(ip.y, n);
   // exact integral over the reference triangle is
   // m!n!/(m+n+2)! = 1/binom(m+n,m)/(m+n+1)/(m+n+2)
}

double apoly2d(const IntegrationPoint &ip, int i, int j, int k)
{
   return pow(1. - ip.x - ip.y, i)*pow(ip.x, j)*pow(ip.y, k);
   // exact integral over the reference triangle is (with p = i+j+k)
   // i!j!k!/(p+2)! = 1/binom(p,i+j)/binom(i+j,i)/(p+1)/(p+2)
}

double poly3d(const IntegrationPoint &ip, int l, int m, int n)
{
   return pow(ip.x, l)*pow(ip.y, m)*pow(ip.z, n);
   // exact integral over the reference tetrahedron is (with p = l+m+n)
   // l!m!n!/(p+3)! = 1/binom(p,l+m)/binom(l+m,l)/(p+1)/(p+2)/(p+3)
}

double poly4d(const IntegrationPoint &ip, int l, int m, int n, int o)
{
   return pow(ip.x, l)*pow(ip.y, m)*pow(ip.z, n)*pow(ip.t, o);
   // exact integral over the reference pentatope is (with p = l+m+n+o)
   // l!m!n!o!/(p+4)! = 1/binom(p,l)/binom(p-l,m)/binom(n+o,o)/(p+1)/(p+2)/(p+3)/(p+4)
}

TEST_CASE("Simplex integration rules", "[SimplexRules]")
{
   //This code is automatically re-executed for all of the sections.
   IntegrationRules my_intrules(0, Quadrature1D::GaussLegendre);
   IntegrationRule single_point_rule(1);

   const int maxn = 32;
   int binom[maxn+1][maxn+1];
   for (int n = 0; n <= maxn; n++)
   {
      binom[n][0] = binom[n][n] = 1;
      for (int k = 1; k < n; k++)
      {
         binom[n][k] = binom[n-1][k] + binom[n-1][k-1];
      }
   }

   SECTION("low triangle integration error on reference element for f=x^m y^n, where m+n <= p")
   {
      for (int order = 0; order <= 25; order++)
      {
         const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, order);

         // using the monomial basis: x^m y^n, 0 <= m+n <= order
         for (int p = 0; p <= order; p++)
         {
            for (int m = p; m >= 0; m--)
            {
               int n = p - m;

               double integral = 0.0;
               for (int i = 0; i < ir.GetNPoints(); i++)
               {
                  const IntegrationPoint &ip = ir.IntPoint(i);
                  integral += ip.weight*poly2d(ip, m, n);
               }

               double exact = 1.0/binom[p][m]/(p + 1)/(p + 2);
               double relerr = 1. - integral/exact;

               //If a test fails any INFO statements preceding the REQUIRE are displayed
               INFO("p=" << p << ", m=" << m << ", n=" << n);
               REQUIRE(fabs(relerr) < 1e-11);
            }
         }
      }
   }

   SECTION("low tet integration error on reference element for f=x^l y^m z^n, where l+m+n <= p")
   {
      for (int order = 0; order <= 21; order++)
      {
         const IntegrationRule &ir = IntRules.Get(Geometry::TETRAHEDRON, order);

         for (int p = 0; p <= order; p++)
         {
            for (int l = p; l >= 0; l--)
            {
               for (int m = p - l; m >= 0; m--)
               {
                  int n = p - l - m;

                  double integral = 0.0;
                  for (int i = 0; i < ir.GetNPoints(); i++)
                  {
                     const IntegrationPoint &ip = ir.IntPoint(i);
                     integral += ip.weight*poly3d(ip, l, m, n);
                  }

                  double exact = 1.0/binom[p][l+m]/binom[l+m][l]/(p+1)/(p+2)/(p+3);
                  double relerr = 1. - integral/exact;

                  //If a test fails any INFO statements preceding the REQUIRE are displayed
                  INFO("p=" << p << ", l=" << l << ", m=" << m << ", n=" << n);
                  REQUIRE(fabs(relerr) < 1e-11);
               }
            }
         }
      }
   }

   SECTION("low pent integration error on reference element for f=x^l y^m z^n t^o, where l+m+n+o <= p")
   {
      for (int order = 0; order <= 10; order++)
      {
         const IntegrationRule &ir = IntRules.Get(Geometry::PENTATOPE, order);

         for (int p = 0; p <= order; p++)
         {
            for (int l = p; l >= 0; l--)
            {
               for (int m = p - l; m >= 0; m--)
               {
                  for (int n = p - l - m; n >= 0; n--)
                  {
                     int o = p - l - m - n;

                     double integral = 0.0;
                     for (int i = 0; i < ir.GetNPoints(); i++)
                     {
                        const IntegrationPoint &ip = ir.IntPoint(i);
                        integral += ip.weight*poly4d(ip, l, m, n, o);
                     }

                     double exact = 1.0/binom[p][l]/binom[p-l][m]/binom[n+o][o]/(p+1)/(p+2)/(p+3)/(p+4);
                     double relerr = 1. - integral/exact;

                     //If a test fails any INFO statements preceding the REQUIRE are displayed
                     INFO("p=" << p << ", l=" << l << ", m=" << m << ", n=" << n << ", o=" << o);
                     REQUIRE(fabs(relerr) < 1e-11);
                  }
               }
            }
         }
      }
   }
}
