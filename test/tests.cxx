#include <linearInterpolation.hxx>
#include <interpolate.hxx>
#include <gtest/gtest.h>
#include <ctime>
#include <iostream>

#define NPOLY (1e03)
#define NPOINTS (1e03)
#define UPPER (0)
#define LOWER (1)

#define COEFFLOWER (0)
#define COEFFUPPER (1)

#define EPSILON (1e-8)

using TYPE =  long double;

using FL::point;

TYPE d_rand(TYPE lower, TYPE upper) {
  TYPE f = (TYPE)rand() / RAND_MAX;

  return lower + f * (upper - lower);
}

template<int DIM>
TYPE getFNValue(const point<DIM, TYPE>& p, TYPE *a) {
  int  n     = DIM;
  int  count = 1;
  TYPE tot   = a[0];

  for (int r = 1; r <= DIM; r++) {
    std::vector<bool> v(n);
    std::fill(v.begin(), v.end() - n + r, true);

    do {
      TYPE val = a[count];

      for (int i = 0; i < n; ++i)
        if (v[i]) val *= p.coord(i);

      count++;
      tot += val;
    } while (std::prev_permutation(v.begin(), v.end()));
  }
  return tot;
}

// ----------ah-hoc (lin/bi/tri-linear) interpolator class TESTS-------
TEST(linearInterpolate, linear) {
  srand(time(0));
  FL::LinearInterpolator<TYPE> li;

  // generate a random rect
  point<1, TYPE> x0;
  point<1, TYPE> x1;

  for (int i = 0; i < NPOLY; ++i) {
    TYPE a = d_rand(COEFFLOWER, COEFFUPPER);
    TYPE b = d_rand(COEFFLOWER, COEFFUPPER);


    x0.coords[0] = d_rand(LOWER, UPPER);
    x0.val       = a * x0.coord(0) + b;

    x1.coords[0] = x0.coord(0) + d_rand(LOWER, UPPER);
    x1.val       = a * x1.coord(0) + b;

    for (int n = 0; n < NPOINTS; ++n) {
      point<1, TYPE> p;
      p.coords[0] = d_rand(x0.coord(0), x1.coord(0));

      point<1, TYPE> intr_val = li.Linear(p, x0, x1);

      TYPE exact = a * p.coord(0) + b;
      ASSERT_NEAR(exact, intr_val.val, EPSILON);
    }
  }
}


TEST(linearInterpolate, bilinear) {
  srand(time(0));
  FL::LinearInterpolator<TYPE> li;


  point<2, TYPE> p0;
  point<2, TYPE> p1;

  point<2, TYPE> p2;
  point<2, TYPE> p3;


  // static test

  TYPE a0 = 1.0;
  TYPE a1 = 1.0;
  TYPE a2 = 1.0;
  TYPE a3 = 1.0;

  TYPE a[4] = { a0, a1, a2, a3 };

  p0.coords[0] = 0;
  p0.coords[1] = 0;
  p0.val       = getFNValue<2>(p0, a);

  p1.coords[0] = 1;
  p1.coords[1] = 0;
  p1.val       = getFNValue<2>(p1, a);

  p2.coords[0] = 0;
  p2.coords[1] = 1;
  p2.val       = getFNValue<2>(p2, a);

  p3.coords[0] = 1;
  p3.coords[1] = 1;
  p3.val       = getFNValue<2>(p3, a);

  point<2, TYPE> p;
  p.coords[0] = 0.4;
  p.coords[1] = 0.5;
  point<2, TYPE> points[4] = { p0, p1, p2, p3 };
  p = li.Bilinear(p, points);
  TYPE exact = getFNValue<2>(p, a);
      ASSERT_NEAR(exact, p.val, EPSILON);


  // dynamic tests
  for (int i = 0; i < NPOLY; ++i) {
    a0 = d_rand(COEFFLOWER, COEFFUPPER);
    a1 = d_rand(COEFFLOWER, COEFFUPPER);
    a2 = d_rand(COEFFLOWER, COEFFUPPER);
    a3 = d_rand(COEFFLOWER, COEFFUPPER);
    TYPE a[4] = { a0, a1, a2, a3 };


    p0.coords[0] = d_rand(LOWER, UPPER);
    p0.coords[1] = d_rand(LOWER, UPPER);
    p0.val       = getFNValue<2>(p0, a);

    p2.coords[0] = p0.coords[0];
    p2.coords[1] = d_rand(LOWER, UPPER);
    p2.val       = getFNValue<2>(p2, a);

    p1.coords[0] = p0.coord(0) + d_rand(0.0, 100);
    p1.coords[1] =  p0.coords[1];
    p1.val       = getFNValue<2>(p1, a);

    p3.coords[0] = p1.coords[0];
    p3.coords[1] =  p2.coords[1];
    p3.val       = getFNValue<2>(p3, a);

    point<2, TYPE> points[4] = { p0, p1, p2, p3 };

    for (int n = 0; n < NPOINTS; ++n) {
      point<2, TYPE> p;
      p.coords[0] = d_rand(p0.coord(0), p1.coord(0));
      p.coords[1] = d_rand(p0.coord(1), p2.coord(1));
      p.val       = 0.0;

      li.Bilinear(p, points);

      exact = getFNValue<2>(p, a);
      ASSERT_NEAR(exact, p.val, EPSILON);
    }
  }
}


TEST(linearInterpolate, trilinear) {
  srand(time(0));
  FL::LinearInterpolator<TYPE> li;

  // 1) generate a random rect
  point<3, TYPE> points[8];
  TYPE s_x = d_rand(UPPER, LOWER);
  TYPE s_y = d_rand(UPPER, LOWER);
  TYPE s_z = d_rand(UPPER, LOWER);
  points[0].coords[0] = d_rand(UPPER, LOWER);
  points[0].coords[1] = d_rand(UPPER, LOWER);
  points[0].coords[2] = d_rand(UPPER, LOWER);


  // poly function
  TYPE poly[8];

  for (int i = 0; i < NPOLY; ++i) {
    // initialize a random 2D poly
    for (int i = 0; i < 8; ++i) {
      poly[i] =   d_rand(LOWER, UPPER);
    }

    // initialize cuboid
    int count = 1;
    points[0].val = getFNValue<3>(points[0], poly);

    for (int z = 0; z < 2; ++z)
      for (int y = 0; y < 2; ++y)
        for (int x = (count == 1) ? 1 : 0; x < 2; ++x, count++) {
          points[count].coords[0] = points[0].coords[0] + x * s_x;
          points[count].coords[1] = points[0].coords[1] + y * s_y;
          points[count].coords[2] = points[0].coords[2] + z * s_z;
          points[count].val       = getFNValue<3>(points[count], poly);
        }

    // pick up a random point in the interior of the cuboid
    for (int n = 0; n < NPOINTS; ++n) {
      point<3, TYPE> p;
      p.coords[0] = d_rand(points[0].coord(0), points[0].coord(0) + s_x);
      p.coords[1] = d_rand(points[0].coord(1), points[0].coord(1) + s_y);
      p.coords[2] = d_rand(points[0].coord(2), points[0].coord(2) + s_z);

      p = li.Trilinear(p, points);

      TYPE exact = getFNValue<3>(p, poly);
      ASSERT_NEAR(exact, p.val, EPSILON);
    }
  }
}


// --------N DIMENSIONAL LINEAR INTERPOLATOR CLASS--------

TEST(NlinearInterpolate, linear) {
  srand(time(0));
  FL::NLinearInterpolator<TYPE> nli;

  // 1) generate a random rect
  point<1, TYPE> x0;
  point<1, TYPE> x1;

  for (int i = 0; i < NPOLY; ++i) {
    TYPE a = d_rand(COEFFLOWER, COEFFUPPER);
    TYPE b = d_rand(COEFFLOWER, COEFFUPPER);


    x0.coords[0] = d_rand(LOWER, UPPER);
    x0.val       = a * x0.coord(0) + b;

    x1.coords[0] = x0.coord(0) + d_rand(LOWER, UPPER);
    x1.val       = a * x1.coord(0) + b;

    for (int n = 0; n < NPOINTS; ++n) {
      point<1, TYPE> p;
      p.coords[0] = d_rand(x0.coord(0), x1.coord(0));

      point<1, TYPE> points[2] = { x0, x1 };
      point<1, TYPE> intr_val  = nli.interpolate(p, points);

      TYPE exact = a * p.coord(0) + b;
      ASSERT_NEAR(exact, intr_val.val, EPSILON);
    }
  }
}
TEST(NlinearInterpolate, bilinear) {
  srand(time(0));

  FL::NLinearInterpolator<TYPE> nli;

  // 1) generate a random rect
  point<2, TYPE> p0;
  point<2, TYPE> p1;

  point<2, TYPE> p2;
  point<2, TYPE> p3;

  // static test

  TYPE a0   = 1.0;
  TYPE a1   = 1.0;
  TYPE a2   = 1.0;
  TYPE a3   = 1.0;
  TYPE a[4] = { a0, a1, a2, a3 };
  p0.coords[0] = 0;
  p0.coords[1] = 0;
  p0.val       = getFNValue<2>(p0, a);

  p1.coords[0] = 1;
  p1.coords[1] = 0;
  p1.val       = getFNValue<2>(p1, a);

  p2.coords[0] = 0;
  p2.coords[1] = 1;
  p2.val       = getFNValue<2>(p2, a);

  p3.coords[0] = 1;
  p3.coords[1] = 1;
  p3.val       = getFNValue<2>(p3, a);

  point<2, TYPE> p;
  p.coords[0] = 0.4;
  p.coords[1] = 0.5;
  point<2, TYPE> points[4] = { p0, p1, p2, p3 };
  p = nli.interpolate(p, points);
  TYPE exact = getFNValue<2>(p, a);
      ASSERT_NEAR(exact, p.val, EPSILON);


  for (int i = 0; i < NPOLY; ++i) {
    a0 = d_rand(COEFFLOWER, COEFFUPPER);
    a1 = d_rand(COEFFLOWER, COEFFUPPER);
    a2 = d_rand(COEFFLOWER, COEFFUPPER);
    a3 = d_rand(COEFFLOWER, COEFFUPPER);
    TYPE a[4] = { a0, a1, a2, a3 };


    p0.coords[0] = d_rand(0, 100);
    p0.coords[1] = d_rand(LOWER, UPPER);
    p0.val       = getFNValue<2>(p0, a);

    p2.coords[0] = p0.coords[0];
    p2.coords[1] = d_rand(LOWER, UPPER);
    p2.val       = getFNValue<2>(p2, a);

    p1.coords[0] = p0.coord(0) + d_rand(0.0, 100);
    p1.coords[1] =  p0.coords[1];
    p1.val       = getFNValue<2>(p1, a);

    p3.coords[0] = p1.coords[0];
    p3.coords[1] =  p2.coords[1];
    p3.val       = getFNValue<2>(p3, a);

    point<2, TYPE> points[4] = { p0, p1, p2, p3 };

    for (int n = 0; n < NPOINTS; ++n) {
      point<2, TYPE> p;
      p.coords[0] = d_rand(p0.coord(0), p1.coord(0));
      p.coords[1] = d_rand(p0.coord(1), p2.coord(1));


      p = nli.interpolate<2>(p, &points[0]);

      TYPE exact = getFNValue<2>(p, a);
      ASSERT_NEAR(exact, p.val, EPSILON);
    }
  }
}


TEST(NlinearInterpolate, trilinear) {
  srand(time(0));
  FL::NLinearInterpolator<TYPE> nli;

  // 1) generate a random rect
  point<3, TYPE> points[8];
  TYPE s_x = d_rand(UPPER, LOWER);
  TYPE s_y = d_rand(UPPER, LOWER);
  TYPE s_z = d_rand(UPPER, LOWER);
  points[0].coords[0] = d_rand(UPPER, LOWER);
  points[0].coords[1] = d_rand(UPPER, LOWER);
  points[0].coords[2] = d_rand(UPPER, LOWER);

  // poly function
  TYPE poly[8];

  for (int i = 0; i < NPOLY; ++i) {
    // initialize a random 2D poly
    for (int i = 0; i < 8; ++i) {
      poly[i] =   d_rand(LOWER, UPPER);
    }

    // initialize cuboid
    int count = 1;
    points[0].val = getFNValue<3>(points[0], poly);

    for (int z = 0; z < 2; ++z)
      for (int y = 0; y < 2; ++y)
        for (int x = (count == 1) ? 1 : 0; x < 2; ++x, count++) {
          points[count].coords[0] = points[0].coords[0] + x * s_x;
          points[count].coords[1] = points[0].coords[1] + y * s_y;
          points[count].coords[2] = points[0].coords[2] + z * s_z;
          points[count].val       = getFNValue<3>(points[count], poly);
        }

    // pick up a random point in the interior of the cuboid
    for (int n = 0; n < NPOINTS; ++n) {
      point<3, TYPE> p;
      p.coords[0] = d_rand(points[0].coord(0), points[0].coord(0) + s_x);
      p.coords[1] = d_rand(points[0].coord(1), points[0].coord(1) + s_y);
      p.coords[2] = d_rand(points[0].coord(2), points[0].coord(2) + s_z);

      p = nli.interpolate(p, points);

      TYPE exact = getFNValue<3>(p, poly);
      ASSERT_NEAR(exact, p.val, EPSILON);
    }
  }
}

TEST(NlinearInterpolate, fourlinear) {
  srand(time(0));
  FL::NLinearInterpolator<TYPE> nli;

  // 1) generate a random rect
  point<4, TYPE> points[16];
  TYPE s_x = d_rand(UPPER, LOWER);
  TYPE s_y = d_rand(UPPER, LOWER);
  TYPE s_z = d_rand(UPPER, LOWER);
  TYPE s_w = d_rand(UPPER, LOWER);

  points[0].coords[0] = d_rand(UPPER, LOWER);
  points[0].coords[1] = d_rand(UPPER, LOWER);
  points[0].coords[2] = d_rand(UPPER, LOWER);
  points[0].coords[3] = d_rand(UPPER, LOWER);

  // poly function
  TYPE poly[16];

  for (int i = 0; i < NPOLY; ++i) {
    // initialize a random 2D poly
    for (int i = 0; i < 16; ++i) poly[i] =   d_rand(LOWER, UPPER);

    // initialize cuboid
    int count = 1;
    points[0].val = getFNValue<4>(points[0], poly);

    for (int w = 0; w < 2; ++w)
      for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
          for (int x = (count == 1) ? 1 : 0; x < 2; ++x, count++) {
            points[count].coords[0] = points[0].coords[0] + x * s_x;
            points[count].coords[1] = points[0].coords[1] + y * s_y;
            points[count].coords[2] = points[0].coords[2] + z * s_z;
            points[count].coords[3] = points[0].coords[3] + w * s_w;
            points[count].val       = getFNValue<4>(points[count], poly);
          }

    // pick up a random point in the interior of the cuboid
    for (int n = 0; n < NPOINTS; ++n) {
      point<4, TYPE> p;
      p.coords[0] = d_rand(points[0].coord(0), points[0].coord(0) + s_x);
      p.coords[1] = d_rand(points[0].coord(1), points[0].coord(1) + s_y);
      p.coords[2] = d_rand(points[0].coord(2), points[0].coord(2) + s_z);
      p.coords[3] = d_rand(points[0].coord(3), points[0].coord(3) + s_w);

      p = nli.interpolate(p, points);

      TYPE exact = getFNValue<4>(p, poly);
      ASSERT_NEAR(exact, p.val, EPSILON);
    }
  }
}


int main(int argc, char **argv) {
  int a = FL::sumPow2<3, int>::value;

  std::cout << "a is " << a << std::endl;

  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();

  return 0;
}
