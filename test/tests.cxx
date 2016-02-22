//#include <interpolate.hxx>
#include <linearInterpolation.hxx>
#include <gtest/gtest.h>
#include <ctime>
#include<iostream>

#define NPOLY (1e04)
#define NPOINTS (2e03)
#define UPPER (-1)
#define LOWER (1)
#define EPSILON (1e-8)

using TYPE =  long double;

using FL::point;
//add abs
TYPE d_rand(TYPE lower, TYPE upper)
{
    TYPE f = (TYPE)rand() / RAND_MAX;
    return lower + f * (upper - lower);
}


TEST(linearInterpolate, linear){
    srand(time(0));
    FL::LinearInterpolator<TYPE> li;
    //1) generate a random rect
    point<1> x0;
    point<1> x1;
    for(int i=0;i <NPOLY ;++i ){

        TYPE a = d_rand(LOWER,UPPER);
        TYPE b = d_rand(LOWER,UPPER);


        x0.coords[0] =d_rand(LOWER,UPPER);
        x0.val=a*x0.coord(0)+b;

        x1.coords[0] = x0.coord(0)+d_rand(LOWER,UPPER);
        x1.val=a*x1.coord(0)+b;

        for(int n = 0; n < NPOINTS ; ++n){
            point<1> p;
            p.coords[0] = d_rand(x0.coord(0),x1.coord(0));
            //std::cout<<"Testing "<<a<<" "<<b<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
            point<1>  intr_val = li.Linear(p, x0,x1);

            TYPE exact = a*p.coord(0)+b;
            ASSERT_NEAR(exact, intr_val.val,EPSILON);

        }
    }
}


TYPE getF2Value(const point<2>& p, TYPE a0,TYPE a1,TYPE a2,TYPE a3){
    return a0 + a1*p.coord(0) + a2*p.coord(1) +a3*p.coord(0)*p.coord(1);
}



TEST(linearInterpolate, bilinear){
    srand(time(0));
    FL::LinearInterpolator<TYPE> li;
    //1) generate a random rect
    point<2> p0;
    point<2> p1;

    point<2> p2;
    point<2> p3;
    for(int i=0;i <NPOLY ;++i ){

        TYPE a0 = d_rand(0.0,1.0);
        TYPE a1 = d_rand(0.0,1.0);
        TYPE a2 = d_rand(0.0,1.0);
        TYPE a3 = d_rand(0.0,1.0);



        p0.coords[0] =d_rand(00,100);
        p0.coords[1] = d_rand(0,100);
        p0.val = getF2Value(p0,a0,a1,a2,a3);

        p2.coords[0] = p0.coords[0];
        p2.coords[1] = d_rand(0,100);
        p2.val= getF2Value(p2,a0,a1,a2,a3);

        p1.coords[0] = p0.coord(0)+d_rand(0.0,100);
        p1.coords[1] =  p0.coords[1];
        p1.val= getF2Value(p1,a0,a1,a2,a3);

        p3.coords[0] = p1.coords[0];
        p3.coords[1] =  p2.coords[1];
        p3.val= getF2Value(p3,a0,a1,a2,a3);

        point<2> points[4] = {p0,p1,p2,p3};
        for(int n = 0; n < NPOINTS ; ++n){
            point<2> p;
            p.coords[0] = d_rand(p0.coord(0),p1.coord(0));
            p.coords[1] = d_rand(p0.coord(1),p2.coord(1));
            // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
            li.Bilinear(p, points);

            TYPE exact = getF2Value(p,a0,a1,a2,a3);
            ASSERT_NEAR(exact, p.val,EPSILON);

        }
    }
}



TYPE getF3Value(const point<3>& p, TYPE a[8]){
    return a[0]+ a[1]*p.coord(0) + a[2]*p.coord(1) +a[3]*p.coord(2) + a[4]*p.coord(0)*p.coord(1)
            +a[5]*p.coord(0)*p.coord(2) + +a[6]*p.coord(1)*p.coord(2)+a[7]*p.coord(0)*p.coord(1)*p.coord(2);
}


TEST(linearInterpolate, trilinear){
    srand(time(0));
    FL::LinearInterpolator<TYPE> li;
    //1) generate a random rect
    point<3> points[8];
    TYPE s_x=d_rand(UPPER,LOWER);
    TYPE s_y=d_rand(UPPER,LOWER);
    TYPE s_z=d_rand(UPPER,LOWER);
    points[0].coords[0] = d_rand(UPPER,LOWER);
    points[0].coords[1] = d_rand(UPPER,LOWER);
    points[0].coords[2] = d_rand(UPPER,LOWER);


    //poly function
    TYPE poly[8];

    for(int i=0;i <NPOLY ;++i ){

        //initialize a random 2D poly
        for (int i = 0; i < 8; ++i) {
            poly[i]=  d_rand(0.0,1.0);
        }

       //initialize cuboid
        int count=1;
        points[0].val = getF3Value(points[0] , poly );
        //points[0].print();
        for (int z = 0 ; z < 2; ++z) {
            for (int y = 0; y < 2; ++y) {
                for (int x = (count==1) ? 1 : 0; x < 2; ++x,count++) {
                    points[count].coords[0] = points[0].coords[0]+ x * s_x;
                    points[count].coords[1] = points[0].coords[1]+ y * s_y;
                    points[count].coords[2] = points[0].coords[2]+ z * s_z;
                    points[count].val = getF3Value( points[count] , poly );
                    //points[count].print();
                }

            }
        }

        //pick up a random point in the interior of the cuboid
        for(int n = 0; n < NPOINTS ; ++n){
            point<3> p;
            p.coords[0] = d_rand(points[0].coord(0),points[0].coord(0) + s_x);
            p.coords[0] = d_rand(points[0].coord(1),points[0].coord(1) + s_y);
            p.coords[1] = d_rand(points[0].coord(2),points[0].coord(2) + s_z);
            // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
            p=li.Trilinear(p, points);

            TYPE exact = getF3Value(p,poly);
            ASSERT_NEAR(exact, p.val,EPSILON);

        }
    }
}



int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

    return 0;
}
