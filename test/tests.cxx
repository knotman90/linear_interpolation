//#include <interpolate.hxx>
#include <linearInterpolation.hxx>
#include <interpolate.hxx>
#include <gtest/gtest.h>
#include <ctime>
#include <iostream>

#define NPOLY (1e03)
#define NPOINTS (3e03)
#define UPPER (-1)
#define LOWER (1)
#define EPSILON (1e-8)

using TYPE =  double;

using FL::point;

TYPE d_rand(TYPE lower, TYPE upper)
{
        TYPE f = (TYPE)rand() / RAND_MAX;
        return lower + f * (upper - lower);
}


TEST(powi,power){
        for (int n = 0; n < 10; ++n) {


                for (int i = 0; i < 10; ++i) {
                        ASSERT_EQ(pow(n,i), FL::powi(n,i));

                }
        }
}

TEST(linearInterpolate, linear){
        srand(time(0));
        FL::LinearInterpolator<TYPE> li;
        //1) generate a random rect
        point<1,TYPE> x0;
        point<1,TYPE>x1;
        for(int i=0; i <NPOLY; ++i ) {

                TYPE a = d_rand(LOWER,UPPER);
                TYPE b = d_rand(LOWER,UPPER);


                x0.coords[0] =d_rand(LOWER,UPPER);
                x0.val=a*x0.coord(0)+b;

                x1.coords[0] = x0.coord(0)+d_rand(LOWER,UPPER);
                x1.val=a*x1.coord(0)+b;

                for(int n = 0; n < NPOINTS; ++n) {
                        point<1,TYPE>p;
                        p.coords[0] = d_rand(x0.coord(0),x1.coord(0));
                        //std::cout<<"Testing "<<a<<" "<<b<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        point<1,TYPE> intr_val = li.Linear(p, x0,x1);

                        TYPE exact = a*p.coord(0)+b;
                        ASSERT_NEAR(exact, intr_val.val,EPSILON);

                }
        }
}



TEST(NlinearInterpolate, linear){
        srand(time(0));
        FL::NLinearInterpolator<TYPE> nli;
        //1) generate a random rect
        point<1,TYPE> x0;
        point<1,TYPE>x1;
        for(int i=0; i <NPOLY; ++i ) {

                TYPE a = d_rand(LOWER,UPPER);
                TYPE b = d_rand(LOWER,UPPER);


                x0.coords[0] =d_rand(LOWER,UPPER);
                x0.val=a*x0.coord(0)+b;

                x1.coords[0] = x0.coord(0)+d_rand(LOWER,UPPER);
                x1.val=a*x1.coord(0)+b;

                for(int n = 0; n < NPOINTS; ++n) {
                        point<1,TYPE>p;
                        p.coords[0] = d_rand(x0.coord(0),x1.coord(0));
                        //std::cout<<"Testing "<<a<<" "<<b<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        point<1,TYPE> points[2] ={x0,x1};
                        point<1,TYPE> intr_val = nli.interpolate(p, points);

                        TYPE exact = a*p.coord(0)+b;
                        ASSERT_NEAR(exact, intr_val.val,EPSILON);

                }
        }
}



TYPE getF2Value(const point<2,TYPE>& p, TYPE a0,TYPE a1,TYPE a2,TYPE a3){
        return a0 + a1*p.coord(0) + a2*p.coord(1) +a3*p.coord(0)*p.coord(1);
}



TEST(linearInterpolate, bilinear){
        srand(time(0));
        FL::LinearInterpolator<TYPE> li;


        point<2,TYPE>p0;
        point<2,TYPE>p1;

        point<2,TYPE>p2;
        point<2,TYPE>p3;

        //static test

        TYPE a0 = 1.0;
        TYPE a1 = 1.0;
        TYPE a2 = 1.0;
        TYPE a3 = 1.0;

        p0.coords[0] =0;
        p0.coords[1] = 0;
        p0.val = getF2Value(p0,a0,a1,a2,a3);

        p1.coords[0] =1;
        p1.coords[1] = 0;
        p1.val = getF2Value(p1,a0,a1,a2,a3);

        p2.coords[0] =0;
        p2.coords[1] = 1;
        p2.val = getF2Value(p2,a0,a1,a2,a3);

        p3.coords[0] =1;
        p3.coords[1] = 1;
        p3.val = getF2Value(p3,a0,a1,a2,a3);

        point<2,TYPE>p;
        p.coords[0] = 0.4;
        p.coords[1] = 0.5;
        point<2,TYPE>points[4] = {p0,p1,p2,p3};
        p=li.Bilinear(p, points);
        TYPE exact = getF2Value(p,a0,a1,a2,a3);
        ASSERT_NEAR(exact, p.val,EPSILON);


        //dynamic tests
        for(int i=0; i <NPOLY; ++i ) {

                a0 = d_rand(0.0,1.0);
                a1 = d_rand(0.0,1.0);
                a2 = d_rand(0.0,1.0);
                a3 = d_rand(0.0,1.0);



                p0.coords[0] =d_rand(0,100);
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

                point<2,TYPE>points[4] = {p0,p1,p2,p3};
                for(int n = 0; n < NPOINTS; ++n) {
                        point<2,TYPE>p;
                        p.coords[0] = d_rand(p0.coord(0),p1.coord(0));
                        p.coords[1] = d_rand(p0.coord(1),p2.coord(1));
                        p.val=0.0;
                        // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        li.Bilinear(p, points);

                        exact = getF2Value(p,a0,a1,a2,a3);
                        ASSERT_NEAR(exact, p.val,EPSILON);

                }
        }
}


TEST(NlinearInterpolate, bilinear){
        srand(time(0));

        FL::NLinearInterpolator<TYPE> nli;
        //1) generate a random rect
        point<2,TYPE>p0;
        point<2,TYPE>p1;

        point<2,TYPE>p2;
        point<2,TYPE>p3;

        //static test

        TYPE a0 = 1.0;
        TYPE a1 = 1.0;
        TYPE a2 = 1.0;
        TYPE a3 = 1.0;

        p0.coords[0] =0;
        p0.coords[1] = 0;
        p0.val = getF2Value(p0,a0,a1,a2,a3);

        p1.coords[0] =1;
        p1.coords[1] = 0;
        p1.val = getF2Value(p1,a0,a1,a2,a3);

        p2.coords[0] =0;
        p2.coords[1] = 1;
        p2.val = getF2Value(p2,a0,a1,a2,a3);

        p3.coords[0] =1;
        p3.coords[1] = 1;
        p3.val = getF2Value(p3,a0,a1,a2,a3);

        point<2,TYPE>p;
        p.coords[0] = 0.4;
        p.coords[1] = 0.5;
        point<2,TYPE>points[4] = {p0,p1,p2,p3};
        p=nli.interpolate(p, points);
        TYPE exact = getF2Value(p,a0,a1,a2,a3);
        ASSERT_NEAR(exact, p.val,EPSILON);




        for(int i=0; i <NPOLY; ++i ) {

                a0 = d_rand(0.0,1.0);
                a1 = d_rand(0.0,1.0);
                a2 = d_rand(0.0,1.0);
                a3 = d_rand(0.0,1.0);



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

                point<2,TYPE>points[4] = {p0,p1,p2,p3};
                for(int n = 0; n < NPOINTS; ++n) {
                        point<2,TYPE>p;
                        p.coords[0] = d_rand(p0.coord(0),p1.coord(0));
                        p.coords[1] = d_rand(p0.coord(1),p2.coord(1));
                        // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;

                        p=nli.interpolate<2>(p, &points[0]);

                        TYPE exact = getF2Value(p,a0,a1,a2,a3);
                        ASSERT_NEAR(exact, p.val,EPSILON);

                }
        }
}



TYPE getF3Value(const point<3,TYPE>& p, TYPE a[8]){
        return a[0]+ a[1]*p.coord(0) + a[2]*p.coord(1) +a[3]*p.coord(2) + a[4]*p.coord(0)*p.coord(1)
               +a[5]*p.coord(0)*p.coord(2) + +a[6]*p.coord(1)*p.coord(2)+a[7]*p.coord(0)*p.coord(1)*p.coord(2);
}


TEST(linearInterpolate, trilinear){
        srand(time(0));
        FL::LinearInterpolator<TYPE> li;
        //1) generate a random rect
        point<3,TYPE>points[8];
        TYPE s_x=d_rand(UPPER,LOWER);
        TYPE s_y=d_rand(UPPER,LOWER);
        TYPE s_z=d_rand(UPPER,LOWER);
        points[0].coords[0] = d_rand(UPPER,LOWER);
        points[0].coords[1] = d_rand(UPPER,LOWER);
        points[0].coords[2] = d_rand(UPPER,LOWER);


        //poly function
        TYPE poly[8];

        for(int i=0; i <NPOLY; ++i ) {

                //initialize a random 2D poly
                for (int i = 0; i < 8; ++i) {
                        poly[i]=  d_rand(0.0,1.0);
                }

                //initialize cuboid
                int count=1;
                points[0].val = getF3Value(points[0], poly );
                //points[0].print();
                for (int z = 0; z < 2; ++z) {
                        for (int y = 0; y < 2; ++y) {
                                for (int x = (count==1) ? 1 : 0; x < 2; ++x,count++) {
                                        points[count].coords[0] = points[0].coords[0]+ x * s_x;
                                        points[count].coords[1] = points[0].coords[1]+ y * s_y;
                                        points[count].coords[2] = points[0].coords[2]+ z * s_z;
                                        points[count].val = getF3Value( points[count], poly );
                                        //points[count].print();
                                }

                        }
                }

                //pick up a random point in the interior of the cuboid
                for(int n = 0; n < NPOINTS; ++n) {
                        point<3,TYPE>p;
                        p.coords[0] = d_rand(points[0].coord(0),points[0].coord(0) + s_x);
                        p.coords[1] = d_rand(points[0].coord(1),points[0].coord(1) + s_y);
                        p.coords[2] = d_rand(points[0].coord(2),points[0].coord(2) + s_z);
                        // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        p=li.Trilinear(p, points);

                        TYPE exact = getF3Value(p,poly);
                        ASSERT_NEAR(exact, p.val,EPSILON);

                }
        }
}



TEST(NlinearInterpolate, trilinear){
        srand(time(0));
        FL::NLinearInterpolator<TYPE> nli;
        //1) generate a random rect
        point<3,TYPE>points[8];
        TYPE s_x=d_rand(UPPER,LOWER);
        TYPE s_y=d_rand(UPPER,LOWER);
        TYPE s_z=d_rand(UPPER,LOWER);
        points[0].coords[0] = d_rand(UPPER,LOWER);
        points[0].coords[1] = d_rand(UPPER,LOWER);
        points[0].coords[2] = d_rand(UPPER,LOWER);


        //poly function
        TYPE poly[8];

        for(int i=0; i <NPOLY; ++i ) {

                //initialize a random 2D poly
                for (int i = 0; i < 8; ++i) {
                        poly[i]=  d_rand(0.0,1.0);
                }

                //initialize cuboid
                int count=1;
                points[0].val = getF3Value(points[0], poly );
                //points[0].print();
                for (int z = 0; z < 2; ++z) {
                        for (int y = 0; y < 2; ++y) {
                                for (int x = (count==1) ? 1 : 0; x < 2; ++x,count++) {
                                        points[count].coords[0] = points[0].coords[0]+ x * s_x;
                                        points[count].coords[1] = points[0].coords[1]+ y * s_y;
                                        points[count].coords[2] = points[0].coords[2]+ z * s_z;
                                        points[count].val = getF3Value( points[count], poly );
                                        //points[count].print();
                                }

                        }
                }

                //pick up a random point in the interior of the cuboid
                for(int n = 0; n < NPOINTS; ++n) {
                        point<3,TYPE>p;
                        p.coords[0] = d_rand(points[0].coord(0),points[0].coord(0) + s_x);
                        p.coords[1] = d_rand(points[0].coord(1),points[0].coord(1) + s_y);
                        p.coords[2] = d_rand(points[0].coord(2),points[0].coord(2) + s_z);
                        // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        p=nli.interpolate(p, points);

                        TYPE exact = getF3Value(p,poly);
                        ASSERT_NEAR(exact, p.val,EPSILON);

                }
        }
}

TYPE getF5Value(const point<5,TYPE>& p, TYPE a[32]){
        return a[0]+ a[1]*p.coord(0) + a[2]*p.coord(1) +a[3]*p.coord(2) + a[4]*p.coord(3)+a[5]*p.coord(4)
               +a[5]*p.coord(0)*p.coord(1)  +a[6]*p.coord(0)*p.coord(2)+a[7]*p.coord(0)*p.coord(3) + a[8]*p.coord(0)*p.coord(4)
               ;
}

TEST(NlinearInterpolate, pentalinear){
        srand(time(0));
        FL::NLinearInterpolator<TYPE> nli;
        //1) generate a random rect
        point<5,TYPE>points[32];
        TYPE s_x=d_rand(UPPER,LOWER);
        TYPE s_y=d_rand(UPPER,LOWER);
        TYPE s_z=d_rand(UPPER,LOWER);
        TYPE s_w=d_rand(UPPER,LOWER);
        TYPE s_r=d_rand(UPPER,LOWER);
        points[0].coords[0] = d_rand(UPPER,LOWER);
        points[0].coords[1] = d_rand(UPPER,LOWER);
        points[0].coords[2] = d_rand(UPPER,LOWER);
        points[0].coords[3] = d_rand(UPPER,LOWER);
        points[0].coords[4] = d_rand(UPPER,LOWER);


        //poly function
        TYPE poly[8];

        for(int i=0; i <NPOLY; ++i ) {

                //initialize a random 2D poly
                for (int i = 0; i < 32; ++i) {
                        poly[i]=  d_rand(0.0,1.0);
                }

                //initialize cuboid
                int count=1;
                points[0].val = getF5Value(points[0], poly );
                //points[0].print();
                for (int r = 0; r < 2; ++r)
                        for (int w = 0; w < 2; ++w)
                                for (int z = 0; z < 2; ++z) {
                                        for (int y = 0; y < 2; ++y) {
                                                for (int x = (count==1) ? 1 : 0; x < 2; ++x,count++) {
                                                        points[count].coords[0] = points[0].coords[0]+ x * s_x;
                                                        points[count].coords[1] = points[0].coords[1]+ y * s_y;
                                                        points[count].coords[2] = points[0].coords[2]+ z * s_z;
                                                        points[count].coords[3] = points[0].coords[3]+ w * s_w;
                                                        points[count].coords[4] = points[0].coords[4]+ r * s_r;
                                                        points[count].val = getF5Value( points[count], poly );
                                                        //points[count].print();
                                                }

                                        }
                                }

                //pick up a random point in the interior of the cuboid
                for(int n = 0; n < NPOINTS; ++n) {
                        point<5,TYPE>p;
                        p.coords[0] = d_rand(points[0].coord(0),points[0].coord(0) + s_x);
                        p.coords[1] = d_rand(points[0].coord(1),points[0].coord(1) + s_y);
                        p.coords[2] = d_rand(points[0].coord(2),points[0].coord(2) + s_z);
                        p.coords[3] = d_rand(points[0].coord(3),points[0].coord(3) + s_w);
                        p.coords[4] = d_rand(points[0].coord(4),points[0].coord(4) + s_r);

                        // std::cout<<"Testing "<<a0<<" "<<a1<<" "<<x0.coord(0)<< " "<<x1.coord(0)<<" "<<p.coord(0)<<std::endl;
                        p=nli.interpolate(p, points);

                        TYPE exact = getF5Value(p,poly);
                        ASSERT_NEAR(exact, p.val,EPSILON);

                }
        }
}





int main(int argc, char** argv){



        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();

        return 0;
}
