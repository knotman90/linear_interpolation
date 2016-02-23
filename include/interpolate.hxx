#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include<assert.h>
#include<linearInterpolation.hxx>
namespace FL{


int powi(int v, int exp){
    if(exp==0)
        return 1;
    int val=v;
    for(int i=2;i<=exp;i++){
        val=val*v;
    }
    return val;
}

template <class T>
class NLinearInterpolator{
    typedef T number_type;
    typedef T* number_pointer;
    typedef T& number_reference;
public:



    template<int DIM>
    point<DIM,T>& interpolate(point<DIM,T>& p, point<DIM,T>* v){

        T ds[DIM];

        int NP=1;
        int sumlog=0;
        for (int d = 0; d < DIM; ++d) {

            ds[d] = (p.coord(d)-v[0].coord(d)) * (1/(v[NP].coord(d) - v[0].coord(d)));
            //printf("dim is %i target is %i ds[%i]=%f\n",d,NP,d,ds[d]);
            sumlog+=NP;
            NP*=2;
        }
        sumlog+=NP;

        //printf("%i,%i\n",NP,sumlog);
        T cs[sumlog];
        for (int i = 0; i < NP; ++i) {
            cs[i]=v[i].val;
            //printf("val %f ",v[i].val);
        }
        //printf("\n");
        int count=NP;
        NP/=2;
        int d =0;
        while(NP>1){
            for (int i = 0; i < NP; ++i) {
                cs[count+i]= cs[count-(NP*2) + i*2]*(1-ds[d]) + cs[count-(NP*2) + i*2 +1] * ds[d];
                //printf("computing %i using %i %i : %f %i\n",(count+i),(count-(NP*2) + i*2),(count-(NP*2) + i*2 +1) ,  cs[count+i],d);

            }
            count+=NP;
            NP/=2;
            d++;
        }


        cs[sumlog-1] = cs[sumlog-3]* (1-ds[DIM-1]) + cs[sumlog-2]* ds[DIM-1];
        //printf("last step= %f %f %f %f\n",cs[sumlog-3] , (1-ds[DIM-1]) ,cs[sumlog-2], ds[DIM-1] ),
        //printf("computing %i using %i %i : %f %i\n",(sumlog-1),(sumlog-3),(sumlog-2),cs[sumlog-1],d);
        p.val=cs[sumlog-1];
        return p;
    }


private:



}; //class interpolator

} //namespace First Light
#endif /* INTERPOLATE */
