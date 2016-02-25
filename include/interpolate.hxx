#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include<assert.h>
#include<linearInterpolation.hxx>
namespace FL{

template<int exp, typename D = int>
struct sumPow2{
    enum { value = ((D)(1<<exp)) + sumPow2<(exp-1),D>::value };
};

template<>
struct sumPow2<0,int>{
    enum { value = 1 };
};

/**
 *Performs Linear interpolation within a N dimensional cuboid
 */
template <class T>
class NLinearInterpolator{

public:

    template<int DIM>
    point<DIM,T>& interpolate(point<DIM,T>& p, point<DIM,T>* v){

        constexpr int sumlog=sumPow2<DIM>::value;
        T cs[sumlog];
#pragma vector always
        for (int i = 0; i < (1<<DIM); ++i)
            cs[i]=v[i].val;

        T ds[DIM];

        int NP=1;

        for (int d = 0; d < DIM; ++d) {
            ds[d] = (p.coord(d)-v[0].coord(d)) / ((v[NP].coord(d) - v[0].coord(d)));
            NP*=2;
        }

        int count=NP;
        NP/=2;
        int d =0;
        //O(log2) iteration
        while(NP>1){

            for (int i = 0; i < NP; ++i){
                const int idx = count-(NP*2) + i*2;
                cs[count+i]= cs[idx]*(1-ds[d]) + cs[idx +1] * ds[d];
            }
            count+=NP;
            NP/=2;
            ++d;
        }

        cs[sumlog-1] = cs[sumlog-3]* (1-ds[DIM-1]) + cs[sumlog-2]* ds[DIM-1];
        p.val=cs[sumlog-1];

        return p;
    }

}; //class interpolator

} //namespace First Light
#endif /* INTERPOLATE */
