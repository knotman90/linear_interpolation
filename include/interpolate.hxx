#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include<assert.h>

namespace FL{

typedef unsigned int uint;

template< int DIM, typename T = double>
class point{

    T coords [DIM] ;
    T val;

    T coord(uint c) {
        assert(c >= 0  && c < DIM);
        return coords[c];
    }
};

template <class T>
class LinearInterpolator{
    typedef T number_type;
    typedef T* number_pointer;
    typedef T& number_reference;
public:
    //(2)^DIM points
    template<int DIM >
    T interpolate(point<DIM>& p, point<DIM>* v){
        T constants[DIM];
        for(int d = 0 ; d < DIM ; d++ ){
            constants[d] = (p.coord[d] );
        }
        for(int d = 0 ; d < DIM -1 ; d++ ){

        }
    }

private:



}; //class interpolator

} //namespace First Light
#endif /* INTERPOLATE */
